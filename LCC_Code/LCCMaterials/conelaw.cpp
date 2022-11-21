#include <stdio.h>
#include <fstream>

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzbndcond.h"

#include "conelaw.h"
#include "commonJC.h"

TConeLaw::TConeLaw(int id) : TConservationLaw(id) {
    fName = "Cone Conservation Law ";
}
TConeLaw::TConeLaw(int id,int type) : TConservationLaw(id) {
}
TConeLaw::TConeLaw(int id,int dim,int type) : TConservationLaw(id) {
}

void TConeLaw::FunctionFlux(TPZVec<STATE> &Ui,TPZFMatrix<STATE> &funcao) {
//    funcao.Resize(1,2); funcao.Zero();
  funcao(0,0) = -1.*Point[1]*Ui[0];
  funcao(0,1) = Point[0]*Ui[0];
}

int TConeLaw::VariableIndex(const std::string &name) {
    int val = TPZMaterial::VariableIndex(name);
    if(val<0)
        return 0;
    return val;
}

void TConeLaw::VariablesName(TPZVec<std::string> &names) {
  names[0] = "state_cone";
}

int TConeLaw::NSolutionVariables(int index) {
  if(index < 1) return 1;
  return TPZMaterial::NSolutionVariables(index);
}

/**returns the solution associated with the var index based on the finite element approximation*/
void TConeLaw::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,
                 int var,TPZVec<STATE> &Solout){
  if(!var) {
    Solout.Resize(1);
    Solout[0] = Sol[var];
    return;
  }
  TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

void TConeLaw::JacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &jacob) {
  jacob(0,0) = -1.*Point[1];
  jacob(0,1) = Point[0];
}
void TConeLaw::JacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal) {
  jacob(0,0) = -1.*Point[1]*normal[0] + normal[1]*Point[0];
}

STATE TConeLaw::MaxEigJacob(TPZVec<STATE> &/*U*/,TPZVec<REAL> &normal) {
  return Max(fabs(Point[0]),fabs(Point[1]));
}

STATE TConeLaw::ValEigJacob(TPZVec<STATE> &u,int order,int dim) {
    return (fabs(Point[0])<fabs(Point[1])?fabs(Point[1]):fabs(Point[0]));
}

void TConeLaw::ValJacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &valjacob) {
  valjacob(0,0) = fabs(Point[1]);
  valjacob(1,0) = fabs(Point[0]);
}
int TConeLaw::ValJacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal) {
  valjacob(0,0) = fabs((-1.)*normal[0]*Point[1]+normal[1]*Point[0]);
	return 0;
}

void TConeLaw::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
    TPZFMatrix<REAL>  &phi = dataleft.phi;
	TPZVec<STATE> &u = dataleft.sol[0];
    int phr = phi.Rows();
    int nvar = bc.Val2().Rows();
    short in,jn,cn,kn;
    TPZManVector<STATE> v2(nvar,0.);
    for(in=0;in<nvar;in++)
        v2[in] = bc.Val2()(in,0);
            
    switch (bc.Type()) {
        case 0 :            // Dirichlet condition
            for(in = 0 ; in < phr; in++) {
                for(cn=0;cn<nvar;cn++)
                    ef(nvar*in+cn,0) += (STATE)(gBigNumber* phi(in,0) * weight) * v2[cn]*fDeltaT*fCoef;
                for (jn = 0 ; jn < phr; jn++) {
                    for(cn=0;cn<nvar;cn++)
                        ek(nvar*in+cn,nvar*jn+cn) += gBigNumber * phi(in,0) * phi(jn,0) * weight*fDeltaT*fCoef;
                }
            }
        break;
        case 1 :            // Neumann condition
            for(in = 0 ; in < phi.Rows(); in++) {
                for(jn=0;jn<nvar;jn++)
                    ef(nvar*in+jn,0) -= data.normal[1]*v2[jn] * (STATE)(phi(in,0) * weight*fDeltaT*fCoef);
            }
        break;
        case 2 :        // mixed condition
            for(in = 0 ; in < phi.Rows(); in++) {
                for(cn=0;cn<nvar;cn++)
                    ef(nvar*in+cn, 0) += v2[cn] * (STATE)(phi(in, 0) * weight);
                for (jn = 0 ; jn < phi.Rows(); jn++) {
                    for(cn=0;cn<nvar;cn++)
                        for(kn=0;kn<nvar;kn++)   // Acho que esto no funciona
                            ek(nvar*in+cn,nvar*jn+kn) += bc.Val1()(cn,kn) * (STATE)(phi(in,0) * phi(jn,0) * weight);     // peso de contorno => integral de contorno
                }
            }
            break;
        case 3:
			 STATE velxnorm = (-Point[1])*data.normal[0]+(Point[0])*data.normal[1];
			 v2[0] = velxnorm*u[0];
			 for(in = 0 ; in < phi.Rows(); in++) {
				 for(jn=0;jn<nvar;jn++)
					 ef(nvar*in+jn,0) -= v2[jn] * (STATE)(phi(in,0) * weight*fDeltaT*fCoef);
			 }

            break;
    }
}
/**To evaluate true-error, L2-error and estimate error*/
void TConeLaw::Errors(TPZVec<REAL> &x,TPZVec<STATE> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,TPZVec<STATE> &flux,
     TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {

//  SetPoint(x);
  REAL udif = sol[0] - u_exact[0];

  values[0] = fabs(udif);
  values[1] = udif * udif;
  values[2] = fabs(udif);
}

int TConeLaw::IdBC(double *x) {
  return -1;
}

