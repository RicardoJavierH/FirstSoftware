#include "stdio.h"
#include "conelawSI.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include <fstream>
#include "commonJC.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"

TConeLawSI::TConeLawSI(int id) : TConservationLawSI(id) {
    fName = "Cone Conservation Law ";
}
TConeLawSI::TConeLawSI(int id,int type) : TConservationLawSI(id) {
}
TConeLawSI::TConeLawSI(int id,int dim,int type) : TConservationLawSI(id) {
}

void TConeLawSI::FunctionFlux(TPZVec<STATE> &Ui,TPZFMatrix<STATE> &funcao) {
    funcao.Resize(1,2); funcao.Zero();
  funcao(0,0) = -1.*Point[1]*Ui[0];
  funcao(0,1) = Point[0]*Ui[0];
}

int TConeLawSI::VariableIndex(const std::string &name) {
  return 0;
}

void TConeLawSI::VariablesName(TPZVec<std::string> &names) {
  names[0] = "state_cone";
}

int TConeLawSI::NSolutionVariables(int index) {
  if(index < 1) return 1;
  return TPZMaterial::NSolutionVariables(index);
}

/**returns the solution associated with the var index based on the finite element approximation*/
void TConeLawSI::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,
                 int var,TPZVec<STATE> &Solout){
  if(!var) {
    Solout.Resize(1);
    Solout[0] = Sol[var];
    return;
  }
  TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

void TConeLawSI::JacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &jacob) {
  jacob(0,0) = -1.*Point[1];
  jacob(0,1) = Point[0];
}
void TConeLawSI::JacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal) {
  jacob(0,0) = -1.*Point[1]*normal[0] + normal[1]*Point[0];
}

STATE TConeLawSI::MaxEigJacob(TPZVec<STATE> &/*U*/,TPZVec<REAL> &normal) {
  return Max(fabs(Point[0]),fabs(Point[1]));
}

STATE TConeLawSI::ValEigJacob(TPZVec<STATE> &u,int order,int dim) {
#ifndef NOTDEBUG
  if(order) PZError << "TConeLawSI::ValEigJacob, undefined order.\n";
#endif
  if(dim==1) 
    return -1.*Point[1];
  return Point[0];
}

void TConeLawSI::ValJacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &valjacob) {
  valjacob(0,0) = fabs(Point[1]);
  valjacob(1,0) = fabs(Point[0]);
}
int TConeLawSI::ValJacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal) {
  valjacob(0,0) = fabs((-1.)*normal[0]*Point[1]+normal[1]*Point[0]);
	return 0;
}

/**To evaluate true-error, L2-error and estimate error*/
void TConeLawSI::Errors(TPZVec<REAL> &x,TPZVec<STATE> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,TPZVec<STATE> &flux,
     TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {

//  SetPoint(x);
  REAL udif = sol[0] - u_exact[0];

  values[0] = fabs(udif);
  values[1] = udif * udif;
  values[2] = fabs(udif);
}

int TConeLawSI::IdBC(double *x) {
  return -1;
}

