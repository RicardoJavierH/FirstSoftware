/*******       FILE :   conslaw.c

Contains the definitions of the methods for the class TConservationLaw.

*******                               *******/

#include "conslawSI.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzbndcond.h"
#include "commonJC.h"

TConservationLawSI::TConservationLawSI(int id) : TConservationLaw(id) {
  fMaxEigen = 0.;
  fExplicit = 0;
}
TConservationLawSI::TConservationLawSI(TConservationLawSI &law) : TConservationLaw(law) {
  fMaxEigen = 0.;
  fExplicit = 0;
}

STATE TConservationLawSI::MaxValJacobFlux(TPZVec<STATE> &u) {
    TPZFMatrix<STATE> jacobian(NStateVariables(),Dimension());
    JacobFlux(u,jacobian);
    int i,j;
    REAL maxval = 0.;
    for(i=0;i<jacobian.Rows();i++) {
        for(j=0;j<jacobian.Cols();j++) {
            maxval += jacobian(i,j)*jacobian(i,j);
        }
    }
    return sqrt(maxval);
}

/**
* @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
* @param data [in] stores all input data
* @param weight [in] is the weight of the integration rule
* @param ek [out] is the stiffness matrix
* @param ef [out] is the load vector
* @since April 16, 2007
*/
void TConservationLawSI::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
	TPZVec<REAL> x = data.x; 
	TPZVec<STATE> sol = data.sol[0];
	TPZFMatrix<STATE> dsol = data.dsol[0];
	TPZFMatrix<REAL> axes = data.axes;
    TPZFMatrix<REAL> jacinv = data.jacinv;
	TPZFMatrix<REAL> phi = data.phi;
	TPZFMatrix<REAL> dphix = data.dphix;
	int phc, phr, dphc, dphr, efr, efc, ekr, ekc;
	int nvar = NStateVariables();
	phc = phi.Cols();
	phr = phi.Rows();
	dphc = dphix.Cols();
    dphr = dphix.Rows();
	if (!phr || dphr != Dimension()) return;
    efr = ef.Rows();
    efc = ef.Cols();
    ekr = ek.Rows();
    ekc = ek.Cols();
    if (phc != 1 || phr != dphc ||
        ekr != phr * nvar || ekc != ekr ||
        efr != ekr || efc != 1) {
        PZError << "\nTConservationLaw. Inconsistent input data : \n" <<
            "phi.Cols() = " << phi.Cols() << " dphi.Cols() = " << dphix.Cols() <<
            " phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
            dphix.Rows() << "\nek.Rows() = " << ek.Rows() << " ek.Cols() = "
            << ek.Cols() <<
            "\nef.Rows() = " << ef.Rows() << " ef.Cols() = "
            << ef.Cols() << "\n";
    }

    int i, j, k, c;
    // Almacena el punto de integracion para uso de los flujos numericos — otras funciones
	SetPoint(x);

    /* Ajuste del vector convectivo para los ejes del elemento geomŽtrico. */
    /** Calculo de la funcion flujo  */
    TPZFMatrix<STATE> flux(nvar,dphr, 0.); //dphr is a spatial domain dimension
    FunctionFlux(sol, flux);
    
    TPZFMatrix<STATE> FluxAxes(nvar,dphr,0.);
    for(i=0;i<nvar;i++) {
        switch(dphr) {
            case 1:
                FluxAxes(i,0) = axes(0,0)*flux(i,0);
                break;
            case 2:
                FluxAxes(i,0) = axes(0,0)*flux(i,0)+axes(0,1)*flux(i,1);
                FluxAxes(i,1) = axes(1,0)*flux(i,0)+axes(1,1)*flux(i,1);
                break;
            case 3:
                FluxAxes(i,0) = axes(0,0)*flux(i,0)+axes(0,1)*flux(i,1)+axes(0,2)*flux(i,2);
                FluxAxes(i,1) = axes(1,0)*flux(i,0)+axes(1,1)*flux(i,1)+axes(1,2)*flux(i,2);
                FluxAxes(i,2) = axes(2,0)*flux(i,0)+axes(2,1)*flux(i,1)+axes(2,2)*flux(i,2);
                break;
            default:
                PZError << "TPZMatPoisson3d::Contribute dimension error " << dphr << std::endl;
        }
    }

    TPZFMatrix<STATE> FluxDphix(nvar,dphc,0.);
    for(i=0;i<nvar;i++)
        for(j=0;j<nvar;j++)
            for(k=0;k<dphr;k++)
                FluxDphix(i,j) = FluxAxes(i,k)*dphix(k,j);
    
    /** Calculo de la funcion fuente para el montaje de las matrizes. */
    TPZVec<STATE> res(nvar, 0.);
	if (fForcingFunction) {
        fForcingFunction->Execute(x, res);
	}

    STATE Coef = (STATE)(weight * fCoef * fDeltaT);
    for (i = 0; i < phr; i++) {
		/**Calculo de keij, apenas o primeiro valor do bloco*/
		for (j = 0; j < phr; j++) {
			for (c = 0; c < nvar; c++) {
				ek(i*nvar + c, j*nvar + c) += (weight * (phi(j, 0)-fDeltaT*FluxDphix(c,j)) * phi(i, 0));
			}
		}

		/**Calculo de feic*/
		for (c = 0; c < nvar; c++) {
            ef(i*nvar + c, 0) += phi(i, 0) * (Coef * res[c] + weight * (sol[c]/fRKOrder));
           // ef(i*nvar + c, 0) += Coef * FluxDphix(c,i);
		}
	}
    
    // Incremento do termino de difusion implicita
    if(fExplicit==2) {
        REAL Coef = weight * fCoef * fDeltaT * fCFLDifusion;
        TPZFMatrix<STATE> JacobianOrig(nvar,dphr,0.);
        TPZFMatrix<STATE> Jacobian(nvar,dphr,0.);
        TPZFMatrix<STATE> JacobianPhixI(nvar,1,0.);
        TPZFMatrix<STATE> JacobianPhixJ(nvar,1,0.);
        TPZFMatrix<STATE> JacobianInvJacobE(nvar,1,0.);
        TPZFMatrix<STATE> JacobianInvJacobN(nvar,1,0.);
        JacobFlux(sol,JacobianOrig);
        TPZFMatrix<STATE> jacinvxy(dphr, dphr, 0.);
//        jacinvxy = jacinv;
        Jacobian = JacobianOrig;
//        for(i=0;i<nvar;i++) {
//            Jacobian(i,0) = axes(0,0)*JacobianOrig(i,0)+axes(0,1)*JacobianOrig(i,1);
//            Jacobian(i,1) = axes(1,0)*JacobianOrig(i,0)+axes(1,1)*JacobianOrig(i,1);
//        }
        //Computando o jacobiano da transformacao do elemento mestre ao elemento triangular
        for (i = 0; i < dphr; i++)
            for (j = 0; j < dphr; j++)
                for (k = 0; k < dphr; k++)
                    jacinvxy(i, j) += (jacinv(i, k)*axes(k, j));

        REAL temp;
        for(c=0;c<nvar;c++) {
            for(k=0;k<dphr;k++) {
                JacobianInvJacobE(c,0) += Jacobian(c,k)*jacinvxy(0,k);
                temp = fabs(JacobianInvJacobE(c,0));
                JacobianInvJacobE(c,0) = abs(temp);
            }
        }
        for(c=0;c<nvar;c++) {
            for(k=0;k<dphr;k++) {
                JacobianInvJacobN(c,0) += Jacobian(c,k)*jacinvxy(1,k);
                temp = fabs(JacobianInvJacobN(c,0));
                JacobianInvJacobN(c,0) = temp;
            }
        }
        JacobianInvJacobE += JacobianInvJacobN;
        for(c=0;c<nvar;c++)
            JacobianInvJacobN(c,0) = 1./JacobianInvJacobE(c,0);
        
        for(i=0;i<phc;i++) {
            for(c=0;c<nvar;c++) {
                for(k=0;k<dphr;k++) {
                    JacobianPhixI(c,0) += Jacobian(c,k)*dphix(k,i);
                }
            }
            for(j=0;j<phr;j++) {
                for(c=0;c<nvar;c++) {
                    for(k=0;k<dphr;k++) {
                        JacobianPhixJ(c,0) += Jacobian(c,k)*dphix(k,j);
                    }
                }
                for(c=0;c<nvar;c++)
                    ek(i,j) -= Coef*JacobianPhixI(c,0)*JacobianPhixJ(c,0)*JacobianInvJacobN(c,0);
            }
        }
    }
}

void TConservationLawSI::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    TConservationLawSI::ContributeInterface(data,dataleft, dataright, weight, ef);
}

void TConservationLawSI::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,REAL weight,TPZFMatrix<STATE> &ef) {

    TPZVec<STATE> ul = dataleft.sol[0];
    TPZVec<STATE> ur = dataright.sol[0];
    TPZFMatrix<REAL> axesr = dataright.axes;
    TPZFMatrix<REAL> axesl = dataleft.axes;
    TPZFMatrix<REAL> phil = dataleft.phi;
    TPZFMatrix<REAL> phir = dataright.phi;
    TPZVec<REAL> normal = data.normal;

    SetPoint(data.x);

    int phrl,phrr,efr;
    int dim = Dimension();
    int nvar = NStateVariables();
    phrl = phil.Rows();
    phrr = phir.Rows();
    efr = ef.Rows();
#ifndef NOTDEBUG
    if(phrl != phrr || efr != (phrl+phrr)*nvar) {
        PZError << "\nTConservationLaw. Inconsistent input data : \n" <<
        " phi.Rows = " << (phil.Rows()+phir.Rows()) << "\nef.Rows() = " << ef.Rows() << "\n";
    }
#endif
    int i, c;
    TPZFMatrix<STATE> fluxr(nvar,dim,0.);
    TPZFMatrix<STATE> fluxl(nvar,dim,0.);
    FunctionFlux(ur, fluxr);
    FunctionFlux(ul, fluxl);
    TPZFMatrix<STATE> FluxAxesl(nvar,dim,0.);
    FluxAxesl = fluxl;
/*    for(i=0;i<nvar;i++) {
        switch(dim) {
            case 1:
                FluxAxesl(i,0) = axesl(0,0)*fluxl(i,0);
                break;
            case 2:
                FluxAxesl(i,0) = axesl(0,0)*fluxl(i,0)+axesl(0,1)*fluxl(i,1);
                FluxAxesl(i,1) = axesl(1,0)*fluxl(i,0)+axesl(1,1)*fluxl(i,1);
                break;
            case 3:
                FluxAxesl(i,0) = axesl(0,0)*fluxl(i,0)+axesl(0,1)*fluxl(i,1)+axesl(0,2)*fluxl(i,2);
                FluxAxesl(i,1) = axesl(1,0)*fluxl(i,0)+axesl(1,1)*fluxl(i,1)+axesl(1,2)*fluxl(i,2);
                FluxAxesl(i,2) = axesl(2,0)*fluxl(i,0)+axesl(2,1)*fluxl(i,1)+axesl(2,2)*fluxl(i,2);
                break;
            default:
                PZError << "TConservationLaw::ContributeInterface dimension error (l) " << dim << endl;
        }
    }*/
    TPZFMatrix<STATE> FluxAxesr(nvar,dim,0.);
    FluxAxesr = fluxr;
/*    for(i=0;i<nvar;i++) {
        switch(dim) {
            case 1:
                FluxAxesr(i,0) = axesr(0,0)*fluxr(i,0);
                break;
            case 2:
                FluxAxesr(i,0) = axesr(0,0)*fluxr(i,0)+axesr(0,1)*fluxr(i,1);
                FluxAxesr(i,1) = axesr(1,0)*fluxr(i,0)+axesr(1,1)*fluxr(i,1);
                break;
            case 3:
                FluxAxesr(i,0) = axesr(0,0)*fluxr(i,0)+axesr(0,1)*fluxr(i,1)+axesr(0,2)*fluxr(i,2);
                FluxAxesr(i,1) = axesr(1,0)*fluxr(i,0)+axesr(1,1)*fluxr(i,1)+axesr(1,2)*fluxr(i,2);
                FluxAxesr(i,2) = axesr(2,0)*fluxr(i,0)+axesr(2,1)*fluxr(i,1)+axesr(2,2)*fluxr(i,2);
                break;
            default:
                PZError << "TConservationLaw::ContributeInterface dimension error (r) " << dim << endl;
        }
    }*/
    // Obteniendo la media de los fluxos laterales para considerar el flujo no ponto de integracion en la interface
//    FluxAxesl += FluxAxesr;
//    FluxAxesl *= 0.5;
    TPZFMatrix<STATE> FluxNormall(nvar,1,0.);
    TPZFMatrix<STATE> FluxNormalr(nvar,1,0.);
    for(i=0;i<nvar;i++)
        for(c=0;c<dim;c++)
            FluxNormalr(i,0) += FluxAxesr(i,c)*normal[c];
    for(i=0;i<nvar;i++)
        for(c=0;c<dim;c++)
            FluxNormall(i,0) += FluxAxesl(i,c)*normal[c];
    FluxNormall += FluxNormalr;
    REAL aux = MaxEigJacob(ul,normal);
    REAL alfa = MaxEigJacob(ur,normal);
    if(alfa<aux) alfa = aux;

    for(i=0;i<nvar;i++) {
      FluxNormall(i,0) += alfa*(ul[i]-ur[i]);
      FluxNormall *= .5;
    }

    // Obteniendo phi*F(u_0^k) * det(Jacobiano) * w para la derecha y la izquierda
    for(i=0;i<phrl;i++) {
        for(c=0;c<nvar;c++) {
     //       if(FluxNormall(c,0)>0.)
                ef(i*nvar+c,0) -=  weight * fDeltaT * fCoef * FluxNormall(c,0) * phil[i];
       //     else
         //       ef(i*nvar+c,0) +=  weight * fDeltaT * fCoef * FluxNormall(c,0) * phil[i];
        }
    }
    for(i=0;i<phrr;i++) {
        for(c=0;c<nvar;c++) {
//            if(FluxNormalr(c,0)>0.)
                ef(phrl+i*nvar+c,0) += weight * fDeltaT * fCoef * FluxNormall(c,0) * phir[i];
  //          else
    //            ef(phrl+i*nvar+c,0) -= weight * fDeltaT * fCoef * FluxNormall(c,0) * phir[i];
        }
    }
    // Utilizando un valor de flujo dado por el Flujo Numerico de Lax Friedrichs
    //Flux=.5*((*fL).MaxEigVal()*(u[i]-u[i+1])+(Funcao(u[i])+Funcao(u[i+1])))
}

/**
 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
 * @param data [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @param bc [in] is the boundary condition material
 * @since October 07, 2011
 */
void TConservationLawSI::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
        
        TPZFMatrix<REAL>  &phi = data.phi;
        int phr = phi.Rows();
        short in,jn;
        STATE v2[1];
        v2[0] = bc.Val2()(0,0);
        
        if(bc.HasForcingFunction()) {            // phi(in, 0) = phi_in                          // JORGE 2013 01 26
            TPZManVector<STATE,1> res(1);
            bc.ForcingFunction()->Execute(data.x,res);       // dphi(i,j) = dphi_j/dxi
            v2[0] = res[0];
        }

        switch (bc.Type()) {
            case 0 :            // Dirichlet condition
                for(in = 0 ; in < phr; in++) {
                    ef(in,0) += (STATE)(gBigNumber* phi(in,0) * weight) * v2[0];
                    for (jn = 0 ; jn < phr; jn++) {
                        ek(in,jn) += gBigNumber * phi(in,0) * phi(jn,0) * weight;
                    }
                }
                break;
            case 1 :            // Neumann condition
                for(in = 0 ; in < phi.Rows(); in++) {
                    ef(in,0) += v2[0] * (STATE)(phi(in,0) * weight);
                }
                break;
            case 2 :        // mixed condition
                for(in = 0 ; in < phi.Rows(); in++) {
                    ef(in, 0) += v2[0] * (STATE)(phi(in, 0) * weight);
                    for (jn = 0 ; jn < phi.Rows(); jn++) {
                        ek(in,jn) += bc.Val1()(0,0) * (STATE)(phi(in,0) * phi(jn,0) * weight);     // peso de contorno => integral de contorno
                    }
                }
                break;
        }

}

STATE TConservationLawSI::MaxValJacobFlux(TPZCompMesh *mesh) {
    STATE maxval = 0.;
    if(!mesh) return maxval;
    int64_t i, nelem = mesh->NElements();
    for(i=0;i<nelem;i++) {
        TPZInterpolationSpace *cel = (TPZInterpolationSpace *)mesh->ElementVec()[i];
        if(!cel) continue;
        int ncorners = cel->Reference()->NCornerNodes();
        STATE val;
        for(int j=0;j<ncorners;j++) {
            TPZManVector<REAL> qsi(cel->Dimension(),0.), point(mesh->Dimension(),0.);
            TPZVec<REAL> normal;
            TPZMaterialData data;
            data.SetAllRequirements(true);
            cel->Reference()->NodePtr(j)->GetCoordinates(point);
            SetPoint(point);
            cel->Reference()->ComputeXInverse(point,qsi,1.e-8);
            cel->InitMaterialData(data);
            cel->ComputeRequiredData(data,qsi);
            val = MaxValJacobFlux(data.sol[0]);
            if(maxval < val) maxval = val;
        }
    }
    return maxval;
}
