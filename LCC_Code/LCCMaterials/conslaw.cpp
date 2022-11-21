/*******       FILE :   conslaw.c

Contains the definitions of the methods for the class TConservationLaw.

*******                               *******/

#include "pzaxestools.h"

#include "conslaw.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzbndcond.h"
#include "commonJC.h"

TConservationLaw::TConservationLaw(int id) : TPZDiscontinuousGalerkin(id), fCoef(1.), fRKOrder(1) {
  fMaxEigen = 0.;
  fExplicit = 0;
}
TConservationLaw::TConservationLaw(TConservationLaw &law) : TPZDiscontinuousGalerkin(), fCoef(1.), fRKOrder(1) {
  fName = law.Name();
  fMaxEigen = 0.;
  fExplicit = 0;
}

/** Para ser usado por el analisis y proveer el orden del metodo numerico
y el coeficiente que debe multiplicar las contribuciones. Coeficientes del esquema numerico en el tiempo */
void TConservationLaw::SetOrderAndCoeficient(int order, REAL Coef) {
    fRKOrder = order;
    fCoef = Coef;
}

/** After conservation law creating is necessary to fill default values */
void TConservationLaw::SetDefaultData() {
	fAlfa = 1.;
	fCurrentTime = 0.;
	fDeltaT = 0.01;
	fCFLDifusion = 0.;
	fMaxEigen = 1.;
	fCoef = 1.;
    fRKOrder = 1;

}

STATE TConservationLaw::MaxValJacobFlux(TPZVec<STATE> &u) {
    TPZFMatrix<STATE> jacobian(NStateVariables(),NStateVariables()*Dimension());
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

void TConservationLaw::SetName(std::string &name) {
  if(name.empty()) {
    if(Dimension()==1) fName = "One-dimensional Conservation Law";
    else if(Dimension()==2) fName = "Two-dimensional Conservation Law";
    else if(Dimension()==3) fName = "Three-dimensional Conservation Law";
    else fName = "Conservation Law undefined";
    return;
  }
  fName = name;
}

void TConservationLaw::Print(std::ostream &out) {
  out << "CONSERVATION LAW" << std::endl << fName << std::endl;
  out << Id() << "\tOrder = " << NStateVariables() << std::endl;
}

void TConservationLaw::Clean() {
  fMaxEigen = 0.;
  fAlfa = fDeltaT = fCurrentTime = 0.;
  fCoef = 0.;
  Point[0] = Point[1] = Point[2] = 0.;
}

/**
* @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
* @param data [in] stores all input data
* @param weight [in] is the weight of the integration rule
* @param ek [out] is the stiffness matrix
* @param ef [out] is the load vector
* @since April 16, 2007
*/
void TConservationLaw::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
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
//    flux.Print(std::cout);
    
    TPZFMatrix<STATE> FluxAxes(nvar,dphr,0.);
//    for(c=0;c<dphr;c++) {
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
                PZError << "TConservationLaw::Contribute dimension error " << dphr << std::endl;
        }
    }
 //   }
/*    FluxAxes.Print(std::cout);
    TPZFMatrix<STATE> FluxDphix(nvar,dphc,0.);
    dphix.Print(std::cout);
    FluxDphix = FluxAxes * dphix;
    FluxDphix.Print(std::cout);
*/
    /** Calculo de la funcion fuente para el montaje de las matrizes. */
    TPZVec<STATE> res(nvar, 0.);
	if (fForcingFunction) {
        fForcingFunction->Execute(x, res);
	}

    STATE Coeff = (STATE)(weight * fCoef * fDeltaT);
    STATE prod;
    for (i = 0; i < phr; i++) {
		/**Calculo de keij, apenas o primeiro valor do bloco*/
		for (j = 0; j < phr; j++) {
			for (c = 0; c < nvar; c++) {
                    ek(i*nvar + c, j*nvar + c) += (weight * phi(j, 0) * phi(i, 0));
			}
		}

		/**Calculo de feic*/
		for (c = 0; c < nvar; c++) {
            prod = 0.;
            for(int l=0;l<dphr;l++)
                prod += dphix(l,i)*FluxAxes(c,l);
            ef(i*nvar + c, 0) += phi(i, 0) * (Coeff * res[c] + weight * (sol[c]/fRKOrder));
            ef(i*nvar + c, 0) += Coeff * prod;
		}
	}
    
    // Incremento do termino de difusion implicita
    if(fExplicit==2) {
		 TPZFMatrix<REAL> jacinvxy(dphr, dphr, 0.);
		 for (i = 0; i < dphr; i++)
			 for (j = 0; j < dphr; j++)
				 for (k = 0; k < dphr; k++)
					 jacinvxy(i, j) += (jacinv(i, k)*axes(k, j));
		 
		 int dsolr = dsol.Rows();
		 REAL coeff = fCoef * fCFLDifusion*weight*fDeltaT;
		 TPZFMatrix<STATE> vjacob(nvar, nvar);
		 TPZFMatrix<STATE> matprod(nvar, nvar);
		 TPZFMatrix<STATE> tau(nvar, nvar, 0.);
		 TPZFMatrix<STATE> temp(nvar, nvar, 0.);
		 TPZVec<REAL> Beta(dphr);
		 
		 /** Computa a matrix tau */
		 if (dphr == 1) tau.Identity();
		 else Tau(jacinvxy, sol, tau);
		 tau *= coeff;
		 
		 for (i = 0; i < phr; i++) {
			 for (k = 0; k < dphr; k++) Beta[k] = dphix(k, i);
			 JacobFlux(sol, vjacob, Beta);
			 Multiply(vjacob, tau, matprod);
			 for (j = 0; j < phr; j++) {
				 for (k = 0; k < dphr; k++) Beta[k] = dphix(k, j);
				 JacobFlux(sol, vjacob, Beta);
				 Multiply(matprod, vjacob, temp);
				 for (k = 0; k < nvar; k++)
					 for (int p = 0; p < nvar; p++)
						 ek(nvar*i + k, nvar*j + p) += temp(k, p);   ///?? Mudanca de sinal 07/03/2003
			 }
		 }

	}
	return;
	{
        if(IsZero(sol[0])) return;
        REAL Coef = weight * fCoef * fDeltaT * fCFLDifusion;
        TPZFMatrix<STATE> JacobianOrig(nvar,nvar*dphr,0.);
        TPZFMatrix<STATE> Jacobian(nvar,nvar*dphr,0.);
        TPZFMatrix<STATE> tau(nvar,nvar,0.);
        TPZFMatrix<STATE> JacobianPhixI(nvar,1,0.);
        TPZFMatrix<STATE> JacobianPhixJ(nvar,1,0.);
        TPZFMatrix<STATE> Beta(nvar,dphr,0.);
    /*    TPZFMatrix<REAL> jacinvxy(dphr, dphr, 0.);
        for (i = 0; i < dphr; i++)
            for (j = 0; j < dphr; j++)
                for (k = 0; k < dphr; k++)
                    jacinvxy(i, j) += (jacinv(i, k)*axes(k, j));
*/
        JacobFlux(sol,JacobianOrig);

        Jacobian = JacobianOrig;
        Jacobian.Print(std::cout);
 /*       for(i=0;i<nvar;i++) {   // Aqui faltam elementos porque el numero de columnas es nvar*dphr
            Jacobian(i,0) = axes(0,0)*JacobianOrig(i,0)+axes(0,1)*JacobianOrig(i,1);
            Jacobian(i,1) = axes(1,0)*JacobianOrig(i,0)+axes(1,1)*JacobianOrig(i,1);
        }*/
        Tau(jacinv,sol,tau);
        tau.Print(std::cout);
        Beta = tau*Jacobian;
        Beta.Print(std::cout);

        for(i=0;i<phr;i++) {
            for(c=0;c<nvar;c++) {
                for(k=0;k<dphr;k++) {
                    JacobianPhixI(c,0) += Beta(c,k)*dphix(k,i);
                }
            }
            for(j=0;j<phr;j++) {
                for(c=0;c<nvar;c++) {
                    for(k=0;k<dphr;k++) {
                        JacobianPhixJ(c,0) += Jacobian(c,k)*dphix(k,j);
                    }
                }
                for(c=0;c<nvar;c++)
                    for(k=0;k<nvar;k++)
                        ek(i*nvar+c,j*nvar+k) -= Coef*JacobianPhixI(c,0)*JacobianPhixJ(k,0);
            }
        }
    }
}

void TConservationLaw::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    TConservationLaw::ContributeInterface(data,dataleft, dataright, weight, ef);
}

void TConservationLaw::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,REAL weight,TPZFMatrix<STATE> &ef) {

    TPZVec<STATE> ul = dataleft.sol[0];
    TPZVec<STATE> ur = dataright.sol[0];
    TPZFMatrix<REAL> axesr = dataright.axes;
    TPZFMatrix<REAL> axesl = dataleft.axes;
    TPZFMatrix<REAL> axes = data.axes;
    TPZFMatrix<REAL> phil = dataleft.phi;
    TPZFMatrix<REAL> phir = dataright.phi;
    TPZManVector<REAL> normal = data.normal;

    SetPoint(data.x);
    
//    if(IsZero(ul[0]) && IsZero(ur[0]))
//    if((fabs(ul[0]) < 1.e-8 ) && (fabs(ur[0]) < 1.e-8 ))
  //      return;

    REAL Coeff = weight*fCoef*fDeltaT;
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
    TPZFMatrix<STATE> FluxAxes(nvar,dim,0.);
    TPZFMatrix<STATE> FluxAxesr(nvar,dim,0.);
    TPZFMatrix<STATE> FluxAxesl(nvar,dim,0.);
    TPZManVector<STATE> DeltaU(nvar,0.);
    TPZFMatrix<STATE> Flux(nvar,dim,0.);
    TPZFMatrix<STATE> FluxNormal(nvar,1,0.);

    FunctionFlux(ul, FluxAxesl);
    FunctionFlux(ur, FluxAxesr);
    fluxr = FluxAxesr;
    fluxl = FluxAxesl;
    /*for(c=0;c<nvar;c++) {
        switch(dim) {
            case 1:
                fluxl(c,0) = axesl(0,0)*FluxAxesl(c,0);
                break;
            case 2:
                fluxl(c,0) = axesl(0,0)*FluxAxesl(c,0)+axesl(0,1)*FluxAxesl(c,1);
                fluxl(c,1) = axesl(1,0)*FluxAxesl(c,0)+axesl(1,1)*FluxAxesl(c,1);
                break;
            case 3:
                fluxl(c,0) = axesl(0,0)*FluxAxesl(c,0)+axesl(0,1)*FluxAxesl(c,1)+axesl(0,2)*FluxAxesl(c,2);
                fluxl(c,1) = axesl(1,0)*FluxAxesl(c,0)+axesl(1,1)*FluxAxesl(c,1)+axesl(1,2)*FluxAxesl(c,2);
                fluxl(c,2) = axesl(2,0)*FluxAxesl(c,0)+axesl(2,1)*FluxAxesl(c,1)+axesl(2,2)*FluxAxesl(c,2);
                break;
            default:
                PZError << "TConeLaw::ContributeInterface dimension error " << dim << endl;
        }
    }
    for(c=0;c<nvar;c++) {
        switch(dim) {
            case 1:
                fluxr(c,0) = axesl(0,0)*FluxAxesr(c,0);
                break;
            case 2:
                fluxr(c,0) = axesr(0,0)*FluxAxesr(c,0)+axesr(0,1)*FluxAxesr(c,1);
                fluxr(c,1) = axesr(1,0)*FluxAxesr(c,0)+axesr(1,1)*FluxAxesr(c,1);
                break;
            case 3:
                fluxr(c,0) = axesl(0,0)*FluxAxesr(c,0)+axesl(0,1)*FluxAxesr(c,1)+axesl(0,2)*FluxAxesr(c,2);
                fluxr(c,1) = axesl(1,0)*FluxAxesr(c,0)+axesl(1,1)*FluxAxesr(c,1)+axesl(1,2)*FluxAxesr(c,2);
                fluxr(c,2) = axesl(2,0)*FluxAxesr(c,0)+axesl(2,1)*FluxAxesr(c,1)+axesl(2,2)*FluxAxesr(c,2);
                break;
            default:
                PZError << "TConeLaw::ContributeInterface dimension error " << dim << endl;
        }
    }*/
    double aux,alfa;
    for(i=0;i<nvar;i++) {
      FluxNormal(i,0) = 0.;
      for(int j=0;j<dim;j++)
        FluxNormal(i,0) += normal[j]*(fluxl(i,j)+fluxr(i,j));
    }

    alfa = MaxEigJacob(ul,normal);
    aux = MaxEigJacob(ur,normal);
    if(alfa<aux) alfa = aux;

    for(i=0;i<nvar;i++) {
      FluxNormal(i,0) += alfa*(ul[i]-ur[i]);
      FluxNormal(i,0) *=.5;
    }

    // Obteniendo phi*F(u_0^k) * det(Jacobiano) * w para la derecha y la izquierda
    for(i=0;i<phrl;i++) {
        for(c=0;c<nvar;c++) {
//            if((FluxNormal(c,0)>0 && ul[c]>0.) || (FluxNormal(c,0)<0 && ur[c]>0.))
                ef(i*nvar+c,0) -= Coeff * FluxNormal(c,0) * phil[i];
        }
    }
    for(i=0;i<phrr;i++) {
        for(c=0;c<nvar;c++) {
//            if((FluxNormal(c,0)>0 && ul[c]>0.) || (FluxNormal(c,0)<0 && ur[c]>0.))
                ef((phrl+i)*nvar+c,0) += Coeff * FluxNormal(c,0) * phir[i];
        }
    }
    // Utilizando un valor de flujo dado por el Flujo Numerico de Lax Friedrichs
    //Flux=.5*((*fL).MaxEigVal()*(u[i]-u[i+1])+(Funcao(u[i])+Funcao(u[i+1])))
}

void TConservationLaw::SetPoint(TPZVec<REAL> &x) {
  int i, dim = x.NElements();
  for(i=0;i<dim;i++) Point[i] = x[i];
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
void TConservationLaw::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
        
        TPZFMatrix<REAL>  &phi = data.phi;
        int phr = phi.Rows();
    int nvar = bc.Val2().Rows();
        short in,jn,cn,kn;
    TPZManVector<STATE> v2(nvar,0.);
    for(in=0;in<nvar;in++)
        v2[in] = bc.Val2()(in,0);
        
        if(bc.HasForcingFunction()) {            // phi(in, 0) = phi_in                          // JORGE 2013 01 26
            TPZManVector<STATE,1> res(1);
            bc.ForcingFunction()->Execute(data.x,res);       // dphi(i,j) = dphi_j/dxi
            for(in=0;in<nvar;in++)
                v2[in] = res[in];
        }

        switch (bc.Type()) {
            case 0 :            // Dirichlet condition
                for(in = 0 ; in < phr; in++) {
                    for(cn=0;cn<nvar;cn++)
                        ef(nvar*in+cn,0) += (STATE)(gBigNumber* phi(in,0) * weight) * v2[cn]*fDeltaT*fCoef;
                    for (jn = 0 ; jn < phr; jn++) {
                        for(cn=0;cn<nvar;cn++)
//                            for(kn=0;kn<nvar;kn++)
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
        }

}

/** To compute the tau matrix to diffusive term */
void TConservationLaw::Tau(TPZFMatrix<REAL> &jacinv,TPZVec<STATE> &sol,TPZFMatrix<STATE> &tau) {
	tau.Zero();
	int i, j, p, dim=jacinv.Rows(), nvar=NStateVariables();
	TPZVec<REAL> Beta(Dimension(),0.);
	TPZFMatrix<STATE> vjacob(nvar,nvar,0.);
	for(i=0;i<dim;i++) {
		for(j=0;j<dim;j++)
			Beta[j] = jacinv(j,i);
		p = ValJacobFlux(sol,vjacob,Beta);
		tau += vjacob;
	}
	if(!p) InverseJacob(tau);
	else tau.Identity();
}

/** To compute the inverse jacobian matrix */
void TConservationLaw::InverseJacob(TPZFMatrix<STATE> &mat) {
	mat(0,0) = 1./mat(0,0);
}

void TConservationLaw::VariablesName(TPZVec<std::string> &names) {
  names[0] = "state";
}

STATE TConservationLaw::MaxValJacobFlux(TPZCompMesh *mesh) {
    STATE maxval = 0.;
    if(!mesh) return maxval;
    int64_t i, nelem = mesh->NElements();
    for(i=0;i<nelem;i++) {
        TPZInterpolationSpace *cel = (TPZInterpolationSpace *)mesh->ElementVec()[i];
        if(!cel) continue;
        int ncorners = cel->Reference()->NCornerNodes();
        STATE val;
        for(int j=0;j<ncorners;j++) {
            TPZManVector<REAL> qsi(cel->Dimension(),0.), point(3,0.);
            TPZVec<REAL> normal;
            TPZMaterialData data;
            data.SetAllRequirements(true);
            cel->Reference()->NodePtr(j)->GetCoordinates(point);
            SetPoint(point);
            cel->Reference()->ComputeXInverse(point,qsi,1.e-8);
            cel->InitMaterialData(data);
            cel->ComputeRequiredData(data,qsi);
            val = fabs(ValEigJacob(data.sol[0],2,mesh->Dimension()));
            if(maxval < val) maxval = val;
        }
    }
    return maxval;
}

// Retorna el valor del paso de tiempo para satisfazer la condicion CFL: Dt < (Dx / (2p+1)*|maxF`(u)|)
REAL TConservationLaw::CFLCondition(TPZCompMesh *mesh) {
    REAL MaxJacobFlux = MaxValJacobFlux(mesh);
    if(MaxJacobFlux>10.) MaxJacobFlux = 10;
    if(IsZero(MaxJacobFlux)) return 0.0;
    REAL h = MinimumDeltaX(mesh);
    return (fCFL*h)/MaxJacobFlux;
}
