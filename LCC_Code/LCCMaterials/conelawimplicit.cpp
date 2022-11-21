#include "stdio.h"
#include "conelawimplicit.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include <fstream>
#include "commonJC.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"

TConeLawImplicit::TConeLawImplicit(int id) : TConeLaw(id) {
    fName = "Cone Conservation Law ";
}
TConeLawImplicit::TConeLawImplicit(int id,int type) : TConeLaw(id) {
}
TConeLawImplicit::TConeLawImplicit(int id,int dim,int type) : TConeLaw(id) {
}


void TConeLawImplicit::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    TPZVec<REAL> x = data.x;
    TPZFMatrix<REAL> jacinv = data.jacinv;
    TPZVec<STATE> sol = data.sol[0];
    TPZFMatrix<STATE> dsol = data.dsol[0];
    TPZFMatrix<REAL> axes = data.axes;
    TPZFMatrix<REAL> phi = data.phi;
    TPZFMatrix<REAL> dphix = data.dphix;
    int phc, phr, dphc, dphr, efr, efc, ekr, ekc;
    int nvar = NStateVariables();
    phc = phi.Cols();
    phr = phi.Rows();
    dphc = dphix.Cols();
    dphr = dphix.Rows();   // It is dimension.
    if (!phr || dphr != Dimension()) return;

    SetPoint(x);

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

    int i,j,c;
    TPZVec<STATE> res(nvar, 0.);
    TPZFMatrix<STATE> JacobianDphix(1,phr,0.);
    for(i=0;i<phr;i++)
        JacobianDphix(0,i) = -1.*Point[1]*dphix(0,i)+Point[0]*dphix(1,i);

    if (fForcingFunction) {
        /**Observar que res contem a solucao Ul-1, e a funcao s(tk,x,sol(l-1))
           toma a tempo tk da lei de conservacao*/
        res.Fill(0.);
        fForcingFunction->Execute(x, res);
    }

    for (i = 0; i < phr; i++) {
        /**Calculo de keij, apenas o primeiro valor do bloco*/
        for (j = 0; j < phr; j++) {
            for (c = 0; c < nvar; c++) {
                ek(i*nvar + c, j*nvar + c) += (weight*phi(j, 0)*phi(i, 0));
                ek(i*nvar + c, j*nvar + c) -= fDeltaT* fCoef * weight* JacobianDphix(0,i)*phi(j, 0);
            }
        }

        /**Calculo de feic*/
        for (c = 0; c < nvar; c++) {
            ef(i*nvar + c, 0) += weight * phi(i, 0) * (fDeltaT*fCoef*res[c] + (sol[c]/fRKOrder));
        }
    }
}

void TConeLawImplicit::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef) {
    return;
}
void TConeLawImplicit::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {

    TPZVec<STATE> ul = dataleft.sol[0];
    TPZVec<STATE> ur = dataright.sol[0];
    TPZFMatrix<STATE> dsol = data.dsol[0];
    TPZFMatrix<REAL> phil = dataleft.phi;
    TPZFMatrix<REAL> phir = dataright.phi;
    TPZVec<REAL> normal = data.normal;

    SetPoint(data.x);

    int phrl,phrr,efr;
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
    int i, j;
    // Creando vectores F(u^k+1)*Normal
    TPZFMatrix<STATE> FluxNormall(1,phrl,0.);
    TPZFMatrix<STATE> FluxNormalr(1,phrr,0.);
    REAL JacobianNormal = -1.*Point[1]*normal[0]+Point[0]*normal[1];
    for(i=0;i<phrl;i++)
        FluxNormall(0,i) = JacobianNormal*phil(i,0);
    for(i=0;i<phrr;i++)
        FluxNormalr(0,i) = JacobianNormal*phir(i,0);
    REAL Theta = 1.;

    // 1) phi_I_left, phi_J_left
    for(i=0; i<phrl; i++) {
        for(j=0; j<phrl; j++) {
            ek(i,j) += (STATE)(weight * fDeltaT * fCoef * ((-Theta) * FluxNormall(0,i) * phil(j,0) + (Theta) * FluxNormall(0,j) * phil(i,0)));
        }
    }
    // 2) phi_I_right, phi_J_right
    for(i=0; i<phrr; i++) {
        for(j=0; j<phrr; j++) {
            ek(i,j) += (STATE)(weight * fDeltaT * fCoef * ((Theta) * FluxNormalr(0,i) * phir(j,0) + (-Theta) * FluxNormalr(0,j) * phir(i,0)));
        }
    }
    // 3) phi_I_left, phi_J_right
    for(i=0; i<phrl; i++) {
        for(j=0; j<phrr; j++) {
            ek(i,j) += (STATE)(weight * fDeltaT * fCoef * ((-Theta) * FluxNormall(0,i) * phir(j,0) + (-Theta) * FluxNormalr(0,j) * phil(i,0)));
        }
    }
    // 4) phi_I_right, phi_J_left
    for(i=0; i<phrr; i++) {
        for(j=0; j<phrl; j++) {
            ek(i,j) += (STATE)(weight * fDeltaT * fCoef * ((Theta) * FluxNormalr(0,i) * phil(j,0) + (Theta) * FluxNormall(0,j) * phir(i,0)));
        }
    }
}
