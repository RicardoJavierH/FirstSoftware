/*******       File : linlaw.h

Header file for class TLinearLaw1D derived from TConservationLaw to
linear conservation laws.

*******              *******/

#ifndef CONELAWIMPLICITH
#define CONELAWIMPLICITH

#include "conelaw.h"

template<class T>
class TPZVec;

class TConeLawImplicit : public TConeLaw {

 public:
    //Construtores e destrutor
    TConeLawImplicit(int id);
    TConeLawImplicit(int id,int type);
    TConeLawImplicit(int id,int dim,int type);
    ~TConeLawImplicit() { }

    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    /**To compute the contribution over edges of the elements*/
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef) override;

    virtual void FunctionFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &flux) override {
        TConeLaw::FunctionFlux(u,flux);
    }
    virtual void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob) override {
        TConeLaw::JacobFlux(u,jacob);
    }
    virtual void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal) override  {
        TConeLaw::JacobFlux(u,jacob,normal);
    }
    virtual STATE MaxEigJacob(TPZVec<STATE> &u,TPZVec<REAL> &normal) override {
        return TConeLaw::MaxEigJacob(u,normal);
    }
    virtual STATE ValEigJacob(TPZVec<STATE> &u,int order,int dim) override {
        return TConeLaw::ValEigJacob(u,order,dim);
    }

    TPZMaterial *NewMaterial() override;

};

inline TPZMaterial *TConeLawImplicit::NewMaterial() {
  TPZMaterial *newmat = new TConeLawImplicit(Id());
  return newmat;
}

#endif
