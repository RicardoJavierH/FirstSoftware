/*******       File : linlaw.h

Header file for class TLinearLaw1D derived from TConservationLaw to
linear conservation laws.

*******              *******/

#ifndef CONELAWH
#define CONELAWH

#include "conslaw.h"

template<class T>
class TPZVec;

class TConeLaw : public TConservationLaw {

 public:
	 virtual int NStateVariables() const { return 1; }

//    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);

  void FunctionFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &flux);
  void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob);
  void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &);
  void ValJacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &valjacob);
  int ValJacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &Beta);

  //Construtores e destrutor
  TConeLaw(int id);
  TConeLaw(int id,int type);
  TConeLaw(int id,int dim,int type);
  ~TConeLaw() { }

  STATE MaxEigJacob(TPZVec<STATE> &u,TPZVec<REAL> &normal);
  STATE ValEigJacob(TPZVec<STATE> &u,int order,int dim=2);

  /** @brief Returns the integrable dimension of the material */
  virtual int Dimension() const { return 2; }

  TPZMaterial *NewMaterial();

  /**compute the value of the flux function to be used by ZZ error estimator*/
  void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux) { }

  virtual int VariableIndex(const std::string &name);
  /**To return the variables name of the material*/
  virtual void VariablesName(TPZVec<std::string> &names);
  int NSolutionVariables(int index);
  /**returns the solution associated with the var index based on the finite element approximation*/
  void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,
		TPZVec<STATE> &Solout);
  void IncrementDiffusion(TPZVec<REAL> &,TPZFMatrix<STATE> &,TPZVec<REAL> &,TPZFMatrix<STATE> &,REAL ,
        REAL ,TPZFMatrix<STATE> &,TPZFMatrix<STATE> &,TPZFMatrix<STATE> &,TPZFMatrix<STATE> &,TPZFMatrix<STATE> &,REAL );

  /**To evaluate true-error, L2-error and estimate error*/
    void Errors(TPZVec<REAL> &x,TPZVec<STATE> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,TPZVec<STATE> &flux,
        TPZVec<STATE> &uexact,TPZFMatrix<STATE> &duexact,TPZVec<REAL> &val);
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;

  int IdBC(double *x);
};

inline TPZMaterial *TConeLaw::NewMaterial() {
  TPZMaterial *newmat = new TConeLaw(Id());
  return newmat;
}

#endif
