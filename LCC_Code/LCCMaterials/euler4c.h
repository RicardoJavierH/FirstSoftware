#ifndef EULERLAW4CHH
#define EULERLAW4CHH

#include "conslaw.h"


class TEulerLaw4C : public TConservationLaw {

  double fGamma;

 public :
  TEulerLaw4C(int id) : TConservationLaw(id) {
    fGamma=1.4;
  }
  TEulerLaw4C(TEulerLaw4C &law) : TConservationLaw(law) {
    fGamma = law.Gamma();
  }
  ~TEulerLaw4C() {
  }

  virtual int NStateVariables() const { return 4; }

  virtual void FunctionFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &flux);
  virtual void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob);
  virtual void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &);
  void ValJacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &valjacob);
  int ValJacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal) { return 0; }
  STATE MaxEigJacob(TPZVec<STATE> &u,TPZVec<REAL> &normal);
  STATE ValEigJacob(TPZVec<STATE> &u,int order,int dim=1);

  void RoeMatrix(TPZVec<STATE> &U,TPZVec<STATE> &Up1,TPZFMatrix<STATE> &Roe);
  // EigRoe must to have Roe matrix, Roe eigenvalues, Roe matrix R and R inverse
  void EigRoeMatrix(TPZVec<STATE> &U,TPZVec<STATE> &Up1,TPZVec<STATE> &EigRoe);
  void ValRoeMatrix(TPZVec<STATE> &U,TPZVec<STATE> &Up1,TPZFMatrix<STATE> &ValRoe);

  void Print(std::ostream &out= std::cout);
  /** @brief Returns the integrable dimension of the material */
  virtual int Dimension() const { return 1; }
  void FluxGodunov(double *Ui,double *Fluxi);
	
	STATE Pression(TPZVec<STATE> &U);
	
  TPZMaterial *NewMaterial();
  /**Read data of the material from a istream (file data)*/
  void SetData(std::istream &data);

  int VariableIndex(const std::string &name);
  int NSolutionVariables(int index);
  /**returns the solution associated with the var index based on the finite element approximation*/
  void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,
		TPZVec<STATE> &Solout);
  void VariablesName(TPZVec<std::string> &names);

  REAL Gamma() { return fGamma; }

  /**compute the value of the flux function to be used by ZZ error estimator*/
  void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux) { }
};

inline TPZMaterial *TEulerLaw4C::NewMaterial() {
  TPZMaterial *newmat = new TEulerLaw4C(Id());
  return newmat;
}

#endif
