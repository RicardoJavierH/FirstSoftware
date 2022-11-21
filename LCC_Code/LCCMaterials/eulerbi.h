#ifndef EULERBIHH
#define EULERBIHH

#include "conslaw.h"

class TEulerLaw2D : public TConservationLaw {

	double fGamma;
  // Valor final do dominio no eixo X
  double fXend;

   public :
	TEulerLaw2D(int id);
	TEulerLaw2D(TEulerLaw2D &law);
	~TEulerLaw2D() {
	}

	virtual int NStateVariables() const override { return 4; }
  /**Return number of variables to print into the postprocessing file*/
  int AllVariablesPost() override { return 9; }
  /**Return possible number of variables to print into the postprocessing file*/
  int AllVariablesImplemented() override { return 9; }
    void SetGamma(REAL gamma) { fGamma = gamma; }

  virtual void FunctionFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &flux);
  virtual void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob);
  virtual STATE MaxEigJacob(TPZVec<STATE> &u,TPZVec<REAL> &normal);
  virtual STATE ValEigJacob(TPZVec<STATE> &, int, int);
	
  void InverseJacob(TPZFMatrix<STATE> &mat);
  int ValJacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal);
  virtual void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal);

  void RoeMatrix(TPZVec<STATE> &U,TPZVec<STATE> &Up1,TPZFMatrix<STATE> &Roe);
  // EigRoe must to have Roe matrix, Roe eigenvalues, Roe matrix R and R inverse
  void EigRoeMatrix(TPZVec<STATE> &U,TPZVec<STATE> &Up1,TPZVec<STATE> &EigRoe);
  void ValRoeMatrix(TPZVec<STATE> &U,TPZVec<STATE> &Up1,TPZFMatrix<STATE> &ValRoe);

  void Print(std::ostream &out= std::cout);
  /** @brief Returns the integrable dimension of the material */
  virtual int Dimension() const { return 2; }
  void FluxGodunov(double *Ui,double *Fluxi);
	
	STATE Pressure(TPZVec<STATE> &U);
	
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

	/** To compute the tau matrix to diffusive term */
	//void Tau(TPZFMatrix<STATE> &jacinv,TPZFMatrix<STATE> &valjacob,TPZFMatrix<STATE> &tau);
	
  /* 0 -> condicao fronteira periodica, -n depende do valor do ponto x*/
  int IdBC(double *x);

  /**compute the value of the flux function to be used by ZZ error estimator*/
  void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux) { }

    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    
  void ContributeOverInterface(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZVec<STATE> &up1,REAL weight,
     REAL area,int type,TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &phi0,TPZFMatrix<REAL> &phi1,TPZVec<REAL> &normal,
     TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);

  void IncrementDiffusion(TPZVec<REAL> &x,TPZFMatrix<STATE> &,TPZVec<REAL> &sol,TPZFMatrix<STATE> &dsol,REAL weight,
        TPZFMatrix<REAL> &axes,TPZFMatrix<STATE> &phi,TPZFMatrix<STATE> &dphi,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,REAL );

	/** Para calculo do Tau em matrices separadas */
	void MatrixDiff(TPZVec<STATE> &sol,TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &jacinv,TPZFMatrix<STATE>
	&ATauA,TPZFMatrix<STATE> &ATauB,TPZFMatrix<STATE> &BTauA,TPZFMatrix<STATE> &BTauB);
	void InvJacob2d(TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &jacinv);
	void JacobFlux(TPZVec<STATE> &U,TPZFMatrix<STATE> &Ajacob,TPZFMatrix<STATE>
	&Bjacob);

};

inline TPZMaterial *TEulerLaw2D::NewMaterial() {
  TPZMaterial *newmat = new TEulerLaw2D(Id());
  return newmat;
}


#endif
