
#ifndef EULERLAW1DH
#define EULERLAW1DH

#include "conslaw.h"


class TEulerLaw1D : public TConservationLaw {

	REAL fGamma;

   public :
	TEulerLaw1D(int id);
	TEulerLaw1D(TEulerLaw1D &law);
	~TEulerLaw1D() {
	}

	virtual  int NStateVariables() const { return 3; }

  virtual void FunctionFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &flux);
  virtual void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob);
	virtual void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal);

  int ValJacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &Beta);
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
	
	STATE Pressure(TPZVec<STATE> &U);
	
	double Phi(double w);
	void SolRiemannProblem(double *Ul,double *Ur,double *Result);

  TPZMaterial *NewMaterial();
  /**Read data of the material from a istream (file data)*/
  void SetData(std::istream &data);

  int VariableIndex(const std::string &name);
  int NSolutionVariables(int index);
  /**returns the solution associated with the var index based on the finite element approximation*/
  void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,
		TPZVec<STATE> &Solout);
  void VariablesName(TPZVec<std::string> &names);

	void Saida();
  REAL Gamma() { return fGamma; }

  /**compute the value of the flux function to be used by ZZ error estimator*/
  void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux) { }

  /** 0 -> condicao fronteira periodica, -n depende do valor do ponto x*/
  int IdBC(double *x);
  /** To compute the contribution over edges of the elements*/
  void ContributeOverInterface(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZVec<STATE> &up1,
    REAL weight,REAL area,int type,TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &phileft,TPZFMatrix<REAL> &phiright,TPZVec<REAL> &normal,
    TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);

  /** Return number of variables to print in postprocessing file*/
  int AllVariablesPost() { return 1; }
  /** Return possible number of variables to print into the postprocessing file*/
  int AllVariablesImplemented() { return 7; }

  void IncrementDiffusion(TPZVec<REAL> &,TPZFMatrix<STATE> &,TPZVec<REAL> &,TPZFMatrix<STATE> &,REAL ,
        TPZFMatrix<STATE> &,TPZFMatrix<STATE> &,TPZFMatrix<STATE> &,TPZFMatrix<STATE> &,TPZFMatrix<STATE> &,REAL );
	void InverseJacob(TPZFMatrix<STATE> &jac);
};

inline TPZMaterial *TEulerLaw1D::NewMaterial() {
  TPZMaterial *newmat = new TEulerLaw1D(Id());
  return newmat;
}

/*class TEulerLaw2D1D : public TConservationLaw1D {

	double vetorX,vetorY;
	double fGamma;
   public :
	TEulerLaw2D1D(FILE *type) : TConservationLaw1D(type) {
		printf("\nVetor direcao : ");
		scanf("%lf%lf",&vetorX,&vetorY);
		fGamma=1.4;
		fOrdem=4;
	}
	TEulerLaw2D1D(char *filein) : TConservationLaw1D(filein) {
		printf("\nVetor direcao : ");
		scanf("%lf%lf",&vetorX,&vetorY);
		fGamma=1.4;
		fOrdem=4;
	}
	TEulerLaw2D1D(double *U,double *U1,char *Name,int NCel) :
		TConservationLaw1D(U,U1,Name,4,NCel) {
		printf("\nVetor direcao : ");
		scanf("%lf%lf",&vetorX,&vetorY);
		fGamma=1.4;
		fOrdem=4;
	}
	~TEulerLaw2D1D() { }

	void FuncaoFlux(double *Ui,double *funcao);
	void JacobFlux(double *Ui,double *jacob);
	void GeraRoeMatrix(double *Ul,double *Ur,double *Roe);

	double MaxEigValJacob(double *Ui);\

	void Saida() { }
};
*/

#endif

