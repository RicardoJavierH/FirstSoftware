/*******       File : linlaw.h

Header file for class TLinearLaw1D derived from TConservationLaw to
linear conservation laws.

*******              *******/

#ifndef LINLAWHHH
#define LINLAWHHH

#include "conslaw.h"

#define MAXORDER 10
template<class T> 
class TPZFMatrix;
template<class T>
class TPZVec;

class TLinearLaw : public TConservationLaw {

 protected:
  int fOrder;                   // Order of the equations system

  TPZFMatrix<STATE> fA;   // A, B and/or C -> convection matrix
  TPZFMatrix<STATE> fValA;// |A| -> |A|=R|EigVal|(InvR)
  REAL fMaxEigVal;              // Valor do maximo autovalor da matrix A
  TPZVec<STATE> fEigVal;         // Vetor de autovalores
  TPZFMatrix<STATE> fEigVect;        // R -> matrix de autovetores
  TPZFMatrix<STATE> fInvEigVect;     // Inversa de R

  TPZFMatrix<STATE> fInverseQrt;
  
 public:
	 virtual int NStateVariables() const { return fOrder; }

  TPZFMatrix<STATE> &A() { return fA; }
  TPZFMatrix<STATE> &ValA() { return fValA; }
  TPZVec<STATE> &EigVal(){ return fEigVal; }
  TPZFMatrix<STATE> &EigVect() { return fEigVect; }
  TPZFMatrix<STATE> &InvEigVect() { return fInvEigVect; }

  virtual void FunctionFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &flux);
//  void Flux(TPZVec<REAL> &u,TPZVec<REAL> &flux);
  void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob);
  void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &);
  void ValJacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &valjacob);
  int ValJacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &Beta);
	/** To compute the inverse jacobian matrix */
	void InverseJacob(TPZFMatrix<STATE> &mat);

  //Construtores e destrutor
  TLinearLaw(int id,int order);
  TLinearLaw(int id,int order,int type);
  ~TLinearLaw() { }

  virtual void SetMatrix(TPZFMatrix<STATE> &A,TPZFMatrix<STATE> &EigVect,TPZFMatrix<STATE> &InvEigVect,TPZVec<STATE> &EigVal);
  virtual void SetMatrix(std::istream &input);
  virtual void RequireMatrix();
  virtual void Print(std::ostream &out= std::cout);
  STATE MaxEigJacob(TPZVec<STATE> &u,TPZVec<REAL> &normal);
  STATE ValEigJacob(TPZVec<STATE> &u,int order,int dim=1);

  /** @brief Returns the integrable dimension of the material */
  virtual int Dimension() const { return 1; }

  TPZMaterial *NewMaterial();
  /**Read data of the material from a istream (file data)*/
  void SetData(std::istream &data);

  /**compute the value of the flux function to be used by ZZ error estimator*/
  void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux) { }

  virtual int VariableIndex(const std::string &name);
  int NSolutionVariables(int index);
  /**returns the solution associated with the var index based on the finite element approximation*/
  void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,
		TPZVec<STATE> &Solout);
  void VariablesName(TPZVec<std::string> &names);
  /**Compute contribution to the stiffness matrix and right hand side at an integration point */
//  void IncrementDiffusion(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix &dsol,
//        REAL weight,REAL area,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);

  void InverseQrt(TPZVec<REAL> &sol,TPZFMatrix<STATE> &C);
  /**Compute contribution to the stiffness matrix and right hand side at an integration point
     void ContributeRhs(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ef);*/

  /**To evaluate true-error, L2-error and estimate error*/
  void Errors(TPZVec<REAL> &x,TPZVec<STATE> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,TPZVec<STATE> &flux,
        TPZVec<STATE> &uexact,TPZFMatrix<STATE> &duexact,TPZVec<REAL> &val);

  void ApplyPeriodic(TPZCompMesh &cmesh);
  /* 0 -> condicao fronteira periodica, -n depende do valor do ponto x*/
  int IdBC(double *x);
  /**Return number of variables to print into the postprocessing file*/
  int AllVariablesPost() { return 1; }
  /**Return number of variables to print into the postprocessing file*/
  int AllVariablesImplement() { return fOrder; }
};

inline TPZMaterial *TLinearLaw::NewMaterial() {
  TPZMaterial *newmat = new TLinearLaw(Id(),fOrder);
  return newmat;
}

#endif
