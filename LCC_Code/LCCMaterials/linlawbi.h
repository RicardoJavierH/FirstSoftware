/*******       File : LinLaw2.h       *******/

#ifndef LINLAWBIH
#define LINLAWBIH

#include "linlaw.h"

class TLinearLaw2D : public TLinearLaw {
 protected:
    TPZFMatrix<STATE> fB;   // B and/or C -> convection matrix
    TPZFMatrix<STATE> fValB;// |B| -> |B|=R|EigVal|(InvR)
    REAL fMaxEigValB;              // Valor do maximo autovalor da matrix B
    TPZVec<STATE> fEigValB;         // Vetor de autovalores
    TPZFMatrix<STATE> fEigVectB;        // R -> matrix de autovetores
    TPZFMatrix<STATE> fInvEigVectB;     // Inversa de R

 public:
  TPZFMatrix<STATE> &B() { return fB; }
  TPZFMatrix<STATE> &ValB() { return fValB; }
  TPZVec<STATE> &EigValB(){ return fEigValB; }
  TPZFMatrix<STATE> &EigVectB() { return fEigVectB; }
  TPZFMatrix<STATE> &InvEigVectB() { return fInvEigVectB; }

  void FunctionFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &flux);
  void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob);
	void JacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal);
  void ValJacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &valjacob);
  int ValJacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &Beta);

  TLinearLaw2D(int id,int ord);
  TLinearLaw2D(int id,int order,int type);
  ~TLinearLaw2D() { }

  virtual void SetMatrix(TPZFMatrix<STATE> &A,TPZFMatrix<STATE> &B,TPZFMatrix<STATE> &EigVect,TPZFMatrix<STATE> &EigVectB,TPZFMatrix<STATE> &InvEigVect,TPZFMatrix<STATE> &InvEigVectB,
  TPZVec<STATE> &EigVal,TPZVec<STATE> &EigValB);
  void SetMatrix(std::istream &input);
  void RequireMatrix();
  void Print(std::ostream &out= std::cout);
  STATE MaxEigJacob(TPZVec<STATE> &u,TPZVec<REAL> &normal);
	STATE ValEigJacob(TPZVec<STATE> &u,int order,int dim);
	/** @brief Returns the integrable dimension of the material */
	virtual int Dimension() const { return 2; }

  TPZMaterial *NewMaterial();
  /**Read data of the material from a istream (file data)*/
  void SetData(std::istream &data);
  
  /* 0 -> condicao fronteira periodica, -n depende do valor do ponto x*/
  int IdBC(double *x);
};

inline TPZMaterial *TLinearLaw2D::NewMaterial() {
  TPZMaterial *newmat = new TLinearLaw2D(Id(),fOrder);
  return newmat;
}

class TLinearLaw2DCircle : public TLinearLaw2D {
  REAL   fAngle;   // Angulo de giro solicitado pelo usuario
 public:
  TLinearLaw2DCircle(int id);
  TLinearLaw2DCircle(int id,REAL angle,int type);
  ~TLinearLaw2DCircle() { }

  void SetData(std::istream &input);
  void SetData(REAL angle);
  
  TPZMaterial *NewMaterial();
};

inline TPZMaterial *TLinearLaw2DCircle::NewMaterial() {
  TPZMaterial *newmat = new TLinearLaw2DCircle(Id());
  return newmat;
}

#endif

