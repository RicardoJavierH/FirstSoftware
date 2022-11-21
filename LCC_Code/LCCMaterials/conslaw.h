/*******       FILE : ConsLaw.h   ConsLaw -> Conservation Law

Cont\'em a declara\c c\~ao da classe base "TConservationLaw". 
Esta \'e uma classe abstrata que permite resolver uma lei de conserva\c c\~ao
do tipo hiperb\'olico.
Para esta resolu\c c\~ao emplea-se o m\'etodo do volume finito e/ou
o m\'etodo de elementos finitos.
Aplica tamb\'em um algoritmo de multiresolu\c c\~ao sobre as medias celulares 
para detectar as singularidades da fun\c c\~ao solu\c c\~ao.
A lei de conserva\c c\~ao pode ser uni- ou bi-dimensional e linear ou n\~ao linear.
*******                                                *******/

#ifndef CONSERVATIONLAWH
#define CONSERVATIONLAWH

#include <stdlib.h>
#include <string.h>

#include "pzdiscgal.h"
#include "pzerror.h"

template<class T>
class TPZVec;

class TPZCompMesh;
class TPZInterpolatedElement;

/**
 * \f$ u_t + fC * div( FF(u) ) = s(t,x)  \f$
 */
class TConservationLaw : public TPZDiscontinuousGalerkin {

  /** Alfa used into numerical flux computations */
  REAL fAlfa;
  /** Current Time = Tk*/
  REAL fCurrentTime;
 protected:
  REAL fDeltaT;
    
	REAL fCFLDifusion;
    REAL fCFL;

  /** To store conservation law name*/
  std::string fName;
  /** Store the maxime jacobian eigenvalue over all solution values */
  STATE  fMaxEigen;


 public:

    /**CONSTRUCTORS and DESTRUTOR */
  TConservationLaw(int id);
  TConservationLaw(TConservationLaw &law);

  /** After conservation law creating is necessary to fill default values */
  void SetDefaultData();

  /** To set and put private data*/
  virtual double SonicPoint();
  std::string Name() override { return fName; }
  void SetName(std::string &name);
  REAL Time() { return fCurrentTime; }
  STATE MaxEigen() { return fMaxEigen; }
    
  REAL Alfa() { return fAlfa; }
  void SetAlfa(REAL alfa) { fAlfa = alfa; }
    
    void SetStepTime(REAL deltaT) { fDeltaT = deltaT; }
  void SetTime(double deltat,double time=0.);
  void SetMaxEigen(REAL maxeigen) { fMaxEigen = maxeigen; }
  
	void SetFactor(REAL cfldif) { fCFLDifusion = cfldif; }

  /** return the order of the system conservation law system*/
  int NFluxes() override { return Dimension(); }
    
    // Para el calculo del paso del tiempo utilizando la condicion CFL, para um u e para a malha toda
    STATE MaxValJacobFlux(TPZVec<STATE> &u);
    STATE MaxValJacobFlux(TPZCompMesh *mesh);

    // Retorna el valor del paso de tiempo para satisfazer la condicion CFL: Dt < (Dx / (2p+1)*|maxF`(u)|)
    REAL CFLCondition(TPZCompMesh *mesh);

  virtual void ValJacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &valjacob);
  virtual int ValJacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &Beta);

    //u tem m componentes, e flux(jacob, ...) tem N*m componentes
  virtual void FunctionFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &flux) = 0;
  virtual void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob)=0;
    virtual void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal)=0;
  /**Implement the maxime eigenvalue of the jacobian matrix : n*jacob(u) */
  virtual STATE MaxEigJacob(TPZVec<STATE> &u,TPZVec<REAL> &normal)=0;
  virtual STATE ValEigJacob(TPZVec<STATE> &u,int order,int dim)=0;
	
  virtual void RoeMatrix(TPZVec<STATE> &U,TPZVec<STATE> &Up1,TPZFMatrix<STATE> &Roe);
  // EigRoe must to have Roe matrix, Roe eigenvalues, Roe matrix R and R inverse
  virtual void EigRoeMatrix(TPZVec<STATE> &U,TPZVec<STATE> &Up1,TPZVec<STATE> &EigRoe);
  virtual void ValRoeMatrix(TPZVec<STATE> &U,TPZVec<STATE> &Up1,TPZFMatrix<STATE> &ValRoe);

  /**
 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
 * @param data [in] stores all input data
 * @param weight [in] is the weight of the integration rule
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */
  virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

  /**Compute contribution to the stiffness matrix and right hand side at an integration point
  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &dsol,
	REAL weight,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix
	&ek,TPZFMatrix &ef) { }*/
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition material
	 * @since October 07, 2011
	 */
  virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
  /*Compute contribution to the stiffness matrix and right hand side at the integration point of a boundary*/
 
	/** To compute the tau matrix to diffusive term */
	virtual void Tau(TPZFMatrix<REAL> &jacinv,TPZVec<STATE> &sol,TPZFMatrix<STATE> &tau);
	/** To compute the inverse jacobian matrix */
	virtual void InverseJacob(TPZFMatrix<STATE> &mat);

  /**To evaluate true-error, L2-error and estimate error*/
  void Errors(TPZVec<REAL> &x,TPZVec<STATE> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,TPZVec<STATE> &flux,
        TPZVec<STATE> &uexact,TPZFMatrix<STATE> &duexact,TPZVec<REAL> &val) override;

  /**To clean the conservation law data*/
  virtual void Clean();

  /**Print conservation law information*/
  void Print(std::ostream &out) override;

  /**To compute the contribution over edges of the elements*/
  virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef) override;
    /**
     * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
     * @param data [in]
     * @param dataleft [in]
     * @param weight [in]
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition object
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
        
    }
	
    void SetParameters(int Explicit,REAL CFL,REAL DiffusionCoef);
  void SetCoef(REAL coef) { fCoef = coef; }
  virtual void ApplyPeriodic(TPZCompMesh &cmesh) { }
  void SetPoint(TPZVec<REAL> &x);

  /** Return number of variables to print in postprocessing file 	*/
  virtual int AllVariablesPost() { return NStateVariables(); }
  virtual int AllVariablesImplemented() { return NStateVariables(); }
  /**To return the variables name of the material*/
  virtual void VariablesName(TPZVec<std::string> &names);
    
    /** Para ser usado por el analisis y proveer el orden del metodo numerico
    y el coeficiente que debe multiplicar las contribuciones. Coeficientes del esquema numerico en el tiempo */
    void SetOrderAndCoeficient(int order, REAL Coef);

  /** Factor to diffussive term*/
  virtual int IdBC(double *x) { return 5; }
    
    // Para forçar o cálculo da soluçao e das normais
    virtual void FillDataRequirements(TPZMaterialData &data) override
    {
        TPZDiscontinuousGalerkin::FillDataRequirements(data);
        data.fNeedsSol = true;
        data.fNeedsNormal = true;
        data.fNeedsHSize = true;
    }
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data) override
    {
        TPZDiscontinuousGalerkin::FillBoundaryConditionDataRequirement(type,data);
        data.fNeedsSol = true;
        data.fNeedsNormal = true;
    }


 protected:
  int        fExplicit;
  /** Coefficient to RungeKutta method : cl*fDeltaT */
  REAL       fCoef;
    int fRKOrder;
  /** Point to application*/
  double Point[3];
};

inline void TConservationLaw::SetParameters(int Explicit,REAL CFL,REAL DiffusionCoef) {
    fExplicit = Explicit;
    fCFLDifusion = DiffusionCoef;
    fCFL = CFL;
}
inline int TConservationLaw::ValJacobFlux(TPZVec<STATE> &/*u*/,TPZFMatrix<STATE> &/*valjacob*/,TPZVec<REAL> &/*beta*/) {
  PZError << "TConservationLaw::ValJacobFlux is called." << std::endl;
	return 0;
}
inline void TConservationLaw::ValJacobFlux(TPZVec<STATE> &/*u*/,TPZFMatrix<STATE> &/*valjacob*/) {
  PZError << "TConservationLaw::ValJacobFlux (void) is called." << std::endl;
}
inline void TConservationLaw::RoeMatrix(TPZVec<STATE> &/*u*/,TPZVec<STATE> &/*up1*/,TPZFMatrix<STATE> &/*Roe*/) {
  PZError << "TConservationLaw::RoeMatrix is called." << std::endl;
}
inline void TConservationLaw::EigRoeMatrix(TPZVec<STATE> &/*u*/,TPZVec<STATE> &/*up1*/,TPZVec<STATE> &/*EigRoe*/) {
  PZError << "TConservationLaw::EigRoeMatrix is called." << std::endl;
}
inline void TConservationLaw::ValRoeMatrix(TPZVec<STATE> &/*u*/,TPZVec<STATE> &/*up1*/,TPZFMatrix<STATE> &/*ValRoe*/) {
  PZError << "TConservationLaw::ValRoeMatrix is called." << std::endl;
}

inline void TConservationLaw::Errors(TPZVec<REAL> &x,TPZVec<STATE> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,TPZVec<STATE> &flux,
      TPZVec<STATE> &uexact,TPZFMatrix<STATE> &duexact,TPZVec<REAL> &val) {
  PZError << "TConservationLaw::Errors is called." << std::endl;
}
inline double TConservationLaw::SonicPoint() {
  PZError << "TConservationLaw::SonicPoint is called.\n";
  return 0.;
}
inline void TConservationLaw::SetTime(double deltat,double time) {
    fCurrentTime = time;
	fDeltaT = deltat;
}

#endif
