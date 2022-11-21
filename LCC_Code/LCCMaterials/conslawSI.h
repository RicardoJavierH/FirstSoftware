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

#ifndef CONSERVATIONSILAWH
#define CONSERVATIONSILAWH

#include <stdlib.h>
#include <string.h>

#include "conslaw.h"
#include "pzerror.h"

template<class T>
class TPZVec;

class TPZCompMesh;
class TPZInterpolatedElement;

/**
 * \f$ u_t + fC * div( FF(u) ) = s(t,x)  \f$
 */
class TConservationLawSI : public TConservationLaw {


 public:

    /**CONSTRUCTORS and DESTRUTOR */
  TConservationLawSI(int id);
  TConservationLawSI(TConservationLawSI &law);

  /** After conservation law creating is necessary to fill default values */
  void SetDefaultData();

    // Para el calculo del paso del tiempo utilizando la condicion CFL, para um u e para a malha toda
    STATE MaxValJacobFlux(TPZVec<STATE> &u);
    STATE MaxValJacobFlux(TPZCompMesh *mesh);

  /**Implement the maxime eigenvalue of the jacobian matrix : n*jacob(u) */
	
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
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override { }
	
};

#endif
