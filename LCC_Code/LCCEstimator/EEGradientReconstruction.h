#ifndef ERRORESTIMATORGRADIENTRECONSTRUCTIONHPP
#define ERRORESTIMATORGRADIENTRECONSTRUCTIONHPP

#include <iostream>
#include "ErrorEstimator.h"


class TEEGradientReconstructionLS : TErrorEstimator {
	
public:
	TEEGradientReconstructionLS(TPZCompMesh *cmesh, REAL majorc, REAL mediumc, REAL minorc, int var);
	~TEEGradientReconstructionLS();
	
	/* Function as error estimator for each element */
	virtual REAL ErrorEstimator(TPZCompEl *cel,int var);
	virtual REAL ErrorEstimator(TPZCompEl *cel,TPZFMatrix<REAL> &grad,int var=0);
};

#endif
