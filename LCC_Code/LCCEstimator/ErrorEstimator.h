#ifndef ERRORESTIMATORBASEHPP
#define ERRORESTIMATORBASEHPP

#include <iostream>
#include "pzstack.h"


class TErrorEstimator {
	
	REAL fMajorCoeff;
	REAL fMediumCoeff;
	REAL fMinorCoeff;
	
	TPZCompMesh *fCMesh;

protected:
	int fVariable;

public:
	TErrorEstimator(TPZCompMesh *cmesh, REAL majorc, REAL mediumc, REAL minorc, int var);
	~TErrorEstimator();
	
	TPZStack<int64_t> fElementsAboveMajorLimit;
	TPZStack<int64_t> fElementsBetweenMajorMiddleLimit;
	TPZStack<int64_t> fElementsBetweenMiddleMinorLimit;
	TPZStack<int64_t> fElementsBelowMinorLimit;
	
	/* Function to determining the classification of the all computational elements */
	void ApplyingEstimation();
	/* To cleaning all the vectors previously selected */
	void CleaningElementVectors();
	
	/* Function as error estimator for each element */
	virtual REAL ErrorEstimator(TPZCompEl *cel,int var) = 0;
	virtual REAL ErrorEstimator(TPZCompEl *cel,TPZFMatrix<REAL> &grad,int var=0) {
		return 0.;
	}
};

#endif
