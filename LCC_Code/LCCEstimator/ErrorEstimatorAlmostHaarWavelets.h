#ifndef ERRORESTIMATORALMOSTHAARWAVELETSHPP
#define ERRORESTIMATORALMOSTHAARWAVELETSHPP

#include <iostream>
#include "ErrorEstimator.h"


class TErrorEstimatorAlmostHaarWavelets : TErrorEstimator {
	
public:
	TErrorEstimatorAlmostHaarWavelets(TPZCompMesh *cmesh, REAL majorc, REAL mediumc, REAL minorc, int var);
	~TErrorEstimatorAlmostHaarWavelets();
	
	/* Function as error estimator for each element */
	virtual REAL ErrorEstimator(TPZCompEl *cel,int var);
	virtual REAL ErrorEstimator(TPZCompEl *cel,TPZFMatrix<REAL> &grad,int var) {
		return 0.;
	}
};

#endif
