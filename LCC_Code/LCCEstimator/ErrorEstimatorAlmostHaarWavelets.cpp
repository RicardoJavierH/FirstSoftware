#include <stdlib.h>

#include "pzcmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include "pzvec.h"

#include "ErrorEstimatorAlmostHaarWavelets.h"

TErrorEstimatorAlmostHaarWavelets::TErrorEstimatorAlmostHaarWavelets(TPZCompMesh *cmesh, REAL majorc, REAL mediumc, REAL minorc, int var) : TErrorEstimator(cmesh,majorc,mediumc,minorc,var) {
}
TErrorEstimatorAlmostHaarWavelets::~TErrorEstimatorAlmostHaarWavelets() {
}

/* Function as error estimator for each element */
REAL TErrorEstimatorAlmostHaarWavelets::ErrorEstimator(TPZCompEl *cel,int var) {
	TPZGeoEl *father = cel->Reference()->Father();
	if(!father) return 0.;
	int i, nsons = father->NSubElements();
	if(!nsons) return 100.;
	TPZVec<TPZCompEl *> Sons(nsons,0);
	TPZVec<REAL> MeansBySubElement(nsons,0.);
	for(i=0;i<nsons;i++) {
		Sons[i] = father->SubElement(i)->Reference();
		if(!Sons[i])
			DebugStop();
		MeansBySubElement[i] = ((TPZInterpolatedElement *)(Sons[i]))->MeanSolution(fVariable);
	}
	REAL mean = 0., diff = 0.;
	for(i=0;i<nsons;i++) {
		mean += MeansBySubElement[i];
	}
	mean /= nsons;
	for(i=0;i<nsons;i++) {
		REAL diffi = mean-MeansBySubElement[i];
		if(diff < diffi)
			diff = diffi;
	}
	return diff;
}
