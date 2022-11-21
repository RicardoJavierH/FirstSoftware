#include <stdlib.h>

#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzvec.h"

#include "ErrorEstimator.h"

TErrorEstimator::TErrorEstimator(TPZCompMesh *cmesh, REAL majorc, REAL mediumc, REAL minorc, int var) {
	fVariable = var;
	fMajorCoeff = majorc;
	fMediumCoeff = mediumc;
	fMinorCoeff = minorc;
	fCMesh = cmesh;
}
TErrorEstimator::~TErrorEstimator() {
}

/* Function to determining the classification of the all computational elements */
void TErrorEstimator::ApplyingEstimation() {
	int64_t nelems = fCMesh->NElements();
	int64_t i, index;
	REAL estimator;
	CleaningElementVectors();
	TPZVec<int> FlagElements(nelems,0);
	TPZFMatrix<REAL> gradients(fCMesh->Dimension(),1);
	for(i=0;i<nelems;i++) {
		TPZCompEl *cel = fCMesh->ElementVec()[i];
		if(!cel || cel->Dimension() != fCMesh->Dimension() || cel->IsInterface())
			continue;
		index = cel->Index();
		if(FlagElements[index]) continue;
		
		estimator = ErrorEstimator(cel,gradients,0);
		
		FlagElements[index] = 1;
		if(estimator < fMinorCoeff)
			fElementsBelowMinorLimit.Push(cel->Index());
		else if(estimator < fMediumCoeff)
			fElementsBetweenMiddleMinorLimit.Push(cel->Index());
		else if(estimator < fMajorCoeff)
			fElementsBetweenMajorMiddleLimit.Push(cel->Index());
		else
			fElementsAboveMajorLimit.Push(cel->Index());
	}
}
void TErrorEstimator::CleaningElementVectors() {
	fElementsAboveMajorLimit.clear();
	fElementsBelowMinorLimit.clear();
	fElementsBetweenMajorMiddleLimit.clear();
	fElementsBetweenMiddleMinorLimit.clear();
}
/*
void Func() {
	// Computing approximation of gradient
	TPZFMatrix<REAL> gradients(dim,1);
	for(int i=0;i<cmesh->NElements();i++) {
		TPZCompEl *cel = cmesh->ElementVec()[i];
		if(!cel || cel->Dimension()!=dim) continue;
		gradients.Zero();
		//					if(!discont)
		GradReconstructionByLeastSquares(cel,gradients,0);
		//					else
		//						GradReconstructionByLeastSquares_Self(cel,gradients,0);
		cout << "Element " << i << " Grad = (";
		for(int j=0;j<dim;j++) {
			cout << gradients[j];
			if(j<dim-1) cout << ",";
		}
		cout << ")\n";
	}

}
*/
