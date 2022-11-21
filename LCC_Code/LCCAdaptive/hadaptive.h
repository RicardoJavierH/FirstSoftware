#ifndef MYHPADAPTIVEHH
#define MYHPADAPTIVEHH

#include <iostream>
#include "pzcompel.h"
#include "pzstack.h"
#include "ErrorEstimator.h"

class TAdaptive {
	
	TPZCompMesh *fCMesh;
	int fMaxLevelRefinements;
	int fMaxPRefinements;
    bool fGDFEM;

 public:
  TAdaptive(TPZCompMesh *cmesh, int maxlevel, int maxp,bool DGFEM);
	~TAdaptive();

    int GetMaxLevelRefinements() { return fMaxLevelRefinements; }
	void SetMaxLevelRefinements(int maxlevel) { fMaxLevelRefinements = maxlevel; }

	int GetMaxPRefinements() { return fMaxLevelRefinements; }
	void SetMaxPRefinements(int maxp) { fMaxLevelRefinements = maxp; }

	void DecrementOrder(TPZStack<int64_t> &Elements);
  void IncrementOrder(TPZStack<int64_t> &Elements);

  void Refinements(TPZStack<int64_t> &Elements);

  void Coarsing(TPZStack<int64_t> &Elements);
	void Adapting(TErrorEstimator *estimator);
	
};

#endif
