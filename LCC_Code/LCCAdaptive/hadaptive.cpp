#include <stdlib.h>
#include <math.h>

#include "TPZInterfaceEl.h"

#include "commonJC.h"
#include "Refinements.h"

#include "hadaptive.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzintel.h"

TAdaptive::TAdaptive(TPZCompMesh *cmesh, int maxlevel, int maxp,bool DGFEM) {
	fCMesh = cmesh;
    fGDFEM = DGFEM;
	fMaxLevelRefinements = maxlevel;
	fMaxPRefinements = maxp;
}
TAdaptive::~TAdaptive() {
}

/** idgfather has the index of the elements that must to be refined */
void TAdaptive::Refinements(TPZStack<int64_t> &Elements) {
	int nelem = Elements.NElements(), dim = fCMesh->Dimension();
  if(!nelem) return;
	int64_t index;
	TPZCompEl *cel;
  int64_t i;
	TPZManVector<int64_t> subindex;

  for(i=0;i<nelem;i++) {
	  index = Elements.Pop();
	  cel = fCMesh->ElementVec()[index];
    if(!cel || cel->IsInterface() || cel->Dimension() < dim)
      continue;
	  if(cel->Reference()->Level() > fMaxLevelRefinements-1)
		  continue;

	  fCMesh->Divide(index,subindex,1);
  }
	fCMesh->AdjustBoundaryElements();
}

/** Necessita testar nivel de refinamento para fazer elemento coarse(grosso)*/
void TAdaptive::Coarsing(TPZStack<int64_t> &Elements) {
  int64_t i, index, nelem = Elements.NElements();
	int j;
  if(!nelem) return;
	int nsub, order, dim = fCMesh->Dimension(), orderold = fCMesh->GetDefaultOrder();
	bool coarse = false;
	TPZCompEl *cel;
  for(i=0;i<nelem;i++) {
	  index = Elements.Pop();
	  if(index < 0) continue;
	  cel = fCMesh->ElementVec()[index];
    if(!cel || cel->IsInterface() || cel->Dimension() < dim) {
		nsub = 1;
		continue;
    }
	 else {
		 order = orderold+1;
		 TPZGeoEl *gel = cel->Reference()->Father();
		 if(!gel || gel->Level()==1)
			 continue;
		 nsub = gel->NSubElements();
		 TPZVec<int64_t> indexes(nsub,-1);
		 /** Verificando se os elementos sucessivos tem o mesmo elemento pai
		  e averiguando se existe elemento com ordem de interpolacao zero */
		 for(j=0;j<nsub;j++) {
			 if(!(gel->SubElement(j)->Reference()))
				 break;
			 int64_t k, ind = gel->SubElement(j)->Reference()->Index();
			 for(k=nelem-i-1;k>-1;k--)
				 if(ind==Elements[k])
					 break;
			 if(k==-1)
				 break;
			 Elements[k] = -1;
			 indexes[j] = ind;
			 cel = fCMesh->ElementVec()[ind];
			 int nsides = cel->Reference()->NSides();
             int orderpartial;
			 if(!fGDFEM) orderpartial = ((TPZInterpolatedElement *)cel)->PreferredSideOrder(nsides-1);
             else orderpartial = order;
			 /**Verificando se os elementos sucessivos tem o mesmo elemento pai
			  e averiguando se existe elemento com ordem de interpolacao zero*/
			if(order>orderpartial) order = orderpartial;
		 }
		 if(j!=nsub)
			 continue;
		 TPZCompEl::SetgOrder(order);
		 fCMesh->Coarsen(indexes,index,fGDFEM);
		 fCMesh->ExpandSolution();
		 coarse = true;
	 }
  }
	if(coarse) {
			TPZCompEl::SetgOrder(orderold);
		fCMesh->AdjustBoundaryElements();
		fCMesh->CleanUpUnconnectedNodes();
			fCMesh->InitializeBlock();
	}
}

void TAdaptive::Adapting(TErrorEstimator *estimator) {
	estimator->ApplyingEstimation();
	UniformHRefineCoarsenSelection(fCMesh,estimator->fElementsAboveMajorLimit,fMaxLevelRefinements);
	Coarsing(estimator->fElementsBelowMinorLimit);
//    UniformPRefinementSelection(fCMesh,estimator->fElementsBetweenMiddleMinorLimit,fMaxPRefinements);
//    fCMesh->Reference()->ResetConnectivities();
//    fCMesh->Reference()->BuildConnectivity();
    fCMesh->AdjustBoundaryElements();
    fCMesh->ExpandSolution();
    fCMesh->CleanUpUnconnectedNodes();
//    fCMesh->InitializeBlock();

}

/** Before to refine a group of geo elements, decrement its interpolation order*/
void TAdaptive::DecrementOrder(TPZStack<int64_t> &Elements) {
/*  int64_t i, j, nelem = elfathersvec.NElements();
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  int side, order, minorder=0;
  TPZInterpolatedElement *el, *newel;
  TPZBlock<STATE> &block = cmesh.Block();
    int oldgorder = TPZCompEl::GetgOrder();
    int64_t index;
	if(!TPZCompEl::GetgOrder() ) {
	  elfathersvec.Resize(0);
		return 0;
	}
  for(i=0;i<nelem;i++) {
    el = (TPZInterpolatedElement *)elvec[elfathersvec[i]];
    side = el->Reference()->NSides();
    order = el->PreferredSideOrder(side-1);
    if(el->Reference()->Level()< MaxLevelRefinements-1 || !order) continue;
    / * Aqui a ordem tem que ser escolhida: se order>1 entao 1 se order =1 entao 0 * /
    if(order==1) {
      for(j=0;j<side;j++) {
        / * Sera que esta-se limpando os blocos com os valores das variaveis??? * /
        if(el->NSideConnects(j)) {
          PZError << "DecrementOrder. Element has connect continuous.\n";
          break;
        }
      }
      if(j<side) continue;
      minorder = 0;
      int nvar = el->Material()->NStateVariables();
      TPZVec<REAL> values(nvar);
      for(j=0;j<nvar;j++) values[j] = el->MeanSolution(j);
      el->PRefine(0);
      int seqnum = el->Connect(side).SequenceNumber();
      if(nvar!=block.Size(seqnum))
        PZError << "DecrementOrder. Dimension of block internal is uncompatible.\n";
      for(j=0;j<nvar;j++) block(seqnum,0,j,0) = values[j];
    }
    else {
      TPZCompEl::SetgOrder(order-1);
      TPZGeoEl *gel = el->Reference();
      gel->ResetReference();
        cmesh.CreateCompEl(gel,index);
//      gel->CreateCompEl(cmesh,index);
      newel = (TPZInterpolatedElement *)(cmesh.ElementVec()[index]);
      newel->CheckConstraintConsistency();
      cmesh.ExpandSolution();
      newel->InterpolateSolution(*el);
      index = el->Index();
      delete el;
      cmesh.ElementVec()[index] = 0;
      cmesh.ElementVec().SetFree(index);
      gel->SetReference(newel);
    }
  }
  TPZCompEl::SetgOrder(oldgorder);
  return minorder;*/
}

/** After to coarsing a group of geo elements, increment its interpolation order*/
void TAdaptive::IncrementOrder(TPZStack<int64_t> &Elements) {
/*  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  int i, nelem = elsons.NElements();
  int side, order, maxorder = 0;
  TPZInterpolatedElement *el;
  for(i=0;i<nelem;i++) {
    el = (TPZInterpolatedElement *)elvec[elsons[i]];
    side = el->Reference()->NSides()-1;
    order = el->PreferredSideOrder(side)+1;
    if(order <= el->PreferredSideOrder(side)) {
      el->PRefine(order);
      if(maxorder<order) maxorder = order;
    }
  }
  return maxorder;*/
}

