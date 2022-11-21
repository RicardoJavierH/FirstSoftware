//
//  SimilarUniformRefinements.cpp
//  PZ
//
//  Created by labmec on 08/10/17.
//
//

#include <stdio.h>
#include "pzcmesh.h"
#include "pzcompel.h"
#include "Refinements.h"

/* 2. Functions to uniform refinement of the geometric meshes.
 Projects: Poisson3D_Shock
 */
// The function refines all geometric elements (without subelements) for any material or for one matidtodivided material - Only for one material
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh) {
    if(nDiv < 1 || !gmesh || gmesh->Dimension()<0) {
        std::cout << "It is nothing done! (UniformRefinement)." << std::endl;
        return;
    }
    
    TPZManVector<TPZGeoEl*> filhos;
    // Loop over (number of refinements) nDiv - It refines all geometric elements without subelements
    for(int D=0; D<nDiv; D++)
    {
        int64_t nels = gmesh->NElements();
        for(int64_t elem = 0; elem < nels; elem++)
        {
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            if(!gel || gel->HasSubElement() || !gel->Dimension())
                continue;
            else
                gel->Divide(filhos);
        }
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}

// The function refines all geometric elements (without subelements) for any material or for one matidtodivided material - Only for one material
void UniformHRefinementSelection(const int nDiv, TPZGeoMesh* gmesh, TPZVec<int64_t>& indexsubelements) {
	if (!indexsubelements.NElements() || !gmesh) {
		return;
	}
	if (nDiv != 1) {
		std::cout << "\nNeed implementation.\n"; return;
	}

	TPZManVector<TPZGeoEl*> filhos;
	int64_t elem, nels = indexsubelements.NElements();
	for (elem = 0; elem < nels; elem++) {
		int64_t index = indexsubelements[elem];
		if (index < 0) continue;
		TPZCompEl* cel = gmesh->Reference()->ElementVec()[index];
		TPZGeoEl* gel = cel->Reference();
		if (!gel || gel->HasSubElement() || !gel->Dimension())
			continue;
		else
			gel->Divide(filhos);
	}
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}
void CheckNeighbours(TPZCompMesh *cmesh,TPZVec<int64_t> &filhos,bool interpolated) {
    int nsons = filhos.NElements();
    if(!nsons)
        return;
    TPZManVector<int64_t> sonsfilhos;
    for(int i=0;i<nsons;i++) {
        TPZCompEl *cel = cmesh->ElementVec()[filhos[i]];
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
//        for(int j=gel->NCornerNodes();j<gel->NSides()-1;j++) {
        for(int j=0;j<gel->NSides()-1;j++) {
            TPZGeoElSide gels(gel,j), neighs;
            int count = 0;
            neighs = gel->Neighbour(j);
            int level = gel->Level();
            TPZCompElSide neighlower;
            neighlower = gels.LowerLevelCompElementList2(1);
            if(!(neighlower.Element()))
                continue;
            int neighlevel = neighlower.Element()->Reference()->Level();
            if(neighlevel+1<level)
                cmesh->Divide(neighlower.Element()->Index(),sonsfilhos,1);
        }
    }
}
void UniformHRefineCoarsenSelection(TPZCompMesh* cmesh, TPZVec<int64_t>& indexfathers,int maxhref) {
	if (!cmesh) {
		return;
	}
//	TPZMatrix<STATE> *matriz = cmesh->Block().Matrix();
	cmesh->Block().SetMatrix(&(cmesh->Solution()));
	
	TPZManVector<int64_t> filhos;
	int64_t index;
	int64_t elem, nels = indexfathers.NElements();
	for (elem = 0; elem < nels; elem++) {
		index = indexfathers[elem];
		TPZCompEl* cel = cmesh->ElementVec()[index];
		if (!cel || !cel->Mesh() || !cel->Dimension() || cel->IsInterface() || cel->Reference()->Level()>maxhref)
			continue;
		else {
			TPZGeoEl* father = cel->Reference()->Father();
			if(!father) {   // Quando nao tem pai, refina ese elemento (esta no nivel zero
				cmesh->Divide(index, filhos, 1);
                CheckNeighbours(cmesh,filhos,1);
			}
			else {   // Quando tem pai, refina o elemento e todos seus irmaos (filhos do mesmo nivel)
				int i, nsub = father->NSubElements();
				for (i = 0; i < nsub; i++) {
					TPZInterpolationSpace* subcel = (TPZInterpolationSpace*)(father->SubElement(i)->Reference());
					if(!subcel)
						continue;
					index = subcel->Index();
//					if(!(subcel->Reference()->HasSubElement())) {
						cmesh->Divide(index, filhos, 1);
						CheckNeighbours(cmesh,filhos,1);
//					}
				}
			}
		}
	}
/*    cmesh->Reference()->ResetConnectivities();
    cmesh->Reference()->BuildConnectivity();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->InitializeBlock();*/

//	cmesh->Reference()->ResetConnectivities();
//	cmesh->Reference()->BuildConnectivity();
//	TPZFMatrix<STATE> Sol = cmesh->Solution();
//	cmesh->AutoBuild();
//	cmesh->AdjustBoundaryElements();
//	cmesh->CleanUpUnconnectedNodes();
}
/*void UniformHRefineCoarsenSelection(TPZCompMesh* cmesh, TPZVec<TPZCompEl *>& indexfathers,TPZVec<TPZVec<int64_t> >& indexsubelements, bool DGFEM) {
	if (!cmesh) {
		return;
	}
	
	TPZManVector<int64_t> filhos;
	int64_t index;
	int64_t elem, nels = indexfathers.NElements();
	for (elem = 0; elem < nels; elem++) {
		TPZCompEl* cel = indexfathers[index];
		if (!cel || !cel->Dimension())
			continue;
		else {
			TPZGeoEl* father = cel->Reference()->Father();
			if(!father)   // Quando nao tem pai, refina ese elemento (esta no nivel zero
				((TPZInterpolatedElement*)cel)->Divide(index, filhos, 0);
			else {   // Quando tem pai, refina o elemento e todos seus irmaos (filhos do mesmo nivel)
				int i, nsub = father->NSubElements();
				for (i = 0; i < nsub; i++) {
					TPZInterpolatedElement* subcel = (TPZInterpolatedElement*)(father->SubElement(i)->Reference());
					index = subcel->Index();
					subcel->Divide(index, filhos, 0);
				}
			}
		}
	}
	cmesh->Reference()->ResetConnectivities();
	cmesh->Reference()->BuildConnectivity();
	nels = indexsubelements.NElements();
	for (elem = 0; elem < nels; elem++) {
		cmesh->Coarsen(indexsubelements[elem], index, DGFEM);
	}
	cmesh->Reference()->ResetConnectivities();
	cmesh->Reference()->BuildConnectivity();
}*/
void UniformHRefineCoarsenSelection(TPZCompMesh* cmesh, TPZVec<int64_t>& indexfathers,TPZVec<TPZVec<int64_t> >& indexsubelements, bool DGFEM) {
	if (!cmesh) {
		return;
	}

	TPZManVector<int64_t> filhos;
	int64_t index;
	int64_t elem, nels = indexfathers.NElements();
	for (elem = 0; elem < nels; elem++) {
		index = indexfathers[elem];
		if (index < 0) continue;
		TPZCompEl* cel = cmesh->ElementVec()[index];
		if (!cel || !cel->Dimension())
			continue;
		else {
			TPZGeoEl* father = cel->Reference()->Father();
			if(!father)   // Quando nao tem pai, refina ese elemento (esta no nivel zero
				((TPZInterpolatedElement*)cel)->Divide(index, filhos, 1);
			else {   // Quando tem pai, refina o elemento e todos seus irmaos (filhos do mesmo nivel)
				int i, nsub = father->NSubElements();
				for (i = 0; i < nsub; i++) {
					TPZInterpolatedElement* subcel = (TPZInterpolatedElement*)(father->SubElement(i)->Reference());
					index = subcel->Index();
					subcel->Divide(index, filhos, 1);
				}
			}
		}
	}
	cmesh->Reference()->ResetConnectivities();
	cmesh->Reference()->BuildConnectivity();
	nels = indexsubelements.NElements();
	for (elem = 0; elem < nels; elem++) {
		cmesh->Coarsen(indexsubelements[elem], index, DGFEM);
	}
	cmesh->Reference()->ResetConnectivities();
	cmesh->Reference()->BuildConnectivity();
}
void UniformPRefinementSelection(TPZCompMesh* cmesh,TPZVec<int64_t>& indextoprefine,int MaxPRefine) {
	if (!indextoprefine.NElements() || !cmesh) {
		return;
	}
	int64_t elem, nels = indextoprefine.NElements();
	for (elem = 0; elem < nels; elem++) {
		if (elem > 3) break;
		int64_t index = indextoprefine[elem];
		if (index < 0) continue;
		TPZInterpolatedElement* cel = (TPZInterpolatedElement*)cmesh->ElementVec()[index];
		if (!cel || cel->IsInterface())
			continue;
		int order = cel->GetPreferredOrder();
        if(order>MaxPRefine-1) continue;
		cel->PRefine(order + 1);
	}
/*	cmesh->CleanUpUnconnectedNodes();
	cmesh->InitializeBlock();
    cmesh->ExpandSolution();*/
}
// The function refines all geometric elements (without subelements) for any material or for material ids in MatIdsVec. If dim < 0 will be refining all elements of the mesh (all dimension not zero)
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, TPZVec<int> *MatIdsVec) {
    if(!nDiv || !gmesh || !dim || dim>gmesh->Dimension()) {
        std::cout << "It is nothing done! (UniformRefinement)." << std::endl;
        return;
    }
    // If dim is negative, the user want to refine geometrical elements not points
    if(dim<0) {
        UniformRefinement(nDiv,gmesh);
        return;
    }
    
    // If MatIdsVector is null, the refinement will be made for all materials
    bool allmaterial = false;
    if(!MatIdsVec || !MatIdsVec->NElements()) {
        allmaterial = true;
    }
    
    int matidtodivided;
    TPZManVector<TPZGeoEl*> filhos;
    // Loop over (number of refinements) nDiv - It refines all geometric elements without subelements
    for(int D=0; D<nDiv; D++)
    {
		int64_t nels = gmesh->NElements();
        for(int64_t elem = 0; elem < nels; elem++)
        {
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            int geldim = gel->Dimension();
            if(!gel || gel->HasSubElement() || !geldim || geldim != dim)
                continue;
            if(!allmaterial) {
                int matidcurrent = gel->MaterialId();
                for(int count=0;count<MatIdsVec->NElements();count++) {
                    matidtodivided = (*MatIdsVec)[count];
                    if(matidcurrent == matidtodivided) {
                        gel->Divide(filhos);
                        break;
                    }
                }
            }
            else{
                gel->Divide(filhos);
            }
        }
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}


void RegularizeMesh(TPZGeoMesh *gmesh, int dimension)
{
    //Control flag
    bool changed = true;
    // If exists something wrong
    if(!gmesh || gmesh->Dimension() < 0)
        DebugStop();

    while (changed)
    {
        changed = false;
		int64_t nel = gmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->ElementVec()[el];
            if (gel->HasSubElement()) {
                continue;
            }
            int dim = gel->Dimension();
            if (dim != dimension)
            {
                continue;
            }
            int nsides = gel->NSides();
			int64_t nrefined = 0;
			int nsidedim = 0;
            for (int is=0; is<nsides; is++) {
                int sidedim = gel->SideDimension(is);
                if (sidedim != dim-1) {
                    continue;
                }
                nsidedim++;
            }
            for (int is=0; is<nsides; is++) {
                int sidedim = gel->SideDimension(is);
                if (sidedim != dim-1) {
                    continue;
                }
                TPZGeoElSide thisside(gel,is);
                TPZGeoElSide neighbour = thisside.Neighbour();
                if (neighbour != thisside) {
                    TPZStack<TPZGeoElSide> subelements;
                    neighbour.GetSubElements2(subelements);
                    int nsub = subelements.size();
                    if (nsub > 0) {
                        nrefined++;
                    }
                    for (int isub=0; isub<nsub; isub++) {
                        TPZGeoElSide sub = subelements[isub];
                        if (sub.Dimension() != dim-1) {
                            continue;
                        }
                        if (sub.HasSubElement()) {
                            TPZManVector<TPZGeoEl *> newsub;
                            gel->Divide(newsub);
                            changed = true;
                            break;
                        }
                    }
                }
                if (gel->HasSubElement()) {
                    break;
                }
            }
            if (nrefined >= nsidedim-1) {
                TPZManVector<TPZGeoEl *> newsub;
                gel->Divide(newsub);
                changed = true;
            }
        }
    }
	gmesh->CleanUp();
	gmesh->BuildConnectivity();
}

