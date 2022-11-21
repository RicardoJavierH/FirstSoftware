
#ifndef Refinements_hpp
#define Refinements_hpp

#include <time.h>
#include <stdio.h>
#include <math.h>

#include <fstream>
#include <cmath>

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"
#include "pzgeotetrahedra.h"
#include "TPZGeoCube.h"
#include "pzgeopyramid.h"
#include "pzgeoprism.h"

#include "pzgeoelbc.h"

#include "pzlog.h"
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzcheckgeom.h"
#include "pzcheckmesh.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZFrontStructMatrix.h"
#include "TPZParFrontStructMatrix.h"

#include "TPZParSkylineStructMatrix.h"
#include "pzsbstrmatrix.h"
#include "pzfstrmatrix.h"

#include "TPZMaterial.h"
#include "pzbndcond.h"

#include "pzfunction.h"

#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include "pzshapelinear.h"

#include "TPZRefPatternTools.h"


#include "pzbuildmultiphysicsmesh.h"
#include "pzcondensedcompel.h"

using namespace std;
using namespace pzshape;
using namespace pzgeom;



/* 2. Functions to uniform refinement of the geometric meshes.
 Projects: Poisson3D_Shock
 */

/** Fucntions to apply refinement. */
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh);
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, TPZVec<int> *MatIdsVec=NULL);
void UniformHRefinementSelection(const int nDiv, TPZGeoMesh* gmesh, TPZVec<int64_t>& indextohrefine);
void UniformHRefineCoarsenSelection(TPZCompMesh* cmesh, TPZVec<int64_t>& indextorefine, TPZVec<TPZVec<int64_t> >& indexestocoarse,bool DGFEM);
//void UniformHRefineCoarsenSelection(TPZCompMesh* cmesh, TPZVec<TPZCompEl *>& indextorefine, TPZVec<TPZVec<int64_t> >& indexestocoarse,bool DGFEM);
void UniformHRefineCoarsenSelection(TPZCompMesh* cmesh, TPZVec<int64_t>& indextorefine,int maxhref);
void UniformPRefinementSelection(TPZCompMesh* cmesh,TPZVec<int64_t>& indextoprefine,int MaxPRefine);

void RegularizeMesh(TPZGeoMesh *gmesh,int dimension);
void PrintNRefinementsByType(int nref, int64_t nels,int64_t newnels,TPZVec<int64_t> &counter,std::ostream &out = std::cout);

#endif
