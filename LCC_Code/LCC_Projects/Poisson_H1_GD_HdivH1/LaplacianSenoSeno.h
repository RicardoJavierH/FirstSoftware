#ifndef STEKLOVPROBLEM
#define STEKLOVPROBLEM

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZMultiphysicsCompMesh.h"

#include "pzvec.h"
#include "pzfmatrix.h"

// Relative to the meshes
// To create geometrical and computational meshes.
TPZGeoMesh* CreatingGeometricMesh(int elementIdMat, int bcdirichlet, int typeel, int &DimProblem);
TPZGeoMesh* CreatingGeometricMeshFromGIDFile(std::string &GeoGridFile, int ElementIDMat, int bcdirichlet, int typeel, int &DimProblem, int &nref);
void CreatingComputationalMesh(TPZCompMesh* cMesh, int elementIdMat, int bcdirichlet, bool DGFEM, int DimProblem,int porder, TPZAutoPointer<TPZFunction<STATE> > &source, TPZAutoPointer<TPZFunction<STATE> > &solExata);

void CMeshFlux(TPZCompMesh* cmesh, int elementIdMat, int DimProblem, int pOrder, TPZAutoPointer<TPZFunction<STATE> > &solExata);
void CMeshPressure(TPZCompMesh *cmesh, int elementIdMat, bool DGFEM, int DimProblem, int pOrder);
void MalhaCompMultifisica(TPZMultiphysicsCompMesh* mphysics, int ElementIDMat, int bcdirichlet, bool DGFEM, int DimProblem, int pOrder, TPZAutoPointer<TPZFunction<STATE> > &source, TPZAutoPointer<TPZFunction<STATE> > &solExata);

// Exact solution of the differential equation  (Laplacian)
// To Model with exact solution Sin*Sin*Sin
void SolExactSeno(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SolExactSenoSeno(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SolExactSenoSenoSeno(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
// Exact solution of the differential equation  (Laplacian)
void SourceFunctionSin1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void SourceFunctionSin2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void SourceFunctionSin3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionSin1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionSin2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionSin3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
// To Model with exact solution ArcTg - It has strong gradient
void SolExactArcTg1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SolExactArcTg2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SolExactArcTg3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
// Exact solution of the differential equation  (Laplacian)
void SourceFunctionArcTg1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void SourceFunctionArcTg2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void SourceFunctionArcTg3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionArcTg1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionArcTg2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionArcTg3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
// To Model with exact solution Sin*Cos*ArcTg (It is strong oscillatory)
void SolExactStrongOsc1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SolExactStrongOsc2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SolExactStrongOsc3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
// Exact solution of the differential equation  (Laplacian)
void SourceFunctionStrongOsc1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void SourceFunctionStrongOsc2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void SourceFunctionStrongOsc3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionStrongOsc1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionStrongOsc2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionStrongOsc3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);

// To exact solution in mixed formulation
void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out, TPZVec<STATE> &errorHDiv, void(*Exact)(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du));

void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &res);
void PermeabilityTensor(const TPZVec<REAL> &pt, TPZVec<STATE> &kabs, TPZFMatrix<STATE> &tensorK);
void ReactionTerm(const TPZVec<REAL> &pt, TPZVec<STATE> &alpha, TPZFMatrix<STATE> &disp);

// Choosing model or problem 
int ChooseEquation(int Equation, int nummethod, int Dim, TPZAutoPointer<TPZFunction<STATE> > &SourceFunc, TPZAutoPointer<TPZFunction<STATE> > &solExata);

#endif
