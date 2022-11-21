#include "TPZMaterial.h"
#include "TPZMatLaplacian.h"
#include "mixedpoisson.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"

#include "TPZCompElDisc.h"

#include "pzbuildmultiphysicsmesh.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"

#include "LaplacianSenoSeno.h"

#include "pzgmesh.h"
#include "IdMaterialsToGIDFiles.h"

// To create simple geometrical mesh
TPZGeoMesh* CreatingGeometricMesh(int elementIdMat, int bcdirichlet, int typeel, int &DimProblem) {
	TPZGeoMesh* gMesh = new TPZGeoMesh();

	if (typeel < 1 || typeel>15) {
		DimProblem = 0; delete gMesh; gMesh = 0; return NULL;
	}
	else if (typeel < 3) {
		DimProblem = 1;
		const int nnodes = 3;
		double coord[nnodes][1] = { {-1.},{0.},{1.} };
		for (int i = 0; i < nnodes; i++) {
			int nodind = gMesh->NodeVec().AllocateNewElement();
			TPZManVector<REAL, 3> nodeCoord(3, 0.);
			nodeCoord[0] = coord[i][0];
			nodeCoord[1] = 0.;
			nodeCoord[2] = 0.;
			gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
		}

		///Criando elementos
		const int nel = 2;
		int els[nel][2] = { {0,1},{1,2} };
		for (int iel = 0; iel < nel; iel++) {
			TPZManVector<int64_t, 2> nodind(2);
			int64_t index;
			nodind[0] = els[iel][0];
			nodind[1] = els[iel][1];
			gMesh->CreateGeoElement(EOned, nodind, elementIdMat, index);
		}

		///Criando elementos de contorno
		const int nelbc = 2;
		int bcels[nelbc][2] = { {0,bcdirichlet},{2,bcdirichlet} };
		for (int iel = 0; iel < nelbc; iel++) {
			TPZManVector<int64_t, 2> nodind(1);
			int64_t index;
			nodind[0] = bcels[iel][0];
			int matid = bcels[iel][1];
			gMesh->CreateGeoElement(EPoint, nodind, matid, index);
		}
	}
	else if (typeel < 8) {
		DimProblem = 2;
		if (typeel == 3) {
			///Criando nós
			const int nnodes = 9;
			double coord[nnodes][2] = { {-1,-1},{0,-1},{0,0},{-1,0}, {1.,-1.},{1,0}, {1,1},{0,1},{-1,1} };
			for (int i = 0; i < nnodes; i++) {
				int nodind = gMesh->NodeVec().AllocateNewElement();
				TPZManVector<REAL, 3> nodeCoord(3);
				nodeCoord[0] = coord[i][0];
				nodeCoord[1] = coord[i][1];
				nodeCoord[2] = 0.;
				gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
			}

			///Criando elementos
			const int nel = 8;
			int els[nel][3] = { { 0,1,2 },{0,2,3},{ 1,5,2 },{ 1,4,5 },{2,5,6},{2,6,7},{3,2,7},{3,7,8} };
			for (int iel = 0; iel < nel; iel++) {
				TPZManVector<int64_t, 3> nodind(3);
				int64_t index;
				nodind[0] = els[iel][0];
				nodind[1] = els[iel][1];
				nodind[2] = els[iel][2];
				gMesh->CreateGeoElement(ETriangle, nodind, elementIdMat, index);
			}

		}
		else {
			///Criando nós
			const int nnodes = 9;
			double coord[nnodes][2] = { {-1,-1},{0,-1},{0,0},{-1,0}, {1.,-1.},{1,0}, {1,1},{0,1},{-1,1} };
			for (int i = 0; i < nnodes; i++) {
				int nodind = gMesh->NodeVec().AllocateNewElement();
				TPZManVector<REAL, 3> nodeCoord(3);
				nodeCoord[0] = coord[i][0];
				nodeCoord[1] = coord[i][1];
				nodeCoord[2] = 0.;
				gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
			}

			///Criando elementos
			const int nel = 4;
			int els[nel][4] = { {0,1,2,3},{1,4,5,2},{3,2,7,8},{2,5,6,7} };
			for (int iel = 0; iel < nel; iel++) {
				TPZManVector<int64_t, 4> nodind(4);
				int64_t index;
				nodind[0] = els[iel][0];
				nodind[1] = els[iel][1];
				nodind[2] = els[iel][2];
				nodind[3] = els[iel][3];
				gMesh->CreateGeoElement(EQuadrilateral, nodind, elementIdMat, index);
			}
		}
		///Criando elementos de contorno
		const int nelbc = 8;
		int bcels[nelbc][3] = { {0,1,bcdirichlet},{1,4,bcdirichlet},{4,5,bcdirichlet},{5,6,bcdirichlet},{6,7,bcdirichlet},{7,8,bcdirichlet},{8,3,bcdirichlet},{3,0,bcdirichlet} };
		for (int iel = 0; iel < nelbc; iel++) {
			TPZManVector<int64_t, 4> nodind(2);
			int64_t index;
			nodind[0] = bcels[iel][0];
			nodind[1] = bcels[iel][1];
			int matid = bcels[iel][2];
			gMesh->CreateGeoElement(EOned, nodind, matid, index);
		}
	}
	else {
		DimProblem = 3;
		if (typeel == 8) {   // prism
			///Criando nós
			const int nnodes = 8;
			double coord[nnodes][3] = { {-1.,-1.,-1.},{-1.,1.,-1.},{-1.,1.,1.},{-1.,-1.,1.},{1.,-1.,-1.},{1.,1.,-1.},{1.,1.,1.},{1.,-1.,1.} };
			for (int i = 0; i < nnodes; i++) {
				int nodind = gMesh->NodeVec().AllocateNewElement();
				TPZManVector<REAL, 3> nodeCoord(3);
				nodeCoord[0] = coord[i][0];
				nodeCoord[1] = coord[i][1];
				nodeCoord[2] = coord[i][2];
				gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
			}

			///Criando elementos
			const int nel = 2;
			int els[nel][6] = { { 0,1,5,3,2,6 }, { 0,5,4,3,6,7 } };
			for (int iel = 0; iel < nel; iel++) {
				TPZManVector<int64_t, 6> nodind(6);
				int64_t index;
				nodind[0] = els[iel][0];
				nodind[1] = els[iel][1];
				nodind[2] = els[iel][2];
				nodind[3] = els[iel][3];
				nodind[4] = els[iel][4];
				nodind[5] = els[iel][5];
				gMesh->CreateGeoElement(EPrisma, nodind, elementIdMat, index);
			}
			///Criando elementos de contorno
			const int nelbc = 4;
			int iel;
			int bcels[nelbc][5] = { {0,1,2,3,bcdirichlet},{1,5,6,2,bcdirichlet},{4,5,6,7,bcdirichlet},{0,4,7,3,bcdirichlet} };
			for (iel = 0; iel < nelbc; iel++) {
				TPZManVector<int64_t, 4> nodind(4);
				int64_t index;
				nodind[0] = bcels[iel][0];
				nodind[1] = bcels[iel][1];
				nodind[2] = bcels[iel][2];
				nodind[3] = bcels[iel][3];
				int matid = bcels[iel][4];
				gMesh->CreateGeoElement(EQuadrilateral, nodind, matid, index);
			}
			int bcelst[nelbc][4] = { {0,1,5,bcdirichlet},{0,5,4,bcdirichlet},{3,2,6,bcdirichlet},{3,6,7,bcdirichlet} };
			for (iel = 0; iel < nelbc; iel++) {
				TPZManVector<int64_t, 3> nodind(3);
				int64_t index;
				nodind[0] = bcelst[iel][0];
				nodind[1] = bcelst[iel][1];
				nodind[2] = bcelst[iel][2];
				int matid = bcelst[iel][3];
				gMesh->CreateGeoElement(ETriangle, nodind, matid, index);
			}
		}
		else if (typeel == 10) {   // hexahedra
			///Criando nós
			const int nnodes = 8;
			double coord[nnodes][3] = { {-1.,-1.,-1.},{-1.,1.,-1.},{-1.,1.,1.},{-1.,-1.,1.},{1.,-1.,-1.},{1.,1.,-1.},{1.,1.,1.},{1.,-1.,1.} };
			for (int i = 0; i < nnodes; i++) {
				int nodind = gMesh->NodeVec().AllocateNewElement();
				TPZManVector<REAL, 3> nodeCoord(3);
				nodeCoord[0] = coord[i][0];
				nodeCoord[1] = coord[i][1];
				nodeCoord[2] = coord[i][2];
				gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
			}

			///Criando elementos
			const int nel = 1;
			int els[nel][8] = { { 0,1,2,3,4,5,6,7 } };
			for (int iel = 0; iel < nel; iel++) {
				TPZManVector<int64_t, 8> nodind(8);
				int64_t index;
				nodind[0] = els[iel][0];
				nodind[1] = els[iel][1];
				nodind[2] = els[iel][2];
				nodind[3] = els[iel][3];
				nodind[4] = els[iel][4];
				nodind[5] = els[iel][5];
				nodind[6] = els[iel][6];
				nodind[7] = els[iel][7];
				gMesh->CreateGeoElement(ECube, nodind, elementIdMat, index);
			}
			///Criando elementos de contorno
			const int nelbc = 6;
			int bcels[nelbc][5] = { {0,1,2,3,bcdirichlet},{0,1,5,4,bcdirichlet},{1,5,6,2,bcdirichlet},{4,5,6,7,bcdirichlet},{0,4,7,3,bcdirichlet},{3,2,6,7,bcdirichlet} };
			for (int iel = 0; iel < nelbc; iel++) {
				TPZManVector<int64_t, 4> nodind(4);
				int64_t index;
				nodind[0] = bcels[iel][0];
				nodind[1] = bcels[iel][1];
				nodind[2] = bcels[iel][2];
				nodind[3] = bcels[iel][3];
				int matid = bcels[iel][4];
				gMesh->CreateGeoElement(EQuadrilateral, nodind, matid, index);
			}
		}
		else if (typeel == 12) { // tetrahedra
	///Criando nós
			const int nnodes = 9;
			double coord[nnodes][3] = { {-1.,-1.,-1.},{-1.,1.,-1.},{-1.,1.,1.},{-1.,-1.,1.},{1.,-1.,-1.},{1.,1.,-1.},{1.,1.,1.},{1.,-1.,1.}, {0.,0.,0.} };
			for (int i = 0; i < nnodes; i++) {
				int nodind = gMesh->NodeVec().AllocateNewElement();
				TPZManVector<REAL, 3> nodeCoord(3);
				nodeCoord[0] = coord[i][0];
				nodeCoord[1] = coord[i][1];
				nodeCoord[2] = coord[i][2];
				gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
			}
			///Criando elementos
			const int nel = 12;
			int els[nel][4] = { {0,1,5,8},{0,5,4,8}, {3,2,6,8},{3,6,7,8}, {0,1,2,8},{0,2,3,8}, {4,5,6,8},{4,6,7,8}, {1,5,6,8},{1,6,2,8}, {4,0,3,8},{4,3,7,8} };
			for (int iel = 0; iel < nel; iel++) {
				TPZManVector<int64_t, 4> nodind(4);
				int64_t index;
				nodind[0] = els[iel][0];
				nodind[1] = els[iel][1];
				nodind[2] = els[iel][2];
				nodind[3] = els[iel][3];
				gMesh->CreateGeoElement(ETetraedro, nodind, elementIdMat, index);
			}
			///Criando elementos de contorno
			const int nelbc = 12;
			int bcels[nelbc][4] = { {0,1,5,bcdirichlet},{0,5,4,bcdirichlet},{3,2,6,bcdirichlet},{3,6,7,bcdirichlet},{0,1,2,bcdirichlet},{0,2,3,bcdirichlet}, {4,5,6,bcdirichlet},{4,6,7,bcdirichlet},{1,5,6,bcdirichlet},{1,6,2,bcdirichlet},{0,4,3,bcdirichlet},{3,4,7,bcdirichlet} };
			for (int iel = 0; iel < nelbc; iel++) {
				TPZManVector<int64_t, 3> nodind(3);
				int64_t index;
				nodind[0] = bcels[iel][0];
				nodind[1] = bcels[iel][1];
				nodind[2] = bcels[iel][2];
				int matid = bcels[iel][3];
				gMesh->CreateGeoElement(ETriangle, nodind, matid, index);
			}
		}
		else if (typeel == 14) {   // pyramid
			///Criando nós
			const int nnodes = 9;
			double coord[nnodes][3] = { {-1.,-1.,-1.},{-1.,1.,-1.},{-1.,1.,1.},{-1.,-1.,1.},{1.,-1.,-1.},{1.,1.,-1.},{1.,1.,1.},{1.,-1.,1.}, {0.,0.,0.} };
			for (int i = 0; i < nnodes; i++) {
				int nodind = gMesh->NodeVec().AllocateNewElement();
				TPZManVector<REAL, 3> nodeCoord(3);
				nodeCoord[0] = coord[i][0];
				nodeCoord[1] = coord[i][1];
				nodeCoord[2] = coord[i][2];
				gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
			}

			///Criando elementos
			const int nel = 6;
			int els[nel][5] = { {0,1,5,4,8}, {3,2,6,7,8}, {0,1,2,3,8}, {4,5,6,7,8}, {1,5,6,2,8}, {4,0,3,7,8} };
			for (int iel = 0; iel < nel; iel++) {
				TPZManVector<int64_t, 5> nodind(5);
				int64_t index;
				nodind[0] = els[iel][0];
				nodind[1] = els[iel][1];
				nodind[2] = els[iel][2];
				nodind[3] = els[iel][3];
				nodind[4] = els[iel][4];
				gMesh->CreateGeoElement(EPiramide, nodind, elementIdMat, index);
			}
			///Criando elementos de contorno
			const int nelbc = 6;
			int bcels[nelbc][5] = { {0,1,2,3,bcdirichlet},{1,5,6,2,bcdirichlet},{4,5,6,7,bcdirichlet},{4,0,3,7,bcdirichlet},{3,2,6,7,bcdirichlet},{0,1,5,4,bcdirichlet} };
			for (int iel = 0; iel < nelbc; iel++) {
				TPZManVector<int64_t, 4> nodind(4);
				int64_t index;
				nodind[0] = bcels[iel][0];
				nodind[1] = bcels[iel][1];
				nodind[2] = bcels[iel][2];
				nodind[3] = bcels[iel][3];
				int matid = bcels[iel][4];
				gMesh->CreateGeoElement(EQuadrilateral, nodind, matid, index);
			}
		}
		else
			gMesh = NULL;
	}
	///Construindo conectividade da malha
	gMesh->SetDimension(DimProblem);
	gMesh->BuildConnectivity();
	return gMesh;
}

TPZGeoMesh *CreatingGeometricMeshFromGIDFile(std::string &GeoGridFile, int ElementIDMat, int bcdirichlet, int typeel, int &DimProblem, int &nref) {
	//nref = 1;
	if (typeel < 1 || typeel>15) {
		DimProblem = 0; return 0;
	}
	else if (typeel < 3) DimProblem = 1;
	else if (typeel < 8) DimProblem = 2;
	else DimProblem = 3;

	// IMPORTANDO MALHA GEOMETRICA DESDE ARQUIVO GID
	TPZReadGIDGrid myreader;
	TPZGeoMesh* gmesh;

	if (GeoGridFile.empty())
		return 0;
	/* Note: the generator from file gid creates new geometric mesh, then is necessary to delete and assert the dimension of the new geometrical mesh */
	gmesh = myreader.GeometricGIDMesh(GeoGridFile);
	if (!gmesh) {
		std::cout << "\nGid file is not found.\n";
		return NULL;
	}
	gmesh->SetDimension(DimProblem);

	// Inserting material id
	SetMaterialIdForElements(gmesh, ElementIDMat);

	// With only one boundary condition (Dirichlet with value ZERO.
	SetMatBCIdForElementsSameByHeigh(gmesh, bcdirichlet, 'a');
	return gmesh;

	// EXAMPLE TO USE INITIAL AND LAST POINTS TO DETERMINE BC
	TPZManVector<REAL> P0(3, 0.);
	TPZManVector<REAL> P1(3, 0.);

	// BC IdMat in Bottom
	SetMatBCIdForElementsSameByHeigh(gmesh, -3, 't');
	SetMatBCIdForElementsSameByHeigh(gmesh, -3, 'b');
	SetMatBCIdForElementsSameByHeigh(gmesh, -3, 'v');
	P1[0] = 2.5;
	P0 = P1;
	P1[0] = 5.;
	SetMatBCIdForElementsSameByHeigh(gmesh, -3, 'r');
	SetMatBCIdForElementsSameByHeigh(gmesh, -3, 'f');
	SetMatBCIdForElementsSameByHeigh(gmesh, -3, 'l');

	return gmesh;
}

// Choosing model or problem 
int ChooseEquation(int Equation, int nummethod, int Dim, TPZAutoPointer<TPZFunction<STATE> > &SourceFunc, TPZAutoPointer<TPZFunction<STATE> > &solExata) {
    int nfunc = -1;
	TPZDummyFunction<STATE> *dums=0, *dum=0;
	if (!Equation) {
		if (Dim == 3) {
			if(nummethod==2)
				dums = new TPZDummyFunction<STATE>(MinusSourceFunctionSin3D, 5);
			else
				dums = new TPZDummyFunction<STATE>(SourceFunctionSin3D, 5);
			dum = new TPZDummyFunction<STATE>(SolExactSenoSenoSeno, 5);
            nfunc = 3;
		}
		else if (Dim == 2) {
			if (nummethod == 2)
				dums = new TPZDummyFunction<STATE>(MinusSourceFunctionSin2D, 5);
			else
				dums = new TPZDummyFunction<STATE>(SourceFunctionSin2D, 5);
			dum = new TPZDummyFunction<STATE>(SolExactSenoSeno, 5);
            nfunc = 2;
		}
		else {
			if (nummethod == 2)
				dums = new TPZDummyFunction<STATE>(MinusSourceFunctionSin1D, 5);
			else
				dums = new TPZDummyFunction<STATE>(SourceFunctionSin1D, 5);
			dum = new TPZDummyFunction<STATE>(SolExactSeno, 5);
            nfunc = 1;
		}
	}
	else if (Equation == 1) {
		if (Dim == 3) {
			if (nummethod == 2)
				dums = new TPZDummyFunction<STATE>(MinusSourceFunctionArcTg3D, 5);
			else
				dums = new TPZDummyFunction<STATE>(SourceFunctionArcTg3D, 5);
			dum = new TPZDummyFunction<STATE>(SolExactArcTg3D, 5);
            nfunc = 6;
		}
		else if (Dim == 2) {
			if (nummethod == 2)
				dums = new TPZDummyFunction<STATE>(MinusSourceFunctionArcTg2D, 5);
			else
				dums = new TPZDummyFunction<STATE>(SourceFunctionArcTg2D, 5);
			dum = new TPZDummyFunction<STATE>(SolExactArcTg2D, 5);
            nfunc = 5;
		}
		else {
			if (nummethod == 2)
				dums = new TPZDummyFunction<STATE>(MinusSourceFunctionArcTg1D, 5);
			else
				dums = new TPZDummyFunction<STATE>(SourceFunctionArcTg1D, 5);
			dum = new TPZDummyFunction<STATE>(SolExactArcTg1D, 5);
            nfunc = 4;
		}
	}
	else if (Equation == 2) {
		if (Dim == 3) {
			if (nummethod == 2)
				dums = new TPZDummyFunction<STATE>(MinusSourceFunctionStrongOsc3D, 5);
			else
				dums = new TPZDummyFunction<STATE>(SourceFunctionStrongOsc3D, 5);
			dum = new TPZDummyFunction<STATE>(SolExactStrongOsc3D, 5);
            nfunc = 9;
		}
		else if (Dim == 2) {
			if (nummethod == 2)
				dums = new TPZDummyFunction<STATE>(MinusSourceFunctionStrongOsc2D, 5);
			else
				dums = new TPZDummyFunction<STATE>(SourceFunctionStrongOsc2D, 5);
			dum = new TPZDummyFunction<STATE>(SolExactStrongOsc2D, 5);
            nfunc = 8;
		}
		else {
			if (nummethod == 2)
				dums = new TPZDummyFunction<STATE>(MinusSourceFunctionStrongOsc1D, 5);
			else
				dums = new TPZDummyFunction<STATE>(SourceFunctionStrongOsc1D, 5);
			dum = new TPZDummyFunction<STATE>(SolExactStrongOsc1D, 5);
            nfunc = 7;
		}
	}
	else {
		std::cout << "\nEquation or Problem undefined.\n\n";
//		return false;
	}
	SourceFunc = dums;
	solExata = dum;
	return nfunc;
}
// To computational mesh
void CreatingComputationalMesh(TPZCompMesh* cMesh, int elementIdMat, int bcdirichlet, bool DGFEM,int DimProblem,int pOrder, TPZAutoPointer<TPZFunction<STATE> > &SourceFunc, TPZAutoPointer<TPZFunction<STATE> > &solExata) {
	/// criar materiais
	TPZMatPoisson3d *mat;
	mat = new TPZMatPoisson3d(elementIdMat, DimProblem);
//	TPZMatLaplacian* mat = new TPZMatLaplacian(elementIdMat, DimProblem);
//	const STATE K = 1., F = 0.;
//	mat->SetParameters(K, F);
	mat->SetForcingFunction(SourceFunc);
	mat->SetForcingFunctionExact(solExata);

	if (DGFEM) {
		///Formulacao nao-simetria de Baumann, Oden e Babuska sem penalizacao
		mat->SetNoPenalty();
		mat->SetNonSymmetric();
	}
	cMesh->InsertMaterialObject(mat);

	///Condições de contorno
	//mat->SetForcingFunction();
	TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);
	TPZBndCond* BCondDirichletNulo = mat->CreateBC(mat, bcdirichlet, 0, val1, val2);//0 = Dirichlet
	cMesh->InsertMaterialObject(BCondDirichletNulo);


	TPZBndCond* BCondNeumannZero = mat->CreateBC(mat, -2, 1, val1, val2);//1 = Neumann
	cMesh->InsertMaterialObject(BCondNeumannZero);

	// Para tipo de problema Problem = 1
	val2(0, 0) = 2.0; // val2(1, 0) = 0.5;
	TPZMaterial* BCondNeumannLeft = mat->CreateBC(mat, -7, 1, val1, val2);//1 = Neumann
	cMesh->InsertMaterialObject(BCondNeumannLeft);
	val1(0, 0) = 1.;
	TPZMaterial* BCondNeumann = mat->CreateBC(mat, -8, 1, val1, val2);//1 = Neumann
	cMesh->InsertMaterialObject(BCondNeumann);

	cMesh->SetDefaultOrder(pOrder);
	cMesh->SetDimModel(DimProblem);

	if (DGFEM) {
		///Criando malha de Galerkin descontínuo
		cMesh->SetAllCreateFunctionsDiscontinuous();
	}
	else {
		///cria malha H1
		cMesh->SetAllCreateFunctionsContinuous();
	}

	///Criando elementos computacionais
	cMesh->AutoBuild();

	if (DGFEM) {
		///Cria elementos de interface
		TPZCreateApproximationSpace::CreateInterfaces(*cMesh);
	}
}
void CMeshPressure(TPZCompMesh *cmesh, int elementIdMat, bool DGFEM, int DimProblem, int pOrder)
{
	/// criar materiais
	TPZMixedPoisson *mat;
	mat = new TPZMixedPoisson(elementIdMat, DimProblem);

	 if (DGFEM) {
		 ///Formulacao nao-simetria de Baumann, Oden e Babuska sem penalizacao
		 mat->SetNoPenalty();
		 mat->SetNonSymmetric();
	 }
	 cmesh->InsertMaterialObject(mat);

	cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(DimProblem);

	if(DGFEM)
		cmesh->SetAllCreateFunctionsDiscontinuous();
	else 
		cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();

	cmesh->ApproxSpace().CreateDisconnectedElements(true);

	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
	if (DGFEM) {
		///Cria elementos de interface
		TPZCreateApproximationSpace::CreateInterfaces(*cmesh);
	}

	int64_t n_connects = cmesh->NConnects();
	for (int64_t i = 0; i < n_connects; ++i) {
		cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
	}
}

#include "TPZVecL2.h"
 void CMeshFlux(TPZCompMesh* cmesh, int elementIdMat, int DimProblem, int pOrder, TPZAutoPointer<TPZFunction<STATE> > &solExata) {
	 TPZVecL2 *mat;
	mat = new TPZVecL2(elementIdMat);
	mat->SetDimension(DimProblem);

	cmesh->InsertMaterialObject(mat);
	TPZFNMatrix<1, STATE> val1(1, 1, 0.), val2(1, 1, 1.);
	TPZBndCond *bc = mat->CreateBC(mat, -3, 0, val1, val2);
	bc->TPZMaterial::SetForcingFunction(solExata);
	cmesh->InsertMaterialObject(bc);

	cmesh->SetDefaultOrder(pOrder);
	cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(DimProblem);

	cmesh->SetDimModel(DimProblem);
	cmesh->AutoBuild();//Ajuste da estrutura de dados computacional
	cmesh->InitializeBlock();
}
#include "pzintel.h"
#include "mixedpoisson.h"
#include "TPZCompMeshTools.h"

void MalhaCompMultifisica(TPZMultiphysicsCompMesh* mphysics, int ElementIDMat, int bcdirichlet, bool DGFEM, int DimProblem, int pOrder, TPZAutoPointer<TPZFunction<STATE> > &source, TPZAutoPointer<TPZFunction<STATE> > &solExata) {

	TPZMixedPoisson *material = new TPZMixedPoisson(ElementIDMat, DimProblem);

	//permeabilidade
	TPZFMatrix<REAL> Ktensor(DimProblem+1, DimProblem+1, 0.);
	TPZFMatrix<REAL> InvK(DimProblem+1, DimProblem+1, 0.);
	Ktensor.Identity();
	InvK.Identity();

//	TPZAutoPointer<TPZFunction<STATE> > solexata;
//	solexata = new TPZDummyFunction<STATE>(SolExactProblem, 5);
	material->SetForcingFunctionExact(solExata);
	material->SetPermeability(1.);

	int int_order = 10;
	//funcao do lado direito da equacao do problema
	material->SetForcingFunction(source);

	material->SetPermeabilityTensor(Ktensor, InvK);

	//inserindo o material na malha computacional
	TPZMaterial *mat(material);
	mphysics->InsertMaterialObject(mat);
	mphysics->SetDimModel(DimProblem);

	//Criando condicoes de contorno
	TPZFMatrix<STATE> val1(1, 1, 0.), val2(1, 1, 0.);

	TPZMaterial * BCond0 = material->CreateBC(mat, bcdirichlet, 0, val1, val2);

	///Inserir condicoes de contorno
	mphysics->InsertMaterialObject(BCond0);

	mphysics->SetAllCreateFunctionsMultiphysicElem();

	TPZManVector<int> active(2, 1);
	TPZManVector<TPZCompMesh *> meshvector(2, 0);

	TPZGeoMesh *gmesh = mphysics->Reference();
	meshvector[0] = new TPZCompMesh(gmesh);
	gmesh->ResetReference();
	CMeshFlux(meshvector[0], ElementIDMat, DimProblem, pOrder,solExata);

	meshvector[1] = new TPZCompMesh(gmesh);
//	meshvector[0] = meshvec[0]; //CreateFluxHDivMesh(problem);
	CMeshPressure(meshvector[1], ElementIDMat, DGFEM, DimProblem, pOrder);
//	meshvector[1] = meshvec[1]; // CreatePressureMesh(problem);
	mphysics->BuildMultiphysicsSpace(active, meshvector);
	mphysics->LoadReferences();
	bool keepmatrix = false;
	bool keeponelagrangian = true;
	TPZCompMeshTools::CreatedCondensedElements(mphysics, keeponelagrangian, keepmatrix);
}


// To Exact solution of the differential equation for three formulations
// MODEL WITH SIN*SIN*SIN SOLUTION --> EQUATION = 0
// 1D - 2D - 3D
void SolExactSeno(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	u[0] = sin(M_PI*x);
	du(0, 0) = M_PI * cos(M_PI*x);
//	du(1, 0) = -M_PI*M_PI*sin(M_PI*x);
}
void SourceFunctionSin1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	u.Fill(0.);
	u[0] = -M_PI*M_PI*sin(M_PI*x);
}
void MinusSourceFunctionSin1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionSin1D(loc, u);
	u[0] *= -1.;
}
void SolExactSenoSeno(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	u[0] = sin(M_PI*x)*sin(M_PI*y);
	du(0, 0) = M_PI * cos(M_PI*x)*sin(M_PI*y);
	du(1, 0) = M_PI * cos(M_PI*y)*sin(M_PI*x);
//	du(2, 0) = -2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
}
void SourceFunctionSin2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	u[0] = -2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
}
void MinusSourceFunctionSin2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionSin2D(loc, u);
	u[0] *= -1.;
}
void SolExactSenoSenoSeno(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL z = loc[2];
	u[0] = sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
	du(0, 0) = M_PI * cos(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
	du(1, 0) = M_PI * cos(M_PI*y)*sin(M_PI*x)*sin(M_PI*z);
	du(2, 0) = M_PI * cos(M_PI*z)*sin(M_PI*x)*sin(M_PI*y);
//	du(3, 0) = -3 * M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
}
void SourceFunctionSin3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL z = loc[2];
	u[0] = -3*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
}
void MinusSourceFunctionSin3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionSin3D(loc, u);
	u[0] *= -1.;
}
// MODEL WITH ARCTG SOLUTION --> EQUATION = 1
// 1D - 2D - 3D
REAL CoefSol = 5.;
REAL CoefArgument = 2.;
REAL ArgumentSum = 1.;
REAL CoefRadio = 2.;
REAL CoefPi = 5.;
REAL ValorY = 0.;

void SolExactArcTg1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL r = x * x;
	REAL w = CoefArgument *(ArgumentSum - CoefRadio *r);
	REAL v = 1. + w * w;
	u[0] = CoefSol *(0.5*M_PI + atan(w));
	du(0, 0) = -(2.*CoefArgument*CoefRadio*CoefSol*x)/v;
//	du(1, 0) = -(200./v)*(1.+(80*r*w/v));
}
void SourceFunctionArcTg1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL r = x * x;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL v = 1. + w * w;
	REAL den = CoefRadio / v;
	u[0] = -2.*CoefArgument*CoefSol*den*(1. + 4.*CoefArgument*den*w*r);
}
void MinusSourceFunctionArcTg1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionArcTg1D(loc, u);
	u[0] *= -1.;
}
void SolExactArcTg2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL r = x * x + y * y;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL v = 1. + w * w;
	u[0] = CoefSol * (0.5*M_PI + atan(w));
	du(0, 0) = -(2.*CoefArgument*CoefRadio*CoefSol*x) / v;
	du(1, 0) = -(2.*CoefArgument*CoefRadio*CoefSol*y) / v;
}
void SourceFunctionArcTg2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL r = x * x + y * y;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL v = 1. + w * w;
	REAL den = CoefRadio / v;
	u[0] = -2.*CoefArgument*CoefSol*den*(2. + 4.*CoefArgument*den*w*r);   // 2. = dim because all second derivatives has this somand
}
void MinusSourceFunctionArcTg2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionArcTg2D(loc, u);
	u[0] *= -1.;
}
void SolExactArcTg3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL z = loc[2];
	const REAL r = x * x + y * y + z * z;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL v = 1. + w * w;
	u[0] = CoefSol * (0.5*M_PI + atan(w));
	du(0, 0) = -(2.*CoefArgument*CoefRadio*CoefSol*x) / v;
	du(1, 0) = -(2.*CoefArgument*CoefRadio*CoefSol*y) / v;
	du(2, 0) = -(2.*CoefArgument*CoefRadio*CoefSol*z) / v;
}
void SourceFunctionArcTg3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL z = loc[2];
	const REAL r = x * x + y * y + z * z;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL v = 1. + w * w;
	REAL den = CoefRadio / v;
	u[0] = -2.*CoefArgument*CoefSol*den*(3. + 4.*CoefArgument*den*w*r);   // 3. = dim because all second derivatives has this somand
}
void MinusSourceFunctionArcTg3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionArcTg3D(loc, u);
	u[0] *= -1.;
}
// MODEL WITH SIN*COS*ARCTG SOLUTION --> EQUATION = 2
// 1D - 2D - 3D
void SolExactStrongOsc1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL r = x * x;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL Factor3 = 0.5*M_PI + atan(w);
	REAL v = 1 + w * w;
	REAL den = CoefRadio / v;
	REAL angle = CoefPi * M_PI;
	REAL YY = 1. + cos(angle*ValorY);
	u[0] = YY*sin(angle*x)*Factor3;
	du(0, 0) = M_PI*YY*CoefPi*cos(angle*x)*Factor3-(2*CoefArgument*YY*x*den*sin(angle*x));
}
void SourceFunctionStrongOsc1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL r = x * x;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL Factor3 = 0.5*M_PI + atan(w);
	REAL v = 1 + w * w;
	REAL den = CoefRadio / v;
	REAL angle = CoefPi * M_PI;
	REAL YY = 1. + cos(angle*ValorY);
	u[0] = -M_PI*M_PI*CoefPi*CoefPi*YY*sin(angle*x)*Factor3-(2*CoefArgument*YY*sin(angle*x)*den)-(8*CoefArgument*den*den*YY*r*w*w*sin(angle*x))-(4*M_PI*CoefArgument*CoefPi*den*YY*x*cos(angle*x));
}
void MinusSourceFunctionStrongOsc1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionStrongOsc1D(loc, u);
	u[0] *= -1.;
}
void SolExactStrongOsc2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL r = x * x + y * y;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL Factor3 = 0.5*M_PI + atan(w);
	REAL v = 1 + w * w;
	REAL den = CoefRadio / v;
	REAL angle = CoefPi * M_PI;
	REAL CosY = 1. + cos(angle*y);
	u[0] = sin(angle*x)*CosY*Factor3;
	du(0, 0) = CosY*(M_PI * CoefPi*cos(angle*x)*Factor3 - (2 * CoefArgument*x*den*sin(angle*x)));
	du(1, 0) = (-1.)*sin(angle*x)*(M_PI * CoefPi*sin(angle*y)*Factor3 + (2 * CoefArgument*y*den*CosY));
}
void SourceFunctionStrongOsc2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL r = x * x + y * y;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL Factor3 = 0.5*M_PI + atan(w);
	REAL v = 1 + w * w;
	REAL den = CoefRadio / v;
	REAL angle = CoefPi * M_PI;
	REAL CosY = 1. + cos(angle*y);
	u[0] = (-1.)*(M_PI*M_PI*CoefPi*CoefPi*sin(angle*x)*Factor3*(CosY+cos(angle*y)));
	u[0] += 4.*den*CoefArgument*(angle*(y*sin(angle*x)*sin(angle*y)-(x*CosY*cos(angle*x)))-(CosY*sin(angle*x)));
	u[0] += (-8.) * CoefArgument*CoefArgument*den*den*w*sin(angle*x)*CosY*r;
}
void MinusSourceFunctionStrongOsc2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionStrongOsc2D(loc, u);
	u[0] *= -1.;
}
void SolExactStrongOsc3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL z = loc[2];
	const REAL r = x * x + y * y + z * z;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL Factor3 = 0.5*M_PI + atan(w);
	REAL v = 1 + w * w;
	REAL den = CoefRadio / v;
	REAL angle = CoefPi * M_PI;
	REAL CosY = 1. + cos(angle*y);
	REAL CosZ = 1. + cos(angle*z);
	u[0] = sin(angle*x)*CosY*CosZ*Factor3;
	du(0, 0) = CosZ*CosY * (M_PI * CoefPi*cos(angle*x)*Factor3 - (2 * CoefArgument*x*den*sin(angle*x)));
	du(1, 0) = (-1.)*CosZ*sin(angle*x)*(M_PI * CoefPi*sin(angle*y)*Factor3 + (2 * CoefArgument*y*den*CosY));
	du(2, 0) = (-1.)*CosY*sin(angle*x)*(M_PI * CoefPi*sin(angle*z)*Factor3 + (2 * CoefArgument*z*den*CosZ));
}
void SourceFunctionStrongOsc3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL z = loc[2];
	const REAL r = x * x + y * y + z * z;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL Factor3 = 0.5*M_PI + atan(w);
	REAL v = 1 + w * w;
	REAL den = CoefRadio / v;
	REAL angle = CoefPi * M_PI;
	REAL CosY = 1. + cos(angle*y);
	REAL CosZ = 1. + cos(angle*z);
	u[0] = (-1.)*(M_PI*M_PI*CoefPi*CoefPi*sin(angle*x)*Factor3*(CosZ*CosY + CosZ*cos(angle*y) + CosY*cos(angle*z)));
	u[0] += 4.*M_PI*den*CoefArgument*CoefPi*(sin(angle*x)*CosY*z*sin(angle*z)+ sin(angle*x)*CosZ*y*sin(angle*y) - cos(angle*x)*CosY*x*CosZ);
	u[0] += (-6.)*CoefArgument*den*sin(angle*x)*CosY*CosZ;
	u[0] += (-8.)*CoefArgument*CoefArgument*den*den*w*sin(angle*x)*CosY*r;
}
void MinusSourceFunctionStrongOsc3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionStrongOsc3D(loc, u);
	u[0] *= -1.;
}


// Hdiv ERRORS
void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out, TPZVec<STATE> &errorHDiv, void(*Exact)(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du))
{
	int64_t nel = hdivmesh->NElements();
	int dim = hdivmesh->Dimension();
	TPZManVector<STATE, 10> globerrors(10, 0.);
	TPZStack<REAL> vech;

	for (int64_t el = 0; el < nel; el++) {
		TPZCompEl *cel = hdivmesh->ElementVec()[el];
		if (!cel) {
			continue;
		}

		TPZGeoEl *gel = cel->Reference();
		if (!gel || gel->Dimension() != dim) {
			continue;
		}
		TPZManVector<REAL, 10> elerror(10, 0.);
		cel->EvaluateError(Exact, elerror, 0);
		int nerr = elerror.size();
		for (int i = 0; i < nerr; i++) {
			globerrors[i] += elerror[i] * elerror[i];
		}

	}
	out << "L2 Error Norm for flux = " << sqrt(globerrors[1]) << std::endl;
	errorHDiv.Resize(3, 0.);
	errorHDiv[0] = sqrt(globerrors[1]);
	errorHDiv[1] = sqrt(globerrors[2]);
	errorHDiv[2] = sqrt(globerrors[3]);
	out << "L2 Norm for divergence = " << sqrt(globerrors[2]) << std::endl;
	out << "Hdiv Norm for flux = " << sqrt(globerrors[3]) << std::endl;

}
void PermeabilityTensor(const TPZVec<REAL> &pt, TPZVec<STATE> &kabs, TPZFMatrix<STATE> &tensorK)
{

	REAL x = pt[0];
	//REAL y = pt[1];
	tensorK.Resize(4, 2);
	kabs.Resize(1, 0.);

	//K
	//REAL temp = 1. + 10.*x;
	REAL temp = 1.;
	tensorK(0, 0) = temp;     tensorK(0, 1) = 0.0;
	tensorK(1, 0) = 0.0;      tensorK(1, 1) = temp;

	//Kinv
	tensorK(2, 0) = 1.0 / temp;     tensorK(2, 1) = 0.0;
	tensorK(3, 0) = 0.0;          tensorK(3, 1) = 1.0 / temp;
}

void ReactionTerm(const TPZVec<REAL> &pt, TPZVec<STATE> &alpha, TPZFMatrix<STATE> &disp)
{
	REAL x = pt[0];
	REAL y = pt[1];
	disp.Resize(2, 2);
	disp.Zero();
	alpha.Resize(1, 1.);

	//Termo de reacao
//	REAL temp = 1. - x * x - y * y;
//	alpha = exp(temp);
}
void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &res) {


	double x = pt[0];
	double y = pt[1];
	res[0] = 0.;

	REAL solp = sin(M_PI*x)*sin(M_PI*y);

	REAL temp1 = 1. + 10.*x;
	REAL temp2 = 1. - x * x - y * y;
	REAL temp3 = exp(temp2);

	res[0] = -10.*M_PI*cos(M_PI*x)*sin(M_PI*y) + (temp3 + 2.*M_PI*M_PI*temp1)*solp;
}

/*
// To Exact Solution in mixed formulation
void SolProblema(const TPZVec<REAL> &pt, TPZVec<STATE> &u, TPZFMatrix<STATE> &flux) {

	REAL x = pt[0];
	REAL y = pt[1];
	REAL z = pt[2];

	u.Resize(1, 0.);
	flux.Resize(3, 1);
	flux(0, 0) = 0., flux(1, 0) = 0., flux(2, 0) = 0.;

	//Solucao u
	REAL solp = sin(M_PI*x)*sin(M_PI*y);
	u[0] = solp;

	REAL temp = 1.;// +10.*x;

	//fluxo em x
	flux(0, 0) = M_PI * temp*cos(M_PI*x)*sin(M_PI*y);
	flux(0, 0) *= -1.;

	//fluxo em y
	flux(1, 0) = M_PI * temp*cos(M_PI*y)*sin(M_PI*x);
	flux(1, 0) *= -1.;

	//divergente: -(dux/dx + duy/dy)
	flux(2, 0) = -10.*M_PI*cos(M_PI*x)*sin(M_PI*y) + 2.*M_PI*M_PI*temp*sin(M_PI*x)*sin(M_PI*y);
}
*/
