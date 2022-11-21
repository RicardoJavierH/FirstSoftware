/// Laplacian example, with solution as Steklov function.
/// On 1D, 2D, 3D
/// Problem = 0 represent boundary condition with flux as Steklov function.
/// Problem = 1 has flux on left side and wall condition on the others faces

#include <iostream>
#include <math.h>
#include <time.h>
#include <stdio.h>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"
#include "tpzgeoelrefpattern.h"

#include "pzanalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"

#include "pzlog.h"

#include "LaplacianSenoSeno.h"
#include "IdMaterialsToGIDFiles.h"
#include "CreateGeoMesh.h"
#include "Refinements.h"
#include "commonJC.h"

// extern variables
std::string Partition;

// To Identifier of the material into the domain
int ElementIDMat = 5;
int bcdirichlet = -3;

int DimProblem;
//int Equation = 0;
// 0 - to use FEM with continuous functions, 1 - to use Galerkin discontinuous, 2 - to use mixed formulation
//int NumericMethod;
bool DGFEM = false;
int stepsolver = 0;  // 0 - direct (LU), (metodos iterativos: 1 - GMRES, 2 - CG (simetrico), ...

// Autopointer to source function and exactsolution of the problem
TPZAutoPointer<TPZFunction<STATE> > solExata;
TPZAutoPointer<TPZFunction<STATE> > SourceFunc;
void (*fExact)(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) = 0;

void PutExact(int n);

// PROGRAM
int main(int argc, char* argv[]) {

	std::string sout;
	sout = argv[0];
	PutNameDirectoryForResults(sout);

	/// PRIORITARY -> RELATIVE TO GEOMETRICAL MESH
	// Type of elements
	int typeel;		// 1 - linear(linear) 3 - triangle(line) 5 - quadrilateral(lin) 8 - prism(lin) 10 - hexahedra(lin) 12 - tetrahedra(lin) 14 - pyramid(lin)
	// To create geometric mesh from GID file or not
	std::string fromgid;

	// Pointers to meshes
	TPZGeoMesh* gmesh = 0;
	TPZMultiphysicsCompMesh *mphysics = 0;
	TPZCompMesh *cmesh = 0;

	// Parameter to refining mesh. Whether the mesh file is from GID the nref is made zero generally
	SimulationData simulationdata;
	// Importing user data from filedata.
	Partition = argv[0];
	bool fimport = UserDataImport(Partition,fromgid, simulationdata);
	if (!fimport) {
		cout << "\nUser data is not found.\nThe executable need of the \"LCC_Simulator\\userdata.txt\" file.\n\n";
		return 100;
	}
	if (!CheckTypeEl(simulationdata.typeel)) {
		cout << "\nType of element undefined.\n";
		return 200;
	}

	// To load refinement pattern
	if(simulationdata.Equation < 3 && simulationdata.typeel != 20) LoadUniformPatterns(simulationdata.typeel);
	else gRefDBase.InitializeAllUniformRefPatterns();

	// Using different approximation spaces 0 - FEM classic; 1 - DG; 2 - HdivH1; 3 (All one to one)
	int InitialEquation = 0;
	if (simulationdata.Equation != 3) {
		InitialEquation = simulationdata.Equation;
		simulationdata.Equation = InitialEquation + 1;
	}
	// To print information of processing times, system dimension and calculated errors (H1, GD and HdivH1)
	std::ofstream saida(sout + "InfoRun_Errors.txt", ios::app);
	if (!saida.is_open()) return -3;
	clock_t texec = clock();

	for (int modelo = InitialEquation; modelo < simulationdata.Equation; modelo++) {
		clock_t tmodelo = clock();
		// identifyer of the element types
		int elementtypes[8] = { 1, 3, 5, 8, 10, 12, 14, 20 };
		int nelementtypes = 0;
		if(simulationdata.typeel == 30)
			nelementtypes = 3;
		typeel = elementtypes[nelementtypes];
		if (simulationdata.typeel < 20) {
			typeel = simulationdata.typeel;
			nelementtypes = 8;
		}
		// Iteration over all element types (linear)
		do {
			// To control of the processing time
			clock_t t1 = clock();             // Get current time
			clock_t ttypeel = clock();
			if (nelementtypes != 8) typeel = elementtypes[nelementtypes++];

			/** Creating geometric mesh from GID file or not */
			if (fromgid.empty())
				gmesh = CreatingGeometricMesh(ElementIDMat, bcdirichlet, typeel, DimProblem);
			else
				gmesh = CreatingGeometricMeshFromGIDFile(fromgid, ElementIDMat, bcdirichlet, typeel, DimProblem, simulationdata.nref);
			if (!gmesh || !gmesh->NElements()) return -2;
			saida << "\nMODEL " << modelo << "\t Dimension " << DimProblem << "\t Element Type " << typeel << "\n\n";
			std::cout << "\nMODEL " << modelo << std::endl;
			std::cout << " Dimension " << DimProblem << "\t Element Type " << typeel << "\t Level of Refinement " << simulationdata.nref << std::endl;

			// Refining mesh (uniform)
			UniformRefinement(simulationdata.nref, gmesh);
			// Printing geomesh only
			if (simulationdata.print) {
				std::stringstream gmeshfile;
				gmeshfile << sout << "GMesh_E" << typeel << "_Dim" << DimProblem << "D.vtk";
				std::ofstream fGMesh(gmeshfile.str());
				TPZVTKGeoMesh::PrintGMeshVTK(gmesh, fGMesh);
			}

			// end of time dedicated to geometric mesh processed
			clock_t t2 = clock();
			saida << "Process GMesh: \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2-t1 << std::endl;

			// Analysis
			TPZAnalysis* an = 0;
			TPZStepSolver<STATE> *step = 0;
			TPZManVector<REAL> erro;

			// Using different approximation spaces 0 - FEM classic; 1 - DG; 2 - HdivH1; 3 (All one to one)
			int InitialMethod = 0;
			if (simulationdata.NumericMethod != 3) {
				InitialMethod = simulationdata.NumericMethod;
				simulationdata.NumericMethod = InitialMethod + 1;
			}

			for (int numethod = InitialMethod; numethod < simulationdata.NumericMethod; numethod++) {
				clock_t tmethod = clock();
				t1 = clock();
				// Determining Equation problem and relationated exact solution (if exists)
                int nfunc = ChooseEquation(modelo, numethod, DimProblem, SourceFunc, solExata);
                PutExact(nfunc);
     //           fExact = solExata.operator TPZFunction<double> &();
				if (DimProblem == 1 && numethod == 1) continue;
				saida << "  Solving with Method " << numethod << "\n\n";
				std::cout << "  Solving with Method " << numethod << "\t Linear Solver " << stepsolver << std::endl;
				DGFEM = false;
				erro.Resize(0);

				//rodar formulacao mista
				if (numethod == 2 && typeel < 12) {
					saida << "Saida do erro para formulacão hdiv, com ordem p = " << simulationdata.pOrder << "\n";
					if (!fromgid.empty())
						saida << "Malha irregular from gid file.\n";

					// Criando a malha computacional multifísica
					TPZManVector<TPZCompMesh *, 2> meshvec(2, 0);
					//Creating computational mesh for multiphysic elements
					mphysics = new TPZMultiphysicsCompMesh(gmesh);
					MalhaCompMultifisica(mphysics, ElementIDMat, bcdirichlet, DGFEM, DimProblem, simulationdata.pOrder, SourceFunc, solExata);
					mphysics->InitializeBlock();
					meshvec = mphysics->MeshVector();

					int nDofTotal;
					nDofTotal = meshvec[0]->NEquations() + meshvec[1]->NEquations();

					int nDofCondensed;
					nDofCondensed = mphysics->NEquations();

					saida << "Grau de Liberdade Total = " << nDofTotal << std::endl << "Grau de Liberdade Condensado = " << nDofCondensed << std::endl;
					std::cout << "   Computational mesh created:\n" << "   Grau de Liberdade Total = " << nDofTotal << std::endl << "   Grau de Liberdade Condensado = " << nDofCondensed << std::endl;

					// Resolvendo o sistema linear
					an = new TPZAnalysis(mphysics, false);
					int numthreads = 0;
					TPZSkylineNSymStructMatrix strmat(mphysics);
					if (numthreads) strmat.SetNumThreads(numthreads);
					step = new TPZStepSolver<STATE>;

					step->SetDirect(ELU);
					an->SetSolver(*step);
					// end of time dedicated to geometric mesh processed
					t2 = clock();
					saida << "Process Creating CMesh : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
					t1 = clock();
					an->Run();
					t2 = clock();
					saida << "Process Run(mixed) : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
					t1 = clock();
					saida << "Process GMesh: \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
					saida << "Dimension of the linear system equations\n";
					saida << "K: " << an->Solver().Matrix()->Rows() << " x " << an->Solver().Matrix()->Cols() << "\n";
					//				saida << "K: " << step.Matrix()->Rows() << " x " << an->Solver()->Matrix()->Cols() << "F: " << an->Rhs()->Rows() << std::endl;

					TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);

					TPZVec<STATE> errorsHDiv;

					saida << "Erro da simulacao multifisica  para o Fluxo\n";
					ErrorHDiv2(meshvec[0], saida, errorsHDiv, fExact);
					t2 = clock();
					saida << "Process Computing Hdiv Errors : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
					t1 = clock();
					std::cout << "   Error Hdiv computed" << std::endl;

					an->SetExact(fExact);
					saida << "Erro da simulacao multifisica  para a Pressao\n";
					an->PostProcess(erro, saida);
					std::cout << "   Error L2 computed" << std::endl;
					t2 = clock();
					saida << "Process Computing L2 errors : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
					t1 = clock();
				}
				else if (numethod < 2) {
					// Determining the numeric method to use
					if (numethod) DGFEM = true;

					// To errors
					saida << "\n\nSaida do erro para formulacão ";
					if(DGFEM) saida << "GD,";
					else saida << "H1,";
					saida << " com ordem p  = " << simulationdata.pOrder << "\n";
					if (!fromgid.empty())
						saida << "Malha irregular from gid file.\n";

					/// PRIORITARY -> RELATIVE TO COMPUTATIONAL MESH
					///Indicacao de malha DGFem. Se false, vamos criar malha H1
					cmesh = new TPZCompMesh(gmesh);
					CreatingComputationalMesh(cmesh, ElementIDMat, bcdirichlet, DGFEM, DimProblem, simulationdata.pOrder, SourceFunc, solExata);
					std::cout << "   Computational mesh created:\n" << "   Grau de Liberdade Total = " << cmesh->NEquations() << std::endl;

					// PRIORITARY -> RELATIVE TO SOLVING THE PROBLEM
					an = new TPZAnalysis(cmesh);

					///Criando structural matrix - skyline não simétrica por causa do DGFEM usado
					TPZSkylineNSymStructMatrix matriz(cmesh);
					if (!stepsolver) {
						///Decomposição LU ou LDLt ---> Resolvendo via método direto
						step = new TPZStepSolver<STATE>;
						step->SetDirect(ELU);
					}
					else if (stepsolver == 1) {
						// Resolvendo via método iterativo GMRES
						TPZAutoPointer<TPZMatrix<STATE> > matbeingcopied = matriz.Create();
						TPZAutoPointer<TPZMatrix<STATE> > matClone = matbeingcopied->Clone();

						TPZStepSolver<STATE>* precond = new TPZStepSolver<STATE>(matClone);
						step = new TPZStepSolver<STATE>(matbeingcopied);
						precond->SetReferenceMatrix(matbeingcopied);
						precond->SetDirect(ELU);
						step->SetGMRES(40, 20, *precond, 1.e-16, 0);
					}
					else if (stepsolver == 2) {
						// Resolvendo via método iterativo - Gradiente conjugado
						TPZAutoPointer<TPZMatrix<STATE> > matbeingcopied = matriz.Create();
						TPZAutoPointer<TPZMatrix<STATE> > matClone = matbeingcopied->Clone();

						TPZStepSolver<STATE>* precond = new TPZStepSolver<STATE>(matClone);
						step = new TPZStepSolver<STATE>(matbeingcopied);
						precond->SetReferenceMatrix(matbeingcopied);
						precond->SetDirect(ELU);
						step->SetCG(20, *precond, 1.0e-12, 0);
					}
					an->SetSolver(*step);
					t2 = clock();
					saida << "Process Creating CMesh : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
					t1 = clock();

					///Assemble da matriz de rigidez e vetor de carga
					an->Assemble();
					t2 = clock();
					saida << "Process Assembling : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
					t1 = clock();
					saida << "Dimension of the linear system equations\n" << "K: " << an->Solver().Matrix()->Rows() << " x " << an->Solver().Matrix()->Cols() << "\n";
					///Resolução do sistema
					an->Solve();
					t2 = clock();
					saida << "Process Solving : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
					t1 = clock();
					std::cout << "   Approximated solution computed." << std::endl;

					// UNNECESSARY -> RELATIVE TO CALCULATING ERROR AFTER SOLVE PROCESS
					///Calculando erro da aproximacao
					an->SetExact(fExact);
					saida << "\nMODEL " << modelo << "   Element: " << typeel << " - Process " << simulationdata.pOrder;
					an->PostProcess(erro, saida);///calculando erro
					saida << " Erro de aproximacao:\n" << "Norma H1 = " << erro[0] << "\nNorma L2 = " << erro[1] << "\nSeminorma H1 = " << erro[2] << "\n\n";
					std::cout << "   Approximated error computed." << std::endl;
					t2 = clock();
					saida << "Process Computing L2 error : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
				}
				if (step) delete step; step = 0;

				// PRIORITARY -> RELATIVE TO PRINT SOLUTION TO VISUALIZATION BY PARAVIEW
				///Exportando para Paraview
				if (an) {
					std::stringstream ssout;
					ssout << sout;
					ssout << "Model" << modelo << "_" << DimProblem << "D_MESH_E" << typeel << "H" << simulationdata.nref << "_p" << simulationdata.pOrder;
					if (numethod == 2) ssout << "_Hdiv";
					if (DGFEM) ssout << "_GD";
					else ssout << "_H1";
					ssout << ".vtk";
					
					an->DefineGraphMesh(DimProblem, simulationdata.sV, simulationdata.vV, ssout.str());
					int resolution = 1;
					an->PostProcess(resolution);
					
					delete an; an = 0;
					std::cout << "   File to visualization of the approximated variables was saved." << std::endl << std::endl;
				}
				if (mphysics) delete mphysics; mphysics = 0;
				if (cmesh) delete cmesh; cmesh = 0;
				clock_t tmethod2 = clock();
				saida << "Process Numerical Method " << numethod << " : \n\tt_i = " << tmethod << "\tt_f = " << tmethod2 << "\tTime elapsed = " << tmethod2 - tmethod << std::endl;
			}
			if (gmesh) delete gmesh; gmesh = 0;
			saida.close();
			if (simulationdata.typeel == 20 && nelementtypes == 3)
				nelementtypes = 8;
			// end of time dedicated to process this typeel
			clock_t ttypeel2 = clock();
			saida << "Process Typeel " << typeel << " : \n\tt_i = " << ttypeel << "\tt_f = " << ttypeel2 << "\tTime elapsed = " << ttypeel2 - ttypeel << std::endl;
		} while (nelementtypes < 7);
		clock_t tmodelo2 = clock();
		saida << "Process Model " << modelo << " : \n\tt_i = " << tmodelo << "\tt_f = " << tmodelo2 << "\tTime elapsed = " << tmodelo2 - tmodelo << std::endl;
	}
	clock_t texec2 = clock();
	saida << "Process Simulator : \n\tt_i = " << texec << "\tt_f = " << texec2 << "\tTime elapsed = " << texec2 - texec << std::endl;
	saida.close();
	return 0;
}

void PutExact(int n) {
	if(n==1) fExact = &SolExactSeno;
	else if(n==2) fExact = &SolExactSenoSeno;
	else if(n==3) fExact = &SolExactSenoSenoSeno;
	else if(n==4) fExact = &SolExactArcTg1D;
	else if(n==5) fExact = &SolExactArcTg2D;
	else if(n==6) fExact = &SolExactArcTg3D;
	else if(n==7) fExact = &SolExactStrongOsc1D;
	else if(n==8) fExact = &SolExactStrongOsc2D;
	else if(n==9) fExact = &SolExactStrongOsc3D;
	else fExact = 0;
}

/* Computing gradient reconstructed
TPZFMatrix<REAL> gradients;
PosProcessGradientReconstruction(cmesh, gradients);

// Printing to VTK
std::stringstream saida;
saida << sout;
saida << "Grad_UnifMeshRUnif_p" << pOrder << ".vtk";
PrintDataMeshVTK(cmesh, saida.str(), gradients);
*/

/* Name of the Variables to output
TPZVec<std::string> scalarVars(4), vectorVars(2);
scalarVars[0] = "Solution";       // Uh - aproximated solution calculated
//		scalarVars[1] = "ExactPressure";  // u(x,y) - exact solution given
		scalarVars[1] = "ExactSolution";  // u(x,y) - exact solution given
		scalarVars[2] = "NormKDu";
		scalarVars[3] = "KDuDx";

		vectorVars[0] = "Flux";           // sigma_h - approximated flux computed from Uh
		vectorVars[1] = "ExactFlux";      // sigma - exact flux given
*/
