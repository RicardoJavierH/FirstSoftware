
#ifndef COMMOM_JCN_hpp
#define COMMOM_JCN_hpp

#include <fstream>
#include <cmath>
#include <stdio.h>

#include "pzvec.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"

#define BIGNUMBER 1.e5

#define COMPUTETIME 0

struct SimulationData {
	int Equation;
	int nuzero;
	int typeel;
	int print;
	int nref;
	int NumericMethod;
	int pOrder;
	int stepsolver;
	int nsV;
	TPZVec<std::string> sV;
	int nvV;
	TPZVec<std::string> vV;
	bool applyadaptivity;
	SimulationData();
};

void PrintInfo(TPZVec<int64_t>& indexsubelements, TPZVec<REAL>& means, TPZVec<REAL>& wavcoefs, std::string& sout);

void PrintInfo(int n, TPZVec<int64_t>& indexsubelements, TPZVec<REAL>& means, std::string& sout);

// To load refinement pattern
void LoadUniformPatterns(int typeel);

// Creating necessary directories and return name of directory to results files
int PutNameDirectoryForResults(std::string& sout);

// To chance order into computational mesh with Hdiv space - mixed formulation
void ChangeInternalConnectOrder(TPZCompMesh *mesh, int typeel, int addToOrder);

// To import user data from filedata
std::ifstream* UserDataImport(std::string &fname, std::string& fromgid, SimulationData &simulationdata);
bool CheckTypeEl(int typeel);

// TO GRADIENT RECONSTRUCTION BY LEAST SQUARES
REAL MeanCell(TPZCompEl *cel, int IntOrder);
void GradientReconstructionByLeastSquares(TPZCompEl *cel, TPZManVector<REAL, 3> &center, TPZVec<REAL> &Grad);
// Generate an output of all geomesh to VTK, associating to each one the given data, creates a file with filename given
void PrintDataMeshVTK(TPZCompMesh *cmesh, std::string &filename, TPZFMatrix<REAL> &elData);
void PosProcessGradientReconstruction(TPZCompMesh *cmesh, TPZFMatrix<REAL> &datagradients);

void GetCommentary(std::istream &input,int nlines);
void GetDataCommented(std::istream &input,int &value);
void GetDataCommented(std::istream &input,int const nlines,int &value);
void GetDataCommented(std::istream &input,REAL &value);
void GetDataCommented(std::istream &input,TPZVec<REAL> &vector);
void GetDataCommented(std::istream &input,char *string,int size);
void GetDataCommented(std::istream &input, std::string &str);

void GetDataCommented(std::istream *input,REAL &value);
void GetDataCommented(std::istream *input,int &value);

void TesteNormal(TPZCompMesh *);

void Multiply(TPZMatrix<STATE> &A,TPZMatrix<STATE> &B,TPZMatrix<STATE> &result);
void InvJacob2d(TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &invjac);

int PrintDiagonal(TPZMatrix<REAL> *matrix, int mask = 0);

// Funcoes para o calculo do passo do tempo satisfazendo a condicao CFL -->>Deve ir para a classe TConservationLaw
REAL MinimumDeltaX(TPZCompMesh *mesh);

#endif
