#include <stdlib.h>

#include "pzcmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include "pzvec.h"
#include "TPZMaterial.h"

#include "EEGradientReconstruction.h"

TEEGradientReconstructionLS::TEEGradientReconstructionLS(TPZCompMesh *cmesh, REAL majorc, REAL mediumc, REAL minorc, int var) : TErrorEstimator(cmesh,majorc,mediumc,minorc,var) {
}
TEEGradientReconstructionLS::~TEEGradientReconstructionLS() {
}

/* Function as error estimator for each element */
REAL TEEGradientReconstructionLS::ErrorEstimator(TPZCompEl *cel,int var) {
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

/* Compute gradient with gradiente reconstruction by least square using neighboards */
REAL TEEGradientReconstructionLS::ErrorEstimator(TPZCompEl *cel,TPZFMatrix<REAL> &grad,int var) {
	// Nada sera realizado para elementos con dimension diferente de la dimension del problema
	if(!cel) return 100.;
	int dim = cel->Mesh()->Dimension();
	if(cel->Dimension()!=dim) return 100.;
	TPZStack<TPZCompElSide> neighs;
	int nneighs;
	int nstates = cel->Material()->NSolutionVariables(var);
	int k, side, nsides = cel->Reference()->NSides()-1;   // Desconsiderando o proprio elemento (lado)
	
	TPZManVector<REAL,3> center(3,0.0), centerbeta(3,0.0);
	TPZManVector<STATE> solalfa(nstates,0.0), solbeta(nstates,0.0);
	TPZFMatrix<REAL> A(dim,dim);    // Linear System matrix
	TPZFMatrix<REAL> B(dim,1,0.);   // Linear System vector
	
	// Creando las matrices para aplicar el metodo de los minimos cuadrados
	TPZFMatrix<REAL> DeltaH;
	TPZFMatrix<REAL> DeltaHTranspose;
	TPZFMatrix<REAL> DifSol;
	
	// Limpiando las matrizes
	A.Zero();
	// Encontramos el centro del elemento corriente cel
	TPZGeoEl* gelalfa = cel->Reference();
	TPZManVector<REAL> centerpsi(gelalfa->Dimension(),0.0);
	gelalfa->CenterPoint(gelalfa->NSides()-1,centerpsi);
	center.Fill(0.);
	gelalfa->X(centerpsi,center);
	cel->Solution(centerpsi,var,solalfa);
	neighs.Resize(0);
	// Procuramos todos los elementos vecinos a cel (sobre todos los lados) sin duplicados
	for(side = 0; side < nsides; side++) {
		TPZCompElSide celside(cel,side);
		celside.ConnectedElementList(neighs,0,0);
	}
	nneighs = neighs.NElements();
	if(!nneighs) return 0.;
	
	// If exist neighboard clean duplicated elements and store only index of neighboard
	TPZStack<int> realneighs;
	int kk,jj;
	int id;   // = neighs[0].Element()->Index();
	for(kk=0;kk<nneighs;kk++) {
		id=neighs[kk].Element()->Index();
		for(jj=0;jj<realneighs.NElements();jj++) {
			if(id == realneighs[jj])
				break;
		}
		if(jj==realneighs.NElements() && cel->Mesh()->ElementVec()[id]->Dimension() == dim)
			realneighs.Push(id);
	}
	nneighs = realneighs.NElements();
	// si no hay vecinos continuamos con el siguiente elemento
	if(!nneighs) return 0.;
	
	// si hay vecinos realizamos el proceso de minimos quadrados para calcular una aproximacion del gradiente
	// Para cada vecino calculamos los deltaH (desde su centro al centro del elemento corriente)
	// y el valor de la solucion en su centro solbeta
	DeltaH.Redim(nneighs,dim);
	DeltaHTranspose.Redim(dim,nneighs);
	DifSol.Redim(nneighs,1);
	// Montando la matriz de los deltas DeltaH y de las diferencias de las soluciones DifSol
	for(int ineighs=0;ineighs<nneighs;ineighs++) {
		TPZGeoEl * gelbeta = cel->Mesh()->ElementVec()[realneighs[ineighs]]->Reference();
		if(!gelbeta)
			DebugStop();
		centerpsi.Fill(0.0);
		centerbeta.Fill(0.0);
		gelbeta->CenterPoint(gelbeta->NSides()-1,centerpsi);
		gelbeta->X(centerpsi,centerbeta);
		gelbeta->Reference()->Solution(centerpsi,var,solbeta);
		
		for(k=0;k<dim;k++)
			DeltaH(ineighs,k) = centerbeta[k] - center[k];
		DifSol(ineighs,0) = solbeta[var] - solalfa[var];
	}
	// Resolviendo el sistema por los minimos cuadrados: DeltaH_t * DifSol = DeltaH_t * DeltaH * Grad(u)
	DeltaH.Transpose(&DeltaHTranspose);
	grad = DeltaHTranspose*DifSol;
	A = DeltaHTranspose*DeltaH;
	A.SolveDirect(grad,ELU);

	REAL GradNorm = 0.;
	for(int ij=0;ij<grad.Rows();ij++)
		GradNorm += grad(ij,0)*grad(ij,0);
	return sqrt(GradNorm);
}
