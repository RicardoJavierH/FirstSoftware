/*******       File : tarungek.c

This file contains the method definitions for class TTimeRungeKutta.

*******              *******/

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "pzcmesh.h"

#include "tarungek.h"
#include "pzsolve.h"
#include "pzstepsolver.h"
#include "pzbndmat.h"
#include "pzmvmesh.h"
#include "pzdxmesh.h"
#include "pzv3dmesh.h"
#include "conslaw.h"
#include "pzgnode.h"
#include "pzcompel.h"
#include "pzvec.h"
#include "pzgeoel.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "pzquad.h"

#include "pzstrmatrix.h"

#include "ErrorEstimatorAlmostHaarWavelets.h"

#include "pzelmat.h"
#include "myheader.h"
#include "commonJC.h"
#include "hadaptive.h"

/*******       TTimeRungeKutta, class derived from TPZAnalysis       *******/
TTimeRungeKutta::TTimeRungeKutta(std::istream &input,TPZCompMesh *mesh,int level)
     : TTimeAnalysis(input,mesh,level) {
   //      GetDataCommented(input,fOrder);
    int i;
         fOrder = 1;
    for(i=0;i<4;i++)
        fCoef[i][0] = fCoef[i][1] = 0.;
	fErrorEstimator = 0;
}

TTimeRungeKutta::TTimeRungeKutta(std::ostream &out) : TTimeAnalysis(out) {
  fLaw = 0;
	fErrorEstimator = 0;
  fOrder = 1;
    fRhsAdaptNorm = 0.75;
  int i;
  for(i=0;i<4;i++)
    fCoef[i][0] = fCoef[i][1] = 0.;
  fCoef[0][1] = 1.;
  fStartTime = 0.;
    fNEndTimes = 1;
    fEndTime = 1.;
  MAXROWS = 700;
  fUk = 0;
  fUZero = 0;
    fAdaptive = 0;
//    fCFLDiffussion = 0.;
  MINIMETIMESTEP =  1.e-7;
}

void TTimeRungeKutta::GetOrder(std::ifstream &input) {
    GetDataCommented(&input,fOrder);
    int i;
    for(i=0;i<4;i++)
        fCoef[i][0] = fCoef[i][1] = 0.;
    fCoef[0][1] = 1.;
    switch(fOrder) {
        case 1:
            break;
        case 2: {
            fCoef[1][0] = fCoef[1][1] = 1./2.;
            break;
        }
        case 3: {
            fCoef[1][0] = fCoef[1][1] = .25;
            fCoef[2][0] = 1./6.;
            fCoef[2][1] = 2./3.;
            break;
        }
        case 4: {
            fCoef[1][0] = 0.;
            fCoef[1][1] = 0.5;
            fCoef[2][0] = 0.;
            fCoef[2][1] = 1.;
            fCoef[3][0] = 1./3.;
            fCoef[3][1] = 1./6.;
            break;
        }
    }
}

void TTimeRungeKutta::Run(std::string &plotfile) {
    if(!fCompMesh || !fLaw) {
        PZError << "TTimeAnalysis::Run, hasn't computational mesh or equation associated.\n";
        return;
    }

    /** Put initial time in fTimes and into conservation law*/
    double ctime = 0.;
    ctime = fStartTime;
#if COMPUTETIME
    clock_t t2, t1 = clock();
#endif

    fLaw->SetParameters(fExplicit,fCFL,fCFLDiffussion);

    int counter2print = 0;
    TPZFMatrix<STATE> Solution0;
    while (ctime < fEndTime) {
		 // Applying error estimator
		 if(fErrorEstimator) {
			 //				  fErrorEstimator->ApplyingEstimation();
			 LoadSolution();
			 TPZVec<std::string> scalar = GraphMesh(2)->ScalarNames();
			 TPZVec<std::string> vec = GraphMesh(2)->VecNames();
			 int stepgrahp = GetStep();
			 // Resgatando as informações para os graficos no paraview
			 fAdaptive->Adapting(fErrorEstimator);
			 SetCompMesh(fCompMesh,0);     // nao utilizar optimized pois da erro no GD
			 DefineGraphMesh(fCompMesh->Dimension(),scalar,vec, plotfile);
			 SetStep(stepgrahp);
			 StructMatrix().operator->()->SetMesh(fCompMesh);
			 fCompMesh->Block().SetMatrix(&fSolution);
		 }
        fTimeStep = fLaw->CFLCondition(fCompMesh);
        if(fTimeStep<MINIMETIMESTEP)
            fTimeStep = MINIMETIMESTEP;
        if(fOrder) fTimeStep /= fOrder;
        Solution0 = Solution();
        for(int RKorder=0;RKorder<fOrder;RKorder++) {
            fLaw->SetTime(fTimeStep,ctime);
            REAL Coef = 0;
            Coef = fCoef[RKorder][1];
            ///Assemble da matriz de rigidez e vetor de carga
            fLaw->SetOrderAndCoeficient(RKorder+1, Coef);
            ///Assemble da matriz de rigidez e vetor de carga
            Assemble();
            if(RKorder) {
                TPZFMatrix<STATE> Rhs;
                Rhs = fRhs;
                fSolution = Solution0;
                LoadSolution();
                fLaw->SetOrderAndCoeficient(RKorder+1, fCoef[RKorder][0]);
                Assemble();
                fRhs += Rhs;
            }
#if COMPUTETIME
            t2 = clock();
            out << "Process Assembling : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
            t1 = clock();
            out << "Dimension of the linear system equations\n" << "K: " << Solver().Matrix()->Rows() << " x " << Solver().Matrix()->Cols() << "\n";
#endif
            ///Resolução do sistema
            Solve();
			  			  
#if COMPUTETIME
            t2 = clock();
            out << "Process Solving : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
#endif
            // Actions by iteration
            ctime += fTimeStep;
            counter2print++;
            if(RKorder == fOrder-1) {
                int resolution = 0;
                PostProcess(resolution);
            }
            std::cout << "=> Approximated solution computed. \nCurrent Time = " << ctime << std::endl;
#if COMPUTETIME
            t1 = clock();
#endif
			  fLaw->SetTime(fTimeStep,ctime);
        }
    }
}
void TTimeRungeKutta::Run(std::ostream &out) {
	if(!fCompMesh || !fLaw) {
		PZError << "TTimeAnalysis::Run, hasn't computational mesh or equation associated.\n";
		return;
	}
	
	/** Put initial time in fTimes and into conservation law*/
	double ctime = 0.;
	ctime = fStartTime;
#if COMPUTETIME
	clock_t t2, t1 = clock();
#endif
	
	fLaw->SetParameters(fExplicit,fCFL,fCFLDiffussion);
	
	int counter2print = 0;
	TPZFMatrix<STATE> Solution0;
	while (ctime < fEndTime) {
		fTimeStep = fLaw->CFLCondition(fCompMesh);
		if(fTimeStep<MINIMETIMESTEP)
			fTimeStep = MINIMETIMESTEP;
		if(fOrder) fTimeStep /= fOrder;
		Solution0 = Solution();
		for(int RKorder=0;RKorder<fOrder;RKorder++) {
			fLaw->SetTime(fTimeStep,ctime);
			REAL Coef = 0;
			Coef = fCoef[RKorder][1];
			///Assemble da matriz de rigidez e vetor de carga
			fLaw->SetOrderAndCoeficient(RKorder+1, Coef);
			///Assemble da matriz de rigidez e vetor de carga
			Assemble();
			if(RKorder) {
				TPZFMatrix<STATE> Rhs;
				Rhs = fRhs;
				fSolution = Solution0;
				LoadSolution();
				fLaw->SetOrderAndCoeficient(RKorder+1, fCoef[RKorder][0]);
				Assemble();
				fRhs += Rhs;
			}
#if COMPUTETIME
			t2 = clock();
			out << "Process Assembling : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
			t1 = clock();
			out << "Dimension of the linear system equations\n" << "K: " << Solver().Matrix()->Rows() << " x " << Solver().Matrix()->Cols() << "\n";
#endif
			///Resolução do sistema
			Solve();
			
#if COMPUTETIME
			t2 = clock();
			out << "Process Solving : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
#endif
			// Actions by iteration
			ctime += fTimeStep;
			counter2print++;
			if(RKorder == fOrder-1) {
				int resolution = 0;
				PostProcess(resolution);
			}
			std::cout << "=> Approximated solution computed. \nCurrent Time = " << ctime << std::endl;
#if COMPUTETIME
			t1 = clock();
#endif
			// Applying error estimator
			if(fErrorEstimator) {
				int ij;
				//				  fErrorEstimator->ApplyingEstimation();
				LoadSolution();
				fAdaptive->Adapting(fErrorEstimator);
				//this->SetCompMesh(cmesh,0);
				//			  fRhs.Redim(fCompMesh->Solution().Rows(),fCompMesh->Solution().Cols());
				//			  fSolution = fCompMesh->Solution();
				//				  fSolution.Redim(fCompMesh->Solution().Rows(),fCompMesh->Solution().Cols());
				//				  LoadSolution(fCompMesh->Solution());
				SetCompMesh(fCompMesh,1);
				StructMatrix().operator->()->SetMesh(fCompMesh);
				fCompMesh->Block().SetMatrix(&fSolution);
			}
			fLaw->SetTime(fTimeStep,ctime);
		}
	}
}


REAL TTimeRungeKutta::MinMod(REAL first,REAL second,REAL third) {
  int mask = (first<0) ? -1 : 1;
  if(((second<0)?-1:1)!=mask || ((third<0)?-1:1)!=mask) return 0.;
  REAL vafirst = (first<0) ? -first : first;
  REAL vasecond = (second<0) ? -second : second;
  REAL vathird = (third<0) ? -third : third;
  if(vafirst>vasecond)
    vafirst = vasecond;
  return mask*((vafirst<vathird) ? vafirst : vathird);
}

