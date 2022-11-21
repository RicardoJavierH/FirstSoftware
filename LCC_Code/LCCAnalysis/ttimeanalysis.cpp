/*******       File : tarungek.c

This file contains the method definitions for class TTimeAnalysis.

*******              *******/

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "ttimeanalysis.h"
#include "pzsolve.h"
#include "pzstepsolver.h"
#include "pzbndmat.h"
#include "pzmvmesh.h"
#include "pzdxmesh.h"
#include "pzv3dmesh.h"
#include "pzgnode.h"
#include "pzcompel.h"
#include "pzvec.h"
#include "pzgeoel.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "pzquad.h"
#include "pzbndcond.h"

#include "pzstrmatrix.h"

#include "pzelmat.h"
#include "myheader.h"
#include "commonJC.h"
#include "hadaptive.h"

#include "EEGradientReconstruction.h"

/*******       TTimeAnalysis, class derived from TPZAnalysis       *******/
TTimeAnalysis::TTimeAnalysis(std::istream &input,TPZCompMesh *mesh,int level)
     : TPZAnalysis(mesh,false), fOrder(0), fEndTime(0.), fTimeStep(0.0001) {
    fUk = 0;
    fName = "Time_Analysis";

    GetDataCommented(&input,fStartTime);
    GetDataCommented(&input,fEndTime);
    GetDataCommented(&input,MINIMETIMESTEP);
         
    TPZMaterial *mat = 0;
    std::map<int ,TPZMaterial * > materialmap(mesh->MaterialVec());
    std::map<int ,TPZMaterial * >::iterator it;
    for (it = materialmap.begin(); it != materialmap.end() ; it++)
    {
        mat = it->second;
        if(mat->Id()<0) continue;
        fLaw = (TConservationLaw *)mat;
        break;
    }
    fLaw->SetTime(fTimeStep,fStartTime);

  /**Adaptive or non-adaptive scheme*/
    fAdaptive = 0;
    
    fCFL = 1.;
    fCFLDiffussion = 0.;
    AppendMethodName(Name());
    fUZero = 0;
    fLevel = level+1;
         
#ifdef PARALLELVERSION
    fRank = 0;
    fSize = 1;
#endif
}

TTimeAnalysis::TTimeAnalysis(std::ostream &out) : TPZAnalysis(), fEndTime(0.) {
  fLaw = 0;
  fOrder = 1;
	fRhsAdaptNorm = 0.75;
  fStartTime = 0.;
	fNEndTimes = 1;
	fEndTime = 1.;
  MAXROWS = 700;
  fUk = 0;
  fUZero = 0;
	fAdaptive = 0;
	fCFLDiffussion = 0.;
  MINIMETIMESTEP =  1.e-7;
/*  PartialStiff = 0;
  PartialRhs = 0;
  PartialSol = 0;*/
#ifdef PARALLELVERSION
  fRank = 0;
  fSize = 1;
#endif
}

TTimeAnalysis::~TTimeAnalysis() {
  if(fUk) delete fUk;
//  if(PartialStiff) delete PartialStiff;
//  if(PartialRhs) delete PartialRhs;
//  if(PartialSol) delete PartialSol;
	if(fAdaptive) delete fAdaptive;
}

// Para recuperar el orden del metodo en el tiempo
void TTimeAnalysis::GetOrder(std::ifstream &input) {
    GetCommentary(input, 2);
    fOrder = 1;
}

/**Imprime estado atual da malha geometrica, computacional, solucao, etc*/
void TTimeAnalysis::PrintCurrentState() {
  /**Imprime objeto analysis corrente*/
  std::ofstream anout("analysis.dat");
  anout << "DATA TIME_ANALYSIS" << std::endl;
  anout << "Name analysis" << std::endl;
  anout << fName<< std::endl;
  anout << "Order Runge-Kutta" << std::endl;
  anout << fOrder << std::endl;
  anout << "Initial Time" << std::endl;
  anout << fStartTime << std::endl;
  anout << "End Time" << std::endl;
  anout << fEndTime << std::endl;
  anout << "CFL_Condition" << std::endl;
  anout << fCFL << std::endl;
  anout << "Coeficiente de difussividade" << std::endl;
  anout << fCFLDiffussion << std::endl;
  anout << "Refinement Level" << std::endl;
  anout << fLevel << std::endl;
//  anout << "Use Limiter" << endl;
//  anout << fUseLimiter << endl;
  anout << "Current Time (law)" << std::endl;
  anout << CurrentTime();
  anout.close();
//  anout.open("gmesh.dat");
//  fCompMesh->Reference()->PrintData(anout);
//  anout.close();
//  anout.open("cmesh.dat");
//  fCompMesh->PrintData(anout);
//  anout.close();
}

int TTimeAnalysis::OrderLaw() {
    return fLaw->NStateVariables();
}

double TTimeAnalysis::CurrentTime() {
  if(fLaw) return fLaw->Time();
  return 0.;
}

/**To adequate the stiffness matrix and rhs depending on the fSolution vector.
   It is necessary when the mesh was refined and fSolution was interpolated*/
void TTimeAnalysis::AdequateMatrix() {
  int numeq = fCompMesh->NEquations();
  if(!fCompMesh || !fStructMatrix || !fSolver) {
    PZError << "TTimeAnalysis::AdequateMatrix. Hasn't stiffness or rhs matrix, or mesh.\n";
    return;
  }
	if(fRhs.Rows() != numeq) {
    fRhs.Redim(numeq,1);
    fSolution.Redim(numeq,1);
    fUk->Redim(numeq,1);
    fCompMesh->InitializeBlock();

    if(fSolution.Rows() != numeq)
      PZError << "TTimeAnalysis::AdequateMatrix has incompatible solution.\n";
  }
}


int TTimeAnalysis::MatrixTruncated(int n) {
    int i, band = 10;   /// Jorge ??? ((TPZAutoPointer<TPZFBMatrix<STATE> >)(fSolver->Matrix()))->GetBand();
  int row = n*MAXROWS-1;
  int nvar = OrderLaw(), numeq = fCompMesh->NEquations();
  if((row+nvar*band+1) >= numeq) return numeq;
  int j, k;
  char mask;
	TPZAutoPointer<TPZMatrix<STATE> > matrix = fSolver->Matrix();
  for(i=0;i<(nvar+1)*band;i++) {
    mask = 'n';
    /**Achando linha cujo ultimo valor nao zero esteja na diagonal */
    for(j=0;j<band;j++) {
      if(matrix->GetVal(row,row+j+1)!=0.) { mask = 'y'; break; }
    }
    /**Verificando se estamos na linha final do bloco achado para truncamento */
    k = 1;
    while(mask=='n' && k<band) {
      for(j=0;j<band-k;j++) {
        if(matrix->GetVal(row-k,row+1+j)!=0.) { mask = 'y'; break; }
      }
      k++;
    }
    row++;
    /**Verificando que a proxima linha pertence a um outro bloco */
    if(mask=='n') return row;
  }
  std::cout << "TTimeAnalysis::MatrixTruncated. Insecure partition matrix." << std::endl;
  return numeq;
}

void TTimeAnalysis::GetSchemeType(char *filename) {
  /**We can to choose implicit or explicit scheme to diffusion*/
  char aux[16];
  int index=0;
  aux[index++] = 'R'; aux[index++] = 'K';
  aux[index++] = '1';
  std::cout << "Time Analysis Adaptive - Parallel and Implicit scheme.\n";
  aux[index++] = 'p'; aux[index++] = 'i';
	aux[index++] = '\0';
  strncat(filename,aux,index);
}

void TTimeAnalysis::Run(std::ostream &out) {
    if(!fCompMesh || !fLaw) {
        PZError << "TTimeAnalysis::Run, hasn't computational mesh or equation associated.\n";
        return;
    }
    
    double ctime = 0.;
#if COMPUTETIME
    /** Put initial time in fTimes and into conservation law*/
    ctime = fStartTime;
    clock_t t2, t1 = clock();
#endif
    
    fLaw->SetParameters(fExplicit,fCFL,fCFLDiffussion);

    int counter2print = 0;
    while (ctime < fEndTime) {
        fTimeStep = fLaw->CFLCondition(fCompMesh);
        fLaw->SetTime(fTimeStep,ctime);

        ///Assemble da matriz de rigidez e vetor de carga
        Assemble();

#if COMPUTETIME
        t2 = clock();
        out << "Process Assembling : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
        t1 = clock();
        out << "Dimension of the linear system equations\n" << "K: " << Solver().Matrix()->Rows() << " x " << Solver().Matrix()->Cols() << "\n";
#endif
        ///Resolu‹o do sistema
        Solve();
		 // Applying error estimator
		 fErrorEstimator->ApplyingEstimation();
		 fAdaptive->Adapting(fErrorEstimator);

#if COMPUTETIME
        t2 = clock();
        out << "Process Solving : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
#endif
       // Actions by iteration
       ctime += fTimeStep;
       counter2print++;
       if(!(counter2print%5)) {
       //  DefineGraphMesh(fCompMesh->Dimension(), fScalarNames,fVectorNames, plotfile);
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

int TTimeAnalysis::Adapting(int &stepgraph,int onlyflux) {
	/** Clean object to adaptive work if is the case */
	if(!fAdaptive) return 0;

  /**Here: We can to change the mesh and acerting stiff matrix, solver and fUk*/
  AdequateMatrix();
	return 1;
}

void TTimeAnalysis::ComputeNextTime(double ctime,REAL EndTime) {
  /*Computing time step from CFL Condition*/
  fTimeStep = StabilityCondition();
  /*Adjusting dt whether it is very small*/
  if(fTimeStep<MINIMETIMESTEP) fTimeStep = MINIMETIMESTEP;

  /**Acerting time step near of end time and store current time*/
  if(ctime+fTimeStep > EndTime)
    fTimeStep = EndTime - ctime;
}

void TTimeAnalysis::Print(char *name,std::ostream &out,int all) {
  out << std::endl << name << std::endl << std::endl << "TTimeAnalysis Analysis" << std::endl;
//  out << "\tRungeKutta method order    = " << fOrder << std::endl;
  out << "\tSpatial Analysis used name = " << Name() << std::endl;
  out << "\tInitial Time               = " << fStartTime << std::endl;
  out << "\tFinish Time                = " << fEndTime << std::endl;
  if(all) PrintSolution(out);
}

void TTimeAnalysis::PrintSolution(std::ostream &out) {
  if(!fCompMesh) {
    std::cout << "TTimeAnalysis::PrintSolution hasn't mesh associated.\n";
    return;
  }
  if(LastTimeReached()==fEndTime) {
    out << std::endl << "Time reached = " << fEndTime << std::endl << "Connect Solution :" << std::endl;
    fCompMesh->ConnectSolution(out);
    return;
  }
  out << std::endl << "End time is not reached" << std::endl;
  out << "current time = " << LastTimeReached() << std::endl;
  fCompMesh->ConnectSolution(out);
}

double TTimeAnalysis::StabilityCondition() {
    REAL area;
    ComputeAreas(Mesh(),&area);
//  if(fLaw->Dimension()==2) area = sqrt(area);
  REAL maxeig = fLaw->MaxEigen();
	if(maxeig < 0.) PZError << "TTimeAnalysis::StabilityCondition. Bad maxeig.\n";
  if(IsZero(maxeig*.1)) return MINIMETIMESTEP;
  if(maxeig > BIGMAXEIG) return fEndTime;
  return (fCFL*area)/maxeig;
}

void TTimeAnalysis::ApplyUZero() {
    if(!fUZero) {
        PZError << "TTimeAnalysis::ApplyUZero. Function UZero() is undefined.\n";
        return;
    }
    TPZBlock<STATE> &block = fCompMesh->Block();
    block.SetMatrix(&fSolution);
    block.Matrix()->Zero();

    int64_t i, index, nelem = fCompMesh->NElements();
    TPZInterpolatedElement *cel;
    int jn, kn, dfseq, dfvar = fLaw->NStateVariables();
    TPZVec<STATE> u(dfvar,0.);
    TPZVec<int64_t> MaskConnects(fCompMesh->NConnects(),0);
    int eltype;

    TPZManVector<REAL> point(3,0.);
    for(i=0;i<nelem;i++) {
        cel = (TPZInterpolatedElement *)fCompMesh->ElementVec()[i];
        if(!cel || cel->IsInterface() || cel->Dimension()<fCompMesh->Dimension()) continue;
        int order = cel->GetPreferredOrder();
        int iphi0 = 3, iphi3 = 4, iphi1 = 1, iphi2 = 2;
        eltype = cel->Type();
        if(eltype < 9) {
            // Preenchendo o valor de UZero para os connects esquinas.
            for(jn=0;jn<cel->NCornerConnects();jn++) {
              index = cel->ConnectIndex(jn);
              if(index<0 || MaskConnects[index]) continue;
              cel->Reference()->Node(jn).GetCoordinates(point);
              
              fUZero(point,u,cel);
              dfseq = cel->Connect(jn).SequenceNumber();
              if(dfvar!=block.Size(dfseq))
                PZError << "TTimeAnalysis::ApplyUZero. Bad size." << std::endl;
              for(kn=0;kn<dfvar;kn++)
                block(dfseq,0,kn,0) = u[kn];
              MaskConnects[index] = 1;
            }
        }
        else {  // eltype == 15 (Discontinuous)
            // Preenchendo o valor de UZero para os connects esquinas.
            index = cel->ConnectIndex(0);
            if(index<0 || MaskConnects[index]) continue;
            REAL factor = 0.25;
            int ncorners = cel->Reference()->NCornerNodes();
            TPZConnect conector = cel->Connect(0);
            dfseq = conector.SequenceNumber();
            if(!order) {
                TPZVec<REAL> qsi(cel->Dimension(),0.);
                cel->Reference()->CenterPoint(cel->Reference()->NSides()-1,qsi);
                cel->Reference()->X(qsi,point);
                fUZero(point,u,cel);
                for(kn=0;kn<dfvar;kn++)
                    block(dfseq,0,kn,0) += u[kn];
            }
            else {
                for(jn=0;jn<ncorners;jn++) {
                    cel->Reference()->Node(jn).GetCoordinates(point);
                    fUZero(point,u,cel);
                    for(kn=0;kn<dfvar;kn++) {
                        if(ncorners==4) {
                            if(order==1) {
                                iphi0 = 0; iphi3 = 3;
                            }
                            switch(jn) {
                            case 0:
                                block(dfseq,0,kn+dfvar*iphi0,0) += factor*u[kn];
                                block(dfseq,0,kn+dfvar*iphi1,0) += -1*factor*u[kn];
                                block(dfseq,0,kn+dfvar*iphi2,0) += -1*factor*u[kn];
                                block(dfseq,0,kn+dfvar*iphi3,0) += factor*u[kn];
                            break;
                            case 1:
                                block(dfseq,0,kn+dfvar*iphi0,0) += -1.*factor*u[kn];
                                block(dfseq,0,kn+dfvar*iphi1,0) += -1.*factor*u[kn];
                                block(dfseq,0,kn+dfvar*iphi2,0) += factor*u[kn];
                                block(dfseq,0,kn+dfvar*iphi3,0) += factor*u[kn];
                            break;
                            case 2:
                                block(dfseq,0,kn+dfvar*iphi0,0) += factor*u[kn];
                                block(dfseq,0,kn+dfvar*iphi1,0) += factor*u[kn];
                                block(dfseq,0,kn+dfvar*iphi2,0) += factor*u[kn];
                                block(dfseq,0,kn+dfvar*iphi3,0) += factor*u[kn];
                            break;
                            case 3:
                                block(dfseq,0,kn+dfvar*iphi0,0) += -1.*factor*u[kn];
                                block(dfseq,0,kn+dfvar*iphi1,0) += factor*u[kn];
                                block(dfseq,0,kn+dfvar*iphi2,0) += -1*factor*u[kn];
                                block(dfseq,0,kn+dfvar*iphi3,0) += factor*u[kn];
                            break;
                            default:
                                break;
                            }
                        }
                        else if(ncorners==3) {
                            if(order==1) iphi0 = 0;
                            switch(jn) {
                                case 0:
                                    block(dfseq,0,kn+dfvar*iphi0,0) += -3.*u[kn];
                                    block(dfseq,0,kn+dfvar*iphi1,0) += -2*u[kn];
                                    block(dfseq,0,kn+dfvar*iphi2,0) += -2*u[kn];
                                    break;
                                case 1:
                                    block(dfseq,0,kn+dfvar*iphi0,0) += -3*u[kn];
                                    block(dfseq,0,kn+dfvar*iphi1,0) += -1*u[kn];
                                break;
                                case 2:
                                    block(dfseq,0,kn+dfvar*iphi0,0) += -3.*u[kn];
                                    block(dfseq,0,kn+dfvar*iphi2,0) += -1.*u[kn];
                                break;
                                default:
                                    break;
                            }
                        }
                    }
                }
            }
            MaskConnects[index] = 1;
        }
    }
    LoadSolution();
}

void TTimeAnalysis::Integral(TPZInterpolatedElement *el,
     void (*fp)(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZInterpolatedElement *cel),TPZVec<STATE> &integral) {
  int i,dim = el->Dimension();
  TPZFMatrix<REAL> axes(3,3,0.);
  TPZFMatrix<REAL> jacobian(dim,dim);
  TPZFMatrix<REAL> jacinv(dim,dim);
  REAL detjac;
  TPZVec<REAL> x(3,0.);
  TPZVec<REAL> intpoint(dim,0.);
  REAL weight = 0.;

  TPZIntPoints *rule;
  rule = &(el->GetIntegrationRule());

  int order = el->Material()->NStateVariables();
  TPZVec<STATE> u(order,0.);
  for(i=0;i<order;i++) integral[i] = 0.;

  for(i=0;i<rule->NPoints();++i) {
    rule->Point(i,intpoint,weight);
    el->Reference()->Jacobian(intpoint,jacobian,axes,detjac,jacinv);
    el->Reference()->X(intpoint, x);
    weight *= detjac;

    fp(x,u,el);
    for(i=0;i<order;i++) integral[i] += (u[i]*weight);
  }
}

/**Return maxime area and minime area using the pointer paramenter*/
REAL TTimeAnalysis::ComputeAreas(TPZCompMesh *cmesh,REAL *areamin) {
  int i,nelems = cmesh->NElements();
  REAL area, areamax = 0.;
  if(areamin) (*areamin) = 1.e12;
  TPZCompEl *el;
  int dim = fLaw->Dimension();

  for(i=0;i<nelems;i++) {
    el = fCompMesh->ElementVec()[i];
    if(!el || el->IsInterface()) continue;
    TPZInterpolatedElement *intel = (TPZInterpolatedElement *)el;
    if(intel->Dimension() != dim) continue;
    area = intel->Reference()->Volume();
    areamax = (area < areamax) ? areamax : area;
    (*areamin) = (area < (*areamin)) ? area : (*areamin);
  }
  return areamax;
}

void TTimeAnalysis::CleanToStartRun() {
  fLaw->Clean();
  fLaw->SetTime(fStartTime,fTimeStep);
  fTimeStep = fCFL = 0.;
  if(!fUk) fUk = new TPZFMatrix<STATE>(fSolution);
  else (*fUk) = fSolution;
}

void TTimeAnalysis::ReadData(std::ifstream &input) {
	GetCommentary(input,1);
	int limiter;
	REAL diff;
  GetDataCommented(input,limiter);
	GetDataCommented(input,fSteadyState);
  GetDataCommented(input,fCFL);
  GetDataCommented(input,diff);
  GetDataCommented(input,fRhsAdaptNorm);
  GetDataCommented(input,fNPlots);
}

/*
void TTimeAnalysis::DefineGraphMesh(int dim,TPZVec<char *> &scalnames,
			       TPZVec<char *> &vecnames,TPZGraphMesh *graph) {
  int dim1 = dim-1;
  TPZDrawStyle style = graph->Style();

  if(fGraphMesh[dim1]) delete fGraphMesh[dim1];
  fScalarNames[dim1] = scalnames;
  fVectorNames[dim1] = vecnames;

  switch(style) {
  case 0:
    fGraphMesh[dim1] = new TPZDXGraphMesh(fCompMesh,dim,(TPZDXGraphMesh *)graph,fLaw);
    break;
  case 1:
    fGraphMesh[dim1] = new TPZMVGraphMesh(fCompMesh,dim,(TPZMVGraphMesh *)graph,fLaw);
    break;
  case 2:
    fGraphMesh[dim1] = new TPZV3DGraphMesh(fCompMesh,dim,(TPZV3DGraphMesh *)graph,fLaw);
    break;
  default:
    cout << "grafgrid was not created\n";
    fGraphMesh[dim1] = 0;
    break;
  }
}
*/
/*void TTimeAnalysis::Run(istream &input,int dimension, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const std::string &plotfile,ostream &out) {
    double ctime = 0.;

    GetDataCommented(input,fExplicit);
    fLaw->SetExplicit(fExplicit);

    if(!fCompMesh || !fLaw) {
        PZError << "TTimeAnalysis::Run, hasn't computational mesh or equation associated.\n";
        return;
    }

    // Put initial time in fTimes and into conservation law
    ctime = fStartTime;
    clock_t t2, t1 = clock();
    
    int count2print = 0;   // contador para impressao de solucoes
    while (ctime < fEndTime) {
        fTimeStep = fLaw->CFLCondition(fCompMesh);
        fLaw->SetTime(fTimeStep,ctime);
        if (!(count2print % 3)) {
            DefineGraphMesh(dimension, scalnames, vecnames, plotfile);
            int resolution = 0;
            LoadSolution();
            PostProcess(resolution);
        }
        count2print++;

        ///Assemble da matriz de rigidez e vetor de carga
        Assemble();
        
        t2 = clock();
        out << "Process Assembling : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
        t1 = clock();
        out << "Dimension of the linear system equations\n" << "K: " << Solver().Matrix()->Rows() << " x " << Solver().Matrix()->Cols() << "\n";
        ///Resolu‹o do sistema
        Solve();
        t2 = clock();
        out << "Process Solving : \n\tt_i = " << t1 << "\tt_f = " << t2 << "\tTime elapsed = " << t2 - t1 << std::endl;
        t1 = clock();
        std::cout << "   Approximated solution computed." << std::endl;

        // Actions by iteration
        ctime += fTimeStep;
        std::cout << "Current Time " << fLaw->Time() << std::endl;
    }
}
*/
