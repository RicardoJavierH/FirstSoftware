/*******       File : tarungek.h

Header file for class TTimeAnalysis. This class has a weighted residual analysis
to solve partial differential equations with first order time derivated. An object
of this class apply the L-order Runge-Kutta method.

*******              *******/

#ifndef TIMEANALYSISHH
#define TIMEANALYSISHH

#include <stdlib.h>
#include <string>

#include "pzanalysis.h"
#include "conslaw.h"

#include "pzbndmat.h"
#include "pzfmatrix.h"

#include "ErrorEstimator.h"

#define BIGMAXEIG 3.e2
#define MAXITERATIONS 30000000

class TAdaptive;
template<class T> class TPZVec;

/**
 Implementa el mŽtodo numŽrico que considera para el tiempo el mŽtodo de Euler (impl’cito — expl’cito)
 y para el espacio puede ser GD — H1 */
class TTimeAnalysis : public TPZAnalysis {

protected:
    int          fOrder;        //Order of the Runge-Kutta method
    REAL         fStartTime;    //Initial time of the time domain
    int          fNEndTimes;    //Number of the partial end times to write post-processing files
    REAL        fEndTime;      //End time of the time domain
    std::string  fName;    //Associated name to time analysis
    REAL         fTimeStep;     //Time step intervale :  DeltaT
    REAL         fCFL;          //The spatial CFL to convergence
    REAL         fCFLDiffussion; //The time CFL to increment diffusion
	int          fExplicit;     // 0 (implicito) =1(explicito com difusividade) >1(explicito Runge-Kutta)
	int          fSteadyState;  // 0 se o processo é nao estacionario, 1 se estacionario

	/** Pointer to object does adaptivity */
	TAdaptive *fAdaptive;
	/* Pointer to make error estimation */
	TErrorEstimator *fErrorEstimator;

    REAL MINIMETIMESTEP;    //Minimo passo de tempo permitido
    int MAXROWS;              //Maximo numero de linhas ao particionar as matrices

    TPZFMatrix<STATE>   *fUk;          //Solution vector to store solution in current time
    TConservationLaw  *fLaw;    //Conservation law to solve. Improve to several equations

    /** Norm of the residual vector */
    REAL fRhsNorm;
    /** Threashold value for aplying adaptivity on the norm */
    REAL fRhsAdaptNorm;
    TPZFMatrix<STATE> fPreviousSol;

public:
	void SetErrorEstimator(TErrorEstimator *estimator) {
		fErrorEstimator = estimator;
	}
	void SetAdaptive(TAdaptive *adaptive) {
		fAdaptive = adaptive;
	}
	

    /** To indicate the level until the comp. mesh was refined */
    int fLevel;
    int fNPlots;

    TTimeAnalysis(std::istream &input,TPZCompMesh *mesh,int level=1);
    TTimeAnalysis(std::ostream &out);
    ~TTimeAnalysis();

    /** Pointer to initial function at T=0 */
    void (*fUZero)(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZInterpolatedElement *cel);

    void SetDeltaT(double dt) { fTimeStep = dt; }
    REAL DeltaT() { return fTimeStep; }
    void SetName(std::string &name) { fName = name; }
    void SetLaw(TConservationLaw *law) { fLaw = law; }
    TConservationLaw *GetLaw() { return fLaw; }
    void SetTimeDomain(double StartT,double EndT) {
        fStartTime = StartT;
        fEndTime = EndT;
    }
    double CurrentTime();

	void SetParameters(int formexplicit,REAL CFL,REAL DiffusionCoef);

    std::string &Name() { return fName; }
    double StartTime() { return fStartTime; }
    double EndTime() { return fEndTime; }
    double LastTimeReached() { return CurrentTime(); }
    int OrderLaw();

    /**Read complementary data to analysis */
    virtual void ReadData(std::ifstream &input);
    /**To clean TTimeAnalysis data members */
    virtual void CleanToStartRun();

    /**Compute time step from CFL condition and compute next time adjusting it at end time*/
    void ComputeNextTime(double currenttime,REAL EndTime);
    double StabilityCondition();
        
    // To solve problem
    /** @brief Calls the appropriate sequence of methods to build a solution or a time stepping sequence */
    virtual void Run(std::ostream &out = std::cout);
    int MatrixTruncated(int n);
    /**Return maxime area and return minime area using pointer parameter*/
    REAL ComputeAreas(TPZCompMesh *cmesh,REAL *areamin);

    /**To put real dimensions in stiffness, rhs and solution matrices*/
    void AdequateMatrix();

    /**Applying adaptive scheme*/
    virtual int Adapting(int &step,int onlyflux);

    void Print(char *namevar,std::ostream &out=std::cout,int all=0);
    void PrintSolution(std::ostream &out=std::cout);
    void AppendMethodName(std::string &name);

    /**Compute initial function into fSolution of the mesh*/
    void ApplyUZero();
    void SetUZero(void (*uzero)(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZInterpolatedElement *cel));
    void Integral(TPZInterpolatedElement *el,
            void (*fp)(TPZVec<REAL> &x,TPZVec<STATE> &u,TPZInterpolatedElement *cel),TPZVec<STATE> &integral);

    virtual void GetSchemeType(char *filename);
    /**Imprime estado atual da malha geometrica, computacional, solucao, etc*/
    void PrintCurrentState();

    virtual void GetOrder(std::ifstream &input);

#ifdef PARALLELVERSION
    int    fRank;   // process id
    int    fSize;   // number of process
    MPI_Status fStatus;
    int    fTypePar;
#endif
};

inline void TTimeAnalysis::AppendMethodName(std::string &name) {
    name += " Runge-Kutta.";
}

inline void TTimeAnalysis::SetUZero(void (*uzero)(TPZVec<REAL> &x,TPZVec<STATE> &u,
    TPZInterpolatedElement *cel)) {
    fUZero = uzero;
}

inline void TTimeAnalysis::SetParameters(int Explicit,REAL CFL,REAL DiffusionCoef) {
    fExplicit = Explicit;
    fCFL = CFL;
    fCFLDiffussion = DiffusionCoef;
}

#endif
