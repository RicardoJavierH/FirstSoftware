#ifndef TIMERUNGEKUTTAHH
#define TIMERUNGEKUTTAHH

#include "ttimeanalysis.h"

/* Implementa el método numérico Runge Kutta de orde 1 (Euler), orden 2, 3 y 4 en la variable
 tiempo, siempre de forma explícita. En el espacio el usuário puede escojer GD ó H1 */
class TTimeRungeKutta : public TTimeAnalysis {
	
    REAL         fCoef[4][2];   //Coefficients to operator H(u). Son los gamma_ja y gamma_jb de la tesis de Jorge

 public:
  TTimeRungeKutta(std::istream &input,TPZCompMesh *mesh,int level=1);
  TTimeRungeKutta(std::ostream &out);
  ~TTimeRungeKutta() { }

	void Run(std::string &plotfile);
  void Run(std::ostream &out=std::cout);
//  void Solve();

  /**Operator Cockburn limiter - Generalized slope limiter*/
  /**It is used when the interpolated order is greater than zero : TPZCompEl::gOrder>0*/
//  void CockburnLimiter(TPZAdmChunkVector<TPZCompEl *> &elvec,int var);
//  void CockburnLimiter1d(TPZAdmChunkVector<TPZCompEl *> &elvec,int var);
//  void CockburnLimiter2d(TPZAdmChunkVector<TPZCompEl *> &elvec,int var);

    virtual void GetOrder(std::ifstream &input);

  /**Compute the minime module of the three real values*/
  REAL MinMod(REAL first,REAL second,REAL third);

};

#endif
