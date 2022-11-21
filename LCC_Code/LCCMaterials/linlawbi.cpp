/*******       File : LinLawBi.c       *******/

#include "linlawbi.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include <fstream>
#include "commonJC.h"
#include <math.h>

TLinearLaw2D::TLinearLaw2D(int id,int order) : TLinearLaw(id,order) {
  fName = "Linear conservation law 2D.";
}
TLinearLaw2D::TLinearLaw2D(int id,int order,int type) :
	TLinearLaw(id,order,type) {
}

/* 0 -> condicao fronteira periodica, -n depende do valor do ponto x*/
int TLinearLaw2D::IdBC(double *x) {
  if(x[0]<0.5) return -1;
  else if(x[1] < 1.) return -2;
  return -1;
}
void TLinearLaw2D::FunctionFlux(TPZVec<STATE> &Ui,TPZFMatrix<STATE> &funcao) {
  int i, j;
    funcao.Zero();
  for(i=0;i<fOrder;i++) {
    for(j=0;j<fOrder;j++) {
      funcao(i,0) += (fA(i,j)*Ui[j]);
      funcao(i,1) += (fB(i,j)*Ui[j]);
    }
  }
}

void TLinearLaw2D::JacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &jacob) {
  int i,j;
  for(i=0;i<2*fOrder;i++)
    for(j=0;j<fOrder;j++)
      jacob(i,j)=fA(i,j);
}
void TLinearLaw2D::JacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal) {
  int i,j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      jacob(i,j)=normal[0]*fA(i,j)+normal[1]*fB(i,j);
}

STATE TLinearLaw2D::MaxEigJacob(TPZVec<STATE> &/*U*/,TPZVec<REAL> &/*normal*/) {
  return Max(fabs(fMaxEigVal),fabs(fMaxEigValB));
}

STATE TLinearLaw2D::ValEigJacob(TPZVec<STATE> &u,int order,int dim) {
#ifndef NOTDEBUG
  if(dim!=1 && dim!=2)
    PZError << "TLinearLaw2D::ValEigJacob Bad parameter dim.\n";
#endif
  return fEigVal[(dim-1)*fOrder+order];
}

void TLinearLaw2D::Print(std::ostream &out) {
  int i,j;
  out << "\nMatrix A :";
  for(i=0;i<fOrder;i++) {
    out << std::endl;
    for(j=0;j<fOrder;j++)
      out << fA(i,j) << "\t";
  }
  out << "\nAuto-Valores de A :\n";
  for(i=0;i<fOrder;i++)
    out << fEigVal[i] << "\t";
  out << "\nMatrix dos auto-vetores direitos de A :";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fEigVect(i,j) << "\t";
  }
  out << "\nInversa da matriz dos auto-vetores de A :";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fInvEigVect(i,j) << "\t";
  }
  out << "\nMaximo autovalor = %lf\n" << fMaxEigVal;
  out << "\nMatrix |A| :";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fValA(i,j) << "\t";
  }
  out << "\nMatrix B :";
  for(i=0;i<fOrder;i++) {
    out << std::endl;
    for(j=0;j<fOrder;j++)
      out << fB(i,j) << "\t";
  }
  out << "\nAuto-Valores de B :\n";
  for(i=0;i<fOrder;i++)
    out << fEigValB[i] << "\t";
  out << "\nMatrix dos auto-vetores direitos de B :";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fEigVectB(i,j) << "\t";
  }
  out << "\nInversa da matriz dos auto-vetores de B :";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fInvEigVectB(i,j) << "\t";
  }
  out << "\nMaximo autovalor de B = %lf\n" << fMaxEigValB;
  out << "\nMatrix |B| :";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fValB(i,j) << "\t";
  }
}

void TLinearLaw2D::ValJacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &valjacob) {
  int i,j,m,index;
  double aux,aux1;
    for(i=0;i<fOrder;i++) {
    for(j=0;j<fOrder;j++) {
      aux=0.;
      for(m=0;m<fOrder;m++) {
	     aux1=fabs(fEigVal[m])*fInvEigVect(m,j);
	     aux1*=fEigVect(i,m);
	     aux+=aux1;
      }
      valjacob(i,j)=aux;
    }
  }
  for(i=0;i<fOrder;i++) {
    index=i*fOrder;
    for(j=0;j<fOrder;j++) {
      aux=0.;
      for(m=0;m<fOrder;m++) {
	     aux1=fabs(fEigValB[m])*fInvEigVectB(m,j);
	     aux1*=fEigVectB(i,m);
	     aux+=aux1;
      }
      valjacob(fOrder+i,j)=aux;
    }
  }
}
/** Revisar esta muito ruin */
int TLinearLaw2D::ValJacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal) {
	return 1;
  int i,j,m,index;
  double aux,aux1;
  for(i=0;i<fOrder;i++) {
    index=i*fOrder;
    for(j=0;j<fOrder;j++) {
      aux=0.;
      for(m=0;m<fOrder;m++) {
	     aux1=fabs(normal[0]*fEigVal[m])*fInvEigVect(m,j);
	     aux1*=fEigVectB(i,m);
	     aux+=aux1;
      }
      valjacob(i,j)=aux;
    }
  }
  for(i=0;i<fOrder;i++) {
    index=i*fOrder;
    for(j=0;j<fOrder;j++) {
      aux=0.;
      for(m=0;m<fOrder;m++) {
	     aux1=fabs(normal[1]*fEigValB[m])*fInvEigVectB(m,j);
	     aux1*=fEigVectB(i,m);
	     aux+=aux1;
      }
      valjacob(i,j)=aux;
    }
  }
	return 0;
}

void TLinearLaw2D::SetData(std::istream &input) {
  TConservationLaw::SetData(input);
  SetMatrix(input);
}

void TLinearLaw2D::SetMatrix(std::istream &input) {
  int i,j;
  double temp;
  fMaxEigVal=0.;
  fMaxEigValB = 0.;
  //Limpando para armazenar
    fA.Zero();
    fB.Zero();
    fEigVect.Zero();
    fValA.Zero();
    fInvEigVect.Zero();
  fEigVectB.Zero();
  fValB.Zero();
  fInvEigVectB.Zero();
    fEigVal.Fill(0.);
    fEigValB.Fill(0.);
  //Ingressa dados da matriz A
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
        input >> fA(i,j);
  //Ingressa dados da matriz A
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
        input >> fB(i,j);
  //Ingresso dos autovalores, autovetores e a inversa dos autovetores
  for(i=0;i<fOrder;i++) {
    input >> temp;
    fMaxEigVal=(fMaxEigVal<fabs(temp))?fabs(temp):fMaxEigVal;
    fEigVal[i]=temp;
  }
  for(i=0;i<fOrder;i++) {
    input >> temp;
    fMaxEigValB=(fMaxEigValB<fabs(temp))?fabs(temp):fMaxEigValB;
    fEigValB[i]=temp;
  }
  //Ingressa dados de los autovectores y la inversa de los autovectores
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
        input >> fEigVect(i,j);
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
        input >> fEigVectB(i,j);
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
        input >> fInvEigVect(i,j);
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
        input >> fInvEigVectB(i,j);

  std::cout << "\nVerifique :\n Matrix A \t EigVector\t InvEigVect\t EigValue";
    fA.Print(std::cout);
  fEigVect.Print(std::cout);
  fInvEigVect.Print(std::cout);
    for(i=0;i<fOrder;i++)
        std::cout << "\t" << fEigVal[i];
  std::cout << "\n";
  std::cout << "\nVerifique :\n Matrix B \t EigVector\t InvEigVect\t EigValue";
    fB.Print(std::cout);
    fEigVectB.Print(std::cout);
    fInvEigVectB.Print(std::cout);
    for(i=0;i<fOrder;i++)
        std::cout << "\t" << fEigValB[i];
  std::cout << "\n";
}

void TLinearLaw2D::SetMatrix(TPZFMatrix<STATE> &A,TPZFMatrix<STATE> &B,TPZFMatrix<STATE> &EigVect,TPZFMatrix<STATE> &EigVectB,TPZFMatrix<STATE> &InvEigVect,TPZFMatrix<STATE> &InvEigVectB,
			   TPZVec<STATE> &EigVal,TPZVec<STATE> &EigValB) {
  int i;
  double temp;
  fMaxEigVal=0.;
  fMaxEigValB = 0.;
    
    fA = A;
    fB = B;
    fEigVect = EigVect;
    fEigVectB = EigVectB;
    fInvEigVect = InvEigVect;
    fInvEigVectB = InvEigVectB;
    
  //Ingresso dos autovalores, autovetores e a inversa dos autovetores
  for(i=0;i<fOrder;i++) {
    temp = EigVal[i];
    fMaxEigVal=(fMaxEigVal<fabs(temp))?fabs(temp):fMaxEigVal;
    fEigVal[i]=temp;
  }
  for(i=0;i<fOrder;i++) {
    temp = EigValB[i];
    fMaxEigValB=(fMaxEigValB < fabs(temp))?fabs(temp):fMaxEigValB;
    fEigValB[i]=temp;
  }
}

void TLinearLaw2D::RequireMatrix() {
  int i,j;
  std::cout << "\nUtilizara o teclado? (y/n) ";
  char c[5];
  std::cin >> c[0];
  if(c[0]=='y') {
    fMaxEigVal=0.;
    fMaxEigVal = 0.;
	std::cout << "\nINGRESE DADOS POR LINHA\n";
	std::cout << "1- Matrix A,\n";
    for(i=0;i<fOrder;i++)
        for(j=0;j<fOrder;j++)
            std::cin >> fA(i,j);
	std::cout << "2- Autovalores de A,\n";
    for(i=0;i<fOrder;i++) {
		std::cin >> fEigVal[i];
      fMaxEigVal=(fMaxEigVal < fabs(fEigVal[i]))?fabs(fEigVal[i]):fMaxEigVal;
    }
	std::cout << "3- Matriz R dos autovetores de A,\n";
    for(i=0;i<fOrder;i++)
        for(j=0;j<fOrder;j++)
            std::cin >> fEigVect(i,j);
	std::cout << "4- Matriz inversa de R,\n";
    for(i=0;i<fOrder;i++)
        for(j=0;j<fOrder;j++)
            std::cin >> fInvEigVect(i,j);
	std::cout << "5- Matrix B,\n";
    for(i=0;i<fOrder;i++)
        for(j=0;j<fOrder;j++)
            std::cin >> fB(i,j);
	std::cout << "6- Autovalores de B,\n";
    for(i=0;i<fOrder;i++) {
		std::cin >> fEigValB[i];
      fMaxEigValB=(fMaxEigValB < fabs(fEigValB[i]))?fabs(fEigValB[i]):fMaxEigValB;
    }
	std::cout << "7- Matriz R dos autovetores de B,\n";
    for(i=0;i<fOrder;i++)
        for(j=0;j<fOrder;j++)
            std::cin >> fEigVect(i,j);
	std::cout << "8- Matriz inversa de R,\n";
    for(i=0;i<fOrder;i++)
        for(j=0;j<fOrder;j++)
            std::cin >> fInvEigVectB(i,j);
	std::cout << "Dados Completos.\n";
  }
  else {
	  std::cout << "\nNome do arquivo com dados relacionados a matriz A e B ";
    char Amatrix[20];
	std::cin >> Amatrix;
	std::ifstream in(Amatrix);
    SetMatrix(in);
    in.close();
  }
}

// *************************   TLinearLaw2DCircle   ***************************

TLinearLaw2DCircle::TLinearLaw2DCircle(int id) : TLinearLaw2D(id,1) {
  fName = "Linear conservation law 2D - Circular.";
}
TLinearLaw2DCircle::TLinearLaw2DCircle(int id,REAL angle,int type) :
	TLinearLaw2D(id,0,type) {
  fOrder = 1;
  SetData(angle);
}

inline void TLinearLaw2DCircle::SetData(std::istream &input) {
  TConservationLaw::SetData(input);
  GetDataCommented(input,fAngle);
  SetData(fAngle);
}

void TLinearLaw2DCircle::SetData(REAL angle) {
  fAngle = (angle*M_PI)/180.;
  fMaxEigVal=0.;
  fMaxEigValB = 0.;
  //Limpando para armazenar
    fA.Zero();
    fB.Zero();
    fEigVect.Zero();
    fValA.Zero();
    fInvEigVect.Zero();
  fEigVectB.Zero();
  fValB.Zero();
  fInvEigVectB.Zero();
    fEigVal.Fill(0.);
    fEigValB.Fill(0.);
  //Ingressa dados da matriz A
  fEigVal[0] = fA(0,0) = cos(fAngle);
  fEigVal[1] = fB(0,0) = sin(fAngle);
  //Ingresso dos autovalores, autovetores e a inversa dos autovetores
  fMaxEigVal = fA(0,0);
  fMaxEigValB = fB(0,0);
  fEigVect(0.0) = fEigVectB(0,0) = 1.;
  fInvEigVect(0,0) = fInvEigVectB(0,0) = 1.;

  std::cout << "\nVerifique :\n Matrix A \t EigVector\t InvEigVect\t EigValue";
    fA.Print(std::cout);
    fEigVect.Print(std::cout);
    fInvEigVect.Print(std::cout);
  std::cout << "\n";
  std::cout << "\nVerifique :\n Matrix B \t EigVector\t InvEigVect\t EigValue";
  fB.Print(std::cout);
  fEigVectB.Print(std::cout);
  fInvEigVectB.Print(std::cout);
  std::cout << "\n";
}

/*
void TLinearLaw2D::EigRoeMatrix(TPZVec<REAL> &Ul,TPZVec<REAL> &Ur,TPZVec<REAL> &EigRoe) {
  int i,j;
  for(i=0;i<fOrder;i++) EigRoe[i]=fEigVal[i];
  for(i=0;i<fOrder*fOrder;i++) EigRoe[fOrder+i]=fEigVect[i];
  j=fOrder*(1+fOrder);
  for(i=0;i<fOrder*fOrder;i++) EigRoe[j+i]=fInvEigVect[i];
}

void TLinearLaw2D::RoeMatrix(TPZVec<REAL> &Ul,TPZVec<REAL> &Ur,TPZFMatrix &Roe) {
  int i,j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      Roe(i,j)=fA[i*fOrder+j];
}

void TLinearLaw2D::ValRoeMatrix(TPZVec<REAL> &Ul,TPZVec<REAL> &Ur,TPZFMatrix &ValRoe) {
  int i,j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      ValRoe(i,j)=fValA[i*fOrder+j];
}
*/


