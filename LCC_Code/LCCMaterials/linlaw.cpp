#include <stdio.h>
#include "linlaw.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include <fstream>
#include "commonJC.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"

TLinearLaw::TLinearLaw(int id,int order) : TConservationLaw(id) {
    fName = "Linear Conservation Law 1D";
  while(order < 1) {
    PZError << "TLinearLaw constructor has bad order " << order << std::endl;
    PZError << "New order (Quit=0) ";
	std::cin >> order;
    if(!order) return;
  }
  fOrder = order;
}
TLinearLaw::TLinearLaw(int id,int order,int type) : TConservationLaw(id) {
  while(order < 0) {
    PZError << "TLinearLaw constructor has bad order " << order << std::endl;
    PZError << "New order (Quit=0) ";
	std::cin >> order;
    if(!order) return;
  }
  if(order>0) {
    fOrder = order;
    RequireMatrix();
  }
  else fOrder = 1;
}

int TLinearLaw::IdBC(double *x) { 
//	if(x[0] < .5) 
		return -4;
//	else 
//		return -2;
}

/*void TLinearLaw::IncrementDiffusion(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix &dsol,
        REAL weight,REAL area,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef) {
  int phr = phi.Rows();
  int i,j,k,p,c, nvar = NStateVariables(), dim=Dimension();
  if(dim==2) area = sqrt(area);
  TPZFMatrix dphixy(dim,phr,0.);
  for(i=0;i<phr;i++)
    for(j=0;j<dim;j++)
      for(k=0;k<dim;k++) dphixy(j,i) += dphi(k,i)*axes(k,j);

  int dsolr = dsol.Rows(), dsolc = dsol.Cols();
  REAL module;

  TPZVec<REAL> divflux(nvar,0.);
  TPZVec<REAL> Beta(dim*nvar,0.);   // Beta = Jacob(F(U)) em soma
  TPZVec<REAL> diffvec(nvar,0.);
  TPZFMatrix jacob(nvar*dim,nvar);
  JacobFlux(sol,jacob);
  for(i=0;i<dim;i++) {
    p = i*nvar;
    for(j=0;j<nvar;j++) {
      Beta[p+j] = ValEigJacob(sol,j,i+1);
      for(k=0;k<nvar;k++)
        divflux[j] += jacob(j+p,k)*dsol(i,k);
    }
  }
  for(i=0;i<phr;i++) {
    p = i*nvar;
    for(j=0;j<nvar;j++) diffvec[j] = 0.;
    for(c=0;c<dim;c++)
      for(j=0;j<nvar;j++)
        diffvec[j] += Beta[c*nvar+j]*dphixy(c,i);
    for(j=0;j<nvar;j++) {
      module = 0.;
      for(c=0;c<dim;c++) module += Beta[c*nvar+j]*Beta[c*nvar+j];
      module = sqrt(module);
      if(IsZero(module)) break;
//      for(c=0;c<phr;c++)
//        ek(p+j,c*nvar+j) += weight*fFactor*area*(diffvec[j]/module)*phi(c,0);
      ef(p+j,0) -= weight*fFactor*area*divflux[j]*(diffvec[j]/module);
    }
  }
}*/

void TLinearLaw::FunctionFlux(TPZVec<STATE> &Ui,TPZFMatrix<STATE> &funcao) {
  int i,j;
    funcao.Zero();
    for(i=0;i<fOrder;i++)
        for(j=0;j<fOrder;j++)
            funcao(i,0) += (fA(i,j)*Ui[j]);
}

int TLinearLaw::VariableIndex(const std::string &name) {
  int index = name.size();
  if(index<fOrder) return index;
  return TPZMaterial::VariableIndex(name);
}

void TLinearLaw::VariablesName(TPZVec<std::string> &names) {
  int i, n = AllVariablesPost();
  for(i=0;i<n;i++) {
    switch(i) {
    case 0:
      names[0] = "state_1";
      break;
    case 1:
      names[1] = "state_2";
      break;
    case 2:
      names[2] = "state_3";
      break;
    default:
      names[i] = "state_nonumber";
      break;
    }
  }
}

int TLinearLaw::NSolutionVariables(int index) {
  if(index < AllVariablesImplemented()) return 1;
  return TPZMaterial::NSolutionVariables(index);
}

/**returns the solution associated with the var index based on the finite element approximation*/
void TLinearLaw::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,
                 int var,TPZVec<STATE> &Solout){
  if(var<fOrder) {
    Solout.Resize(1);
    Solout[0] = Sol[var];
    return;
  }
  TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

void TLinearLaw::ApplyPeriodic(TPZCompMesh &cmesh) {
  int i, nelem = cmesh.NElements();
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  TPZCompEl *cel;
  TPZInterpolatedElement *el;
	int seqnum1, seqnum2;

  for(i=0;i<nelem;i++) {
    cel = elvec[i];
    if(!cel || cel->IsInterface()) continue;
    el = (TPZInterpolatedElement *)cel;
    if(el->Material()->Id()>-1) continue;
    if(IsZero(el->Reference()->NodePtr(0)->Coord(0)-1.)) {
 //     if(el->IsConnectContinuous(0))
   //     seqnum2 = el->Connect(0).SequenceNumber();
     // else 
		seqnum2 = el->Connect(1).SequenceNumber();
    }
    else {
     // if(el->IsConnectContinuous(0))
       // seqnum1 = el->Connect(0).SequenceNumber();
     // else 
		seqnum1 = el->Connect(1).SequenceNumber();
    }
  }
  int size = cmesh.Block().Size(seqnum2);
  for(i=0;i<size;i++)
    cmesh.Block()(seqnum1,0,i,0) = cmesh.Block()(seqnum2,0,i,0);
}

void TLinearLaw::Print(std::ostream &out) {
  int i,j;
  //	PrintData1D();
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
  out << "\nInversa da matriz dos auto-vetores de A:";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fInvEigVect(i,j) << "\t";
  }
  out << "\nMaximo autovalor = \n" << fMaxEigVal;
  out << "\nMatrix |A| :";
  for(i=0;i<fOrder;i++) {
    out << "\n";
    for(j=0;j<fOrder;j++)
      out << fValA(i,j) << "\t";
  }
}

void TLinearLaw::JacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &jacob) {
  int i,j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      jacob(i,j)=fA(i,j);
}
void TLinearLaw::JacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal) {
  int i,j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      jacob(i,j) = normal[0]*fA(i,j);
}

STATE TLinearLaw::MaxEigJacob(TPZVec<STATE> &/*U*/,TPZVec<REAL> &/*normal*/) {
  return fMaxEigVal;
}
STATE TLinearLaw::ValEigJacob(TPZVec<STATE> &u,int order,int dim) {
#ifndef NOTDEBUG
  if(dim!=1) PZError << "TLinearLaw::ValEigJacob Bad parameter dim.\n";
#endif
  return fEigVal[order];
}

void TLinearLaw::ValJacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &valjacob) {
  int i,j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      valjacob(i,j)=fValA(i,j);
}
int TLinearLaw::ValJacobFlux(TPZVec<STATE> &/*Ui*/,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal) {
  int i,j;
	if(fOrder>4) return 1;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      valjacob(i,j)=fabs(normal[0])*fValA(i,j);
	return 0;
}

void TLinearLaw::InverseJacob(TPZFMatrix<STATE> &jac) {
	REAL a = jac(0,0);
	REAL det, b, c, d, e, f, g, h, i;
	REAL EIHF, FGDI, HDGE;
	switch(fOrder) {
	case 1:
		jac(0,0) = 1./a;
		break;
	case 2:
		det = a*jac(1,1)-jac(0,1)*jac(1,0);
		det = 1./det;
	  if(IsZero(det)) {
  	  PZError << "TLinearLaw::InverseJacob. Determinante zero, matriz singular.\n";
	  	exit(1);
	  }
		jac(0,0) = jac(1,1)*det;
		jac(1,1) = a*det;
		det *= -1.;
		jac(0,1) *= det;
		jac(1,0) *= det;
		break;
	case 3:
		b=jac.GetVal(0,1); c=jac.GetVal(0,2);
		d=jac.GetVal(1,0); e=jac.GetVal(1,1); f=jac.GetVal(1,2);
		g=jac.GetVal(2,0); h=jac.GetVal(2,1); i=jac.GetVal(2,2);
	  EIHF = e*i-h*f; FGDI = f*g-d*i; HDGE = h*d-g*e;
  	det = a*EIHF + b*FGDI + c*HDGE;
	  if(IsZero(det)) {
  	  PZError << "TLinearLaw::InverseJacob. Determinante zero, matriz singular.\n";
	  	exit(1);
	  }
  
  	det = 1./det;
	  jac(0,0) = det*EIHF;
  	jac(1,0) = det*FGDI;
	  jac(2,0) = det*HDGE;
	  jac(0,1) = det*(c*h-b*i);
  	jac(1,1) = det*(a*i-g*c);
	  jac(2,1) = det*(g*b-a*h);
  	jac(0,2) = det*(b*f-e*c);
	  jac(1,2) = det*(d*c-a*f);
	  jac(2,2) = det*(a*e-d*b);
		break;
	case 4:
	{
		b=jac.GetVal(0,1); c=jac.GetVal(0,2); d=jac.GetVal(0,3);
		e=jac.GetVal(1,0); f=jac.GetVal(1,1); g=jac.GetVal(1,2); h=jac.GetVal(1,3);
		i=jac.GetVal(2,0);
		REAL j=jac.GetVal(2,1), k=jac.GetVal(2,2), l=jac.GetVal(2,3);
		REAL m=jac.GetVal(3,0), n=jac.GetVal(3,1), r=jac.GetVal(3,2), s=jac.GetVal(3,3);
	  REAL BHDF = b*h-d*f, BKCJ = b*k-c*j, BSDN = b*s-d*n, CFBG = c*f-b*g;
		REAL CLDK = c*l-d*k, CNBR = c*n-b*r, DGCH = d*g-c*h, DJBL = d*j-b*l;
		REAL DRCS = d*r-c*s, FRGN = f*r-g*n, GJFK = g*j-f*k, GSHR = g*s-h*r;
		REAL FSHN = f*s-h*n, LNJS = l*n-j*s, HJFL = h*j-f*l, HKGL = h*k-g*l;
		REAL KNJR = k*n-j*r, LRKS = l*r-k*s;
	  det = DJBL*(g*m-e*r)+BHDF*(k*m-i*r)+BSDN*(g*i-e*k);
  	det += (FSHN*(a*k-c*i)+LNJS*(a*g-c*e)+HJFL*(a*r-c*m));
	 	if(IsZero(det)) {
  		PZError << "TEulerLaw2D::InverseJacob. Determinante zero, matriz singular.\n";
	 		exit(1);
	  }
  
  	det = 1./det;
	  jac(0,0) = det*(k*FSHN+g*LNJS+r*HJFL);
  	jac(1,0) = det*(m*HKGL+i*GSHR+e*LRKS);
	  jac(2,0) = (-det)*(m*HJFL+i*FSHN+e*LNJS);
  	jac(3,0) = det*(m*GJFK+i*FRGN+e*KNJR);
	  jac(0,1) = det*(d*KNJR+b*LRKS-c*LNJS);
	  jac(1,1) = det*(m*CLDK+i*DRCS-a*LRKS);
  	jac(2,1) = det*(m*DJBL+i*BSDN+a*LNJS);
	  jac(3,1) = det*(m*BKCJ+i*CNBR-a*KNJR);
  	jac(0,2) = det*(d*FRGN-c*FSHN+b*GSHR);
	  jac(1,2) = det*(m*DGCH-e*DRCS-a*GSHR);
  	jac(2,2) = det*(m*BHDF-e*BSDN+a*FSHN);
	  jac(3,2) = det*(m*CFBG-e*CNBR-a*FRGN);
  	jac(0,3) = det*(d*GJFK-c*HJFL+b*HKGL);
	  jac(1,3) = (-det)*(i*DGCH+e*CLDK+a*HKGL);
  	jac(2,3) = det*(a*HJFL-i*BHDF-e*DJBL);
	  jac(3,3) = (-det)*(i*CFBG+e*BKCJ+a*GJFK);
		break;
	}
	default:
		jac.Identity();
	}
}

void TLinearLaw::SetData(std::istream &input) {
  TConservationLaw::SetData(input);
  SetMatrix(input);
}

void TLinearLaw::SetMatrix(std::istream &input) {
  int i,j,k,index;
  double temp;
  fMaxEigVal=0.;
  //Ingressa dados da matriz A
  GetCommentary(input,1);   // Matrix A, Eigenvalues, Eigenvectors and InverseEigenvectors
  for(i=0;i<fOrder;i++)
      for(j=0;j<fOrder;j++)
          input >> fA(i,j);
  //Ingresso dos autovalores, autovetores e a inversa dos autovetores
  for(i=0;i<fOrder;i++) {
    input >> temp;
    fMaxEigVal=(fMaxEigVal,fabs(temp))?fabs(temp):fMaxEigVal;
    fEigVal[i]=temp;
  }
  for(i=0;i<fOrder;i++)
      for(j=0;j<fOrder;j++)
          input >> fEigVect(i,j);
  for(i=0;i<fOrder;i++)
      for(j=0;j<fOrder;j++)
          input >> fInvEigVect(i,j);

  REAL prod, eigval;
  for(i=0;i<fOrder;i++) {
    prod = 0.;
    for(j=0;j<fOrder;j++) {
      for(k=0;k<fOrder;k++) {
        eigval = fEigVal[k]/fabs(fEigVal[k]);
        prod += (fEigVect(i,k)*eigval*fInvEigVect(k,j));
      }
      fInverseQrt(i,j) = prod;
    }
  }
  double aux;
  for(i=0;i<fOrder;i++) {
    for(j=0;j<fOrder;j++) {
      aux=0.;
      for(k=0;k<fOrder;k++) {
        temp = fabs(fEigVal[k])*fInvEigVect(k,j);
        temp *= fEigVect(i,k);
        aux += temp;
      }
      fValA(i,j) = aux;
    }
  }

  std::cout << "\nVerifique :\n Matrix A \t EigVector\t InvEigVect\t EigValue";
  for(i=0;i<fOrder;i++) {
    index=i*fOrder;
    for(j=0;j<fOrder;j++) {
		std::cout << "\n" << fA(i,j) << "\t" << fEigVect(i,j);
		std::cout << "\t" << fInvEigVect(i,j);
    }
	std::cout << "\t" << fEigVal[i];
  }
  std::cout << "\n";
}

void TLinearLaw::SetMatrix(TPZFMatrix<STATE> &A,TPZFMatrix<STATE> &EigVect,TPZFMatrix<STATE> &InvEigVect,
			   TPZVec<STATE> &EigVal) {
  int i,j,k;
  double temp,eigval;
  fMaxEigVal=0.;
#ifdef DEBUG
  if(fOrder!=A.Rows() || fOrder!=EigVect.Rows() || fOrder!=InvEigVect.Rows() ||
     fOrder!=EigVal.NElements()) {
    PZError << "TLinearLaw1D::SetMatrix. Some matrix is out of the dimension.\n";
    return;
  }
#endif
    fA = A;
  //Ingresso dos autovalores, autovetores e a inversa dos autovetores
  for(i=0;i<fOrder;i++) {
    temp = EigVal[i];
    fMaxEigVal=(fMaxEigVal<fabs(temp))?fabs(temp):fMaxEigVal;
    fEigVal[i]=temp;
  }
    fEigVect = EigVect;
    fInvEigVect = InvEigVect;
  /**Criando a matriz C = A*|A|¯1 */
  for(i=0;i<fOrder;i++) {
    temp = 0.;
    for(j=0;j<fOrder;j++) {
      for(k=0;k<fOrder;k++) {
        eigval = fEigVal[k]/fabs(fEigVal[k]);
        temp += (fEigVect(i,k)*eigval*fInvEigVect(k,j));
      }
      fInverseQrt(i,j) = temp;
    }
  }
  for(i=0;i<fOrder;i++) {
    for(j=0;j<fOrder;j++) {
      temp=0.;
      for(k=0;k<fOrder;k++) {
	     eigval = fabs(fEigVal[k])*fInvEigVect(k,j);
	     eigval *= fEigVect(i,k);
	     temp += eigval;
      }
      fValA(i,j) = temp;
    }
  }
}

void TLinearLaw::RequireMatrix() {
  int i,j;
  std::cout << "\nIntroduz a matriz A pelo teclado? (y/n) ";
  char c[5];
  std::cin >> c[0];
  if(c[0]=='y') {
    fMaxEigVal=0.;
	std::cout << "\nINGRESE DADOS POR LINHA\n";
	std::cout << "1- Matrix A,\n";
    for(i=0;i<fOrder;i++)
        for(j=0;j<fOrder;j++)
            std::cin >> fA(i,j);
	std::cout << "2- Autovalores de A,\n";
    for(i=0;i<fOrder;i++) {
		std::cin >> fEigVal[i];
      fMaxEigVal=(fMaxEigVal,fabs(fEigVal[i]))?fabs(fEigVal[i]):fMaxEigVal;
    }
	std::cout << "3- Matriz R dos autovetores de A,\n";
    for(i=0;i<fOrder;i++)
        for(j=0;j<fOrder;j++)
            std::cin >> fEigVect(i,j);
	std::cout << "4- Matriz inversa de R,\n";
    for(i=0;i<fOrder;i++)
        for(j=0;j<fOrder;j++)
            std::cin >> fInvEigVect(i,j);
	std::cout << "Dados Completos.\n";
  }
  else {
	  std::cout << "\nNome do arquivo com dados relacionados a matriz A ";
    char Amatrix[20];
	std::cin >> Amatrix;
	std::ifstream in(Amatrix);
    SetMatrix(in);
    in.close();
  }
  int k,index;
  REAL temp,eigval;
  for(i=0;i<fOrder;i++) {
    for(j=0;j<fOrder;j++) {
      temp=0.;
      for(k=0;k<fOrder;k++) {
	     eigval = fabs(fEigVal[k])*fInvEigVect(k,j);
	     eigval *= fEigVect(i,k);
	     temp += eigval;
      }
      fValA(i,j) = temp;
    }
  }
}

void TLinearLaw::InverseQrt(TPZVec<REAL> &sol,TPZFMatrix<STATE> &C) {
  int i, j;
  for(i=0;i<fOrder;i++)
    for(j=0;j<fOrder;j++)
      C(i,j) = fInverseQrt(i,j);
}

/**To evaluate true-error, L2-error and estimate error*/
void TLinearLaw::Errors(TPZVec<REAL> &x,TPZVec<STATE> &sol,TPZFMatrix<STATE> &dsol,TPZFMatrix<REAL> &axes,TPZVec<STATE> &flux,
     TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {

  TPZVec<STATE> udif(sol);
  int i, nvar = NStateVariables();
  for(i=0;i<nvar;i++) udif[i] -= u_exact[i];

  values.Fill(0.);

  for(i=0;i<nvar;i++) {
    values[0] += udif[i]*udif[i];
    values[1] += udif[i]*udif[i];
    values[2] += fabs(udif[i]);
  }
}
