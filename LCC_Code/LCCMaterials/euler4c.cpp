
#include "euler4c.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include <math.h>
#include <stdio.h>

/*******   Equacao de Euler Bi-dimensional tratada em uma dimensao   ******/

int TEulerLaw4C::VariableIndex(const std::string &name) {
  if(name == "density") return 0;
  else if(name == "velocity_x") return 5;
  else if(name == "velocity_y") return 6;
  else if(name == "energy") return 7;
  else if(name == "pression") return 4;
  else if(name == "dens_velocity_x") return 1;
  else if(name == "dens_velocity_y") return 2;
  else if(name == "dens_energy") return 3;
  return TPZMaterial::VariableIndex(name);
}

void TEulerLaw4C::VariablesName(TPZVec<std::string> &names) {
  int i;
  for(i=0;i<4;i++) {
    switch(i) {
    case 0:
      names[0] = "density";
      break;
    case 1:
      names[1] = "velocity_x";
      break;
    case 2:
      names[2] = "velocity_y";
      break;
    case 3:
      names[3] = "energy";
      break;
    case 4:
      names[4] = "pression";
      break;
    case 5:
      names[5] = "dens_velocity_x";
      break;
    case 6:
      names[6] = "dens_velocity_y";
      break;
    case 7:
      names[7] = "dens_energy";
      break;
    default:
      names[i] = "state_nonumber";
      break;
    }
  }
}

int TEulerLaw4C::NSolutionVariables(int index) {
  if(index<0) {
    PZError << "TEulerLaw4C::NSolutionVariables. Bad parameter index.\n";
    return -1;
  }
  if(index < 2*NStateVariables()) return 1;
  return TPZMaterial::NSolutionVariables(index);
}

/**returns the solution associated with the var index based on the finite element approximation*/
void TEulerLaw4C::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,
                 int var,TPZVec<STATE> &Solout){
  switch(var) {
  case 0:
    Solout.Resize(1);
    Solout[0] = Sol[0];
    return;
  case 5:
    Solout.Resize(1);
    Solout[0] = Sol[1]/Sol[0];
    return;
  case 6:
    Solout.Resize(1);
    Solout[0] = Sol[2]/Sol[0];
    return;
  case 7:
    Solout.Resize(1);
    Solout[0] = Sol[3]/Sol[0];
    return;
  case 4:
    Solout.Resize(1);
		Solout[0] = Pression(Sol);
    return;
  case 1:
    Solout.Resize(1);
    Solout[0] = Sol[1];
    return;
  case 2:
    Solout.Resize(1);
    Solout[0] = Sol[2];
    return;
  case 3:
    Solout.Resize(1);
    Solout[0] = Sol[3];
    return;
  default:
    PZError << "TEulerLaw4C::Solution. Bad Parameter var.\n";
  }
  TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

STATE TEulerLaw4C::Pression(TPZVec<STATE> &U) {
  REAL velocity = 0.;
  if(!IsZero(U[0])) velocity = sqrt(U[1]*U[1]+U[2]*U[2])/U[0];
  return (fGamma-1.)*(U[3]-(.5*U[0]*velocity*velocity));
}

void TEulerLaw4C::FunctionFlux(TPZVec<STATE> &Ui,TPZFMatrix<STATE> &flux) {
    STATE pressao, velX=0., velY;
	if(IsZero(Ui[0])) {
		if(!IsZero(Ui[1])) {
		    printf("\nERRO : Densidade nula e momento em X nao nulo\n");
		    exit(1);
		}
		if(!IsZero(Ui[2])) {
		    printf("\nERRO : Densidade nula e momento em X nao nulo\n");
		    exit(1);
		}
		pressao=(fGamma-1.)*Ui[3];
	}
	else {
		velX=Ui[1]/Ui[0];
		velY=Ui[2]/Ui[0];
		pressao=(Ui[1]*velX)+(Ui[2]*velY);
		pressao=Ui[3]-(.5*pressao);
		pressao*=(fGamma-1.);
	}
/*	Verificacao para existencia das caracteristicas
	if((2*Ui[0]*Ui[2])<(Ui[1]*Ui[1])) {
		printf("\n ERRO, Na determinacao das caracteristicas (f[u])\n");
		exit(1);
	} */
	flux(0,0)=Ui[1];
	flux(0,1)=Ui[1]*velX;
	flux(0,1)+=pressao;
	flux(0,2)=Ui[2]*velX;
	flux(0,3)=Ui[3]+pressao;
	flux(0,3)*=velX;
}

void TEulerLaw4C::JacobFlux(TPZVec<STATE> &U,TPZFMatrix<STATE> &jacob) {
  double vel,auxiliar;
  if(IsZero(U[0])) {
    if(!IsZero(U[1])) {
      printf("\nERRO : Densidade NULA. Falla jacobiano.\n");
      exit(1);
    }
    vel=0.;
    if(!IsZero(U[2])) {
      printf("\nERRO : Densidade NULA e energia nao nula.\n");
      exit(1);
    }
    auxiliar=0.;
  }
  else {
    vel=U[1]/U[0];
    auxiliar=U[2]/U[0];
  }

  jacob(0,0)=jacob(0,2)=0.;
  jacob(0,1)=1.;
  jacob(1,2)=fGamma-1.;
  jacob(1,0)=.5*(fGamma-3.)*vel*vel;
  jacob(1,1)=(3.-fGamma)*vel;
  jacob(2,0)=(-1.*fGamma*vel*auxiliar)+((fGamma-1.)*vel*vel*vel);
  jacob(2,1)=(fGamma*auxiliar)-(1.5*(fGamma-1.)*vel*vel);
  jacob(2,2)=fGamma*vel;
}
void TEulerLaw4C::JacobFlux(TPZVec<STATE> &U,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal) {
  double vel,auxiliar;
  if(IsZero(U[0])) {
    if(!IsZero(U[1])) {
      printf("\nERRO : Densidade NULA. Falla jacobiano.\n");
      exit(1);
    }
    vel=0.;
    if(!IsZero(U[2])) {
      printf("\nERRO : Densidade NULA e energia nao nula.\n");
      exit(1);
    }
    auxiliar=0.;
  }
  else {
    vel=U[1]/U[0];
    auxiliar=U[2]/U[0];
  }

  jacob(0,0)=jacob(0,2)=0.;
  jacob(0,1)=1.;
  jacob(1,2)=fGamma-1.;
  jacob(1,0)=.5*(fGamma-3.)*vel*vel;
  jacob(1,1)=(3.-fGamma)*vel;
  jacob(2,0)=(-1.*fGamma*vel*auxiliar)+((fGamma-1.)*vel*vel*vel);
  jacob(2,1)=(fGamma*auxiliar)-(1.5*(fGamma-1.)*vel*vel);
  jacob(2,2)=fGamma*vel;
	jacob *= normal[0];
}

void TEulerLaw4C::ValJacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &valjacob) {
}

void TEulerLaw4C::ValRoeMatrix(TPZVec<STATE> &u,TPZVec<STATE> &Up1,TPZFMatrix<STATE> &ValRoe) {
}
void TEulerLaw4C::RoeMatrix(TPZVec<STATE> &U,TPZVec<STATE> &Up1,TPZFMatrix<STATE> &Roe) {
}

void TEulerLaw4C::EigRoeMatrix(TPZVec<STATE> &U,TPZVec<STATE> &Up1,TPZVec<STATE> &Roe) {
	STATE velprom,entalprom,cprom,raizcprom,dif;
	STATE densi, densip1, vel, velp1, ental, entalp1;
	dif=fGamma-1.;
	densi=U[0];
	densip1=Up1[0];
	vel=U[1]/densi;
	velp1=Up1[1]/densip1;
	if((densi<0.) || (densip1<0.)) {
		printf("\nERRO : Densidade negativa\n");
		exit(1);
	}
	//determinando os valores de entalpia a esquerda e direita
	ental=(-.5)*dif*vel*vel;
	ental+=((fGamma*U[2])/densi);
	entalp1=(-.5)*dif*velp1*velp1;
	entalp1+=((fGamma*Up1[2])/densip1);
	densi=sqrt(densi);
	densip1=sqrt(densip1);
	entalprom=densip1*entalp1;
	entalprom+=(densi*ental);
	ental=densi+densip1;
	entalprom/=ental;
	velprom=densip1*velp1;
	velprom+=(densi*vel);
	velprom/=ental;
	cprom=(-.5)*velprom*velprom;
	cprom+=entalprom;
	if(cprom<0) {
		printf("\nERRO : Nao caracteristicas na matriz de Roe\n");
		exit(1);
	}
	cprom*=dif;
	raizcprom=sqrt(cprom);
	//armazenando os autovalores e autovetores da matriz de Roe por linha
	Roe[0]=Roe[6]=(velprom-raizcprom);
	Roe[1]=Roe[7]=velprom;
	Roe[2]=Roe[8]=(velprom+raizcprom);
	Roe[3]=Roe[4]=Roe[5]=1.;
	Roe[9]=(entalprom-(velprom*raizcprom));
	Roe[10]=(.5*velprom*velprom);
	Roe[11]=(entalprom+(velprom*raizcprom));
	//para armazenar os elementos da matriz inversa de R
	densi=1./raizcprom;
	densip1=0.5*dif*densi*densi;
	vel=densip1*velprom;
	velp1=0.5*vel*velprom;
	densi*=0.5;
	ental=densi*velprom;
	Roe[12]=velp1+ental;
	Roe[13]=(-1.)*(vel+densi);
	Roe[14]=Roe[20]=densip1;
	Roe[15]=1.-(2.*velp1);
	Roe[16]=2*vel;
	Roe[17]=(-2.)*densip1;
	Roe[18]=velp1-ental;
	Roe[19]=densi-vel;
}

STATE TEulerLaw4C::MaxEigJacob(TPZVec<STATE> &U,TPZVec<REAL> &/*normal*/) {
	STATE vel;
  if(IsZero(U[0])) return 0.;
	vel=U[1]/U[0];              //Que acontece se Ui[0]= 0.
	STATE velson;
	velson=U[3]/U[0];
	velson=velson-(.5*vel*vel);
	velson*=(fGamma*(fGamma-1.));
	if(velson<0.) {
//		printf("\nERRO, na determinacao de c en MaxEigValJacob.\n");
//		exit(1);
    return 0.;
	}
	velson=sqrt(velson);
	return Max(fabs(vel-velson),fabs(vel+velson));
}
STATE TEulerLaw4C::ValEigJacob(TPZVec<STATE> &u,int order,int dim) {
#ifndef NOTDEBUG
  if(dim!=1) PZError << "TEulerLaw4C::ValEigJacob. Bad parameter dim.\n";
#endif
  STATE vel;
	if(IsZero(u[0])) return 0.;
  vel=u[1]/u[0];              //Que acontece se Ui[0]= 0.
  if(!order) return vel;
  STATE velson;
  velson=u[2]/u[0];
  velson=velson-(.5*vel*vel);
  velson*=(fGamma*(fGamma-1.));
  if(velson<0.) return 0.;
  velson=sqrt(velson);
  if(order==1) return vel+velson;
  return vel-velson;
}

void TEulerLaw4C::SetData(std::istream &input) {
  TConservationLaw::SetData(input);
}

void TEulerLaw4C::Print(std::ostream &out) {
}
