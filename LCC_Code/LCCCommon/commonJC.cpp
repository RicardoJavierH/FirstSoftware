#include "commonJC.h"
#include "myheader.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMaterial.h"

SimulationData::SimulationData() {
	Equation = 1;
	nuzero = 1;
	typeel = 5;
	print = 0;
	nref = 1;
	NumericMethod = 0;
	pOrder = 1;
	stepsolver = 0;
	nsV = 1;
	sV.resize(1);
	sV[0] = "Solution";
	nvV = 0;
	applyadaptivity = 0;
}

// To import user data from filedata : fromgid, typeel, nref, NumericMethod, pOrder, stepsolver, nScalarVars,scalarVars,nVectorVars,vectorVars
std::ifstream *UserDataImport(std::string &Partition, std::string& fromgid, SimulationData &simulationdata) {
	size_t tt;
#ifdef WIN32
	tt = Partition.find(':');     //)('/');
	Partition.erase(tt + 1);
#else
	tt = Partition.find("Documents");
	Partition.erase(tt + 9);
#endif
	Partition += "/LCC_SimulatorIO/InputData/";
	fromgid = Partition;
	Partition += "userdata.txt";

	std::ifstream *input = new std::ifstream(Partition, std::ios::in);
	if (!input || !input->is_open())
		return NULL;
	char chardata;
	int data;
	int count = 0;
	std::string name;
	do {
		input->get(chardata);
		if (chardata == '@') {
			count++;
			(*input) >> data;
			switch (count) {
			case 1:
				if (data) {
					(*input) >> name;
					fromgid += name;
				}
				else fromgid.clear();
				name.clear();
				break;
			case 2:
					simulationdata.typeel = data;
				break;
			case 3:
					simulationdata.print = data;
				break;
			case 4:
					simulationdata.nref = data;
				break;
			case 5:
					simulationdata.Equation = data;
				break;
			case 6:
					simulationdata.nuzero = data;
                break;
			case 7:
					simulationdata.NumericMethod = data;
				break;
			case 8:
					simulationdata.pOrder = data;
				break;
			case 9:
					simulationdata.stepsolver = data;
				break;
			case 10:
			{
				simulationdata.nsV = data;
				if (simulationdata.nsV) {
					simulationdata.sV.Resize(0);
					simulationdata.sV.Resize(simulationdata.nsV);
				}
				for (int c = 0; c < simulationdata.nsV; c++) {
					(*input) >> name;
					simulationdata.sV[c] = name;
				}
			}
			break;
			case 11:
			{
				simulationdata.nvV = data;
				if (simulationdata.nvV) {
					simulationdata.vV.Resize(0);
					simulationdata.vV.Resize(simulationdata.nvV);
				}
				for (int c = 0; c < simulationdata.nvV; c++) {
					(*input) >> name;
					simulationdata.vV[c] = name;
				}
			}
			break;
			case 12:
					simulationdata.applyadaptivity = data;
					break;
			}
		}
	} while (!(input->eof()) && count<11);
	return input;
}
bool CheckTypeEl(int typeel) {
	if (typeel == 1 || typeel == 3 || typeel == 5 || typeel == 8 || typeel == 10 || typeel == 12 || typeel == 14 || typeel == 20 || typeel == 30 || typeel == 40)
		return true;
	return false;
}
void PrintInfo(TPZVec<int64_t>& indexsubelements, TPZVec<REAL>& means, TPZVec<REAL>& wavcoefs, std::string& sout) {
	std::string file = sout + "infowavelets.txt";
	std::ofstream infowavelets(file, std::ios::app);
	if (!infowavelets.is_open()) return;
	int64_t ndata;
	indexsubelements.Print(infowavelets);
	for (ndata = 0; ndata < indexsubelements.NElements(); ndata++)
		infowavelets << std::endl << indexsubelements[ndata] << "\t" << means[ndata] << "\t" << wavcoefs[ndata];
	infowavelets.close();
}
void PrintInfo(int n, TPZVec<int64_t>& indexsubelements, TPZVec<REAL>& means, std::string& sout) {
	std::string file = sout + "infoerrors.txt";
	std::ofstream infoerrors(file, std::ios::app);
	if (!infoerrors.is_open()) return;
	int64_t ndata;
	indexsubelements.Print(infoerrors);
	for (ndata = 0; ndata < indexsubelements.NElements(); ndata++) {
		infoerrors << std::endl << indexsubelements[ndata] << "\t";
		for (int nn = 0; nn < n; nn++)
			infoerrors << means[ndata + nn] << "\t";
	}
	infoerrors.close();
}

// To load refinement pattern
void LoadUniformPatterns(int typeel) {
	// TO WASTE
	if (typeel == 1) {
		gRefDBase.InitializeUniformRefPattern(EOned);
	}
	else if (typeel == 3) {
		gRefDBase.InitializeUniformRefPattern(EOned);
		gRefDBase.InitializeUniformRefPattern(ETriangle);
	}
	else if(typeel == 5) {
		gRefDBase.InitializeUniformRefPattern(EOned);
		gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	}
	else {
		gRefDBase.InitializeAllUniformRefPatterns();
	}
}

// Creating necessary directories and return name of directory to results files
int PutNameDirectoryForResults(std::string& sout) {
	struct tm* timeinfo;
	std::string command;
	time_t tempo; time(&tempo);
	timeinfo = localtime(&tempo);
#ifdef WIN32
	command = "mkdir \"";   // command = "mkdir -p \"";
#else
	command = "mkdir -p \"";   // command = "mkdir -p \"";
#endif

	size_t tt;
#ifdef WIN32
	tt = sout.find(':');     //)('/');
	sout.erase(tt + 1);
#else
	tt = sout.find("Documents");
	sout.erase(tt + 9);
#endif
	sout += "/LCC_SimulatorIO/SimulationResults/R";
	std::stringstream ssout;
	ssout << sout;
    ssout << timeinfo->tm_year + 1900 << "_" << std::setfill('0') << std::setw(2)  << timeinfo->tm_mon + 1 << "_";
    ssout << std::setfill('0') << std::setw(2) << timeinfo->tm_mday << "_" << std::setfill('0') << std::setw(2) << timeinfo->tm_hour << "_";
    ssout << std::setfill('0') << std::setw(2)<< timeinfo->tm_min << "_" << std::setfill('0') << std::setw(2) << timeinfo->tm_sec << "/";
	char lixo[256];
	ssout.getline(lixo, 256);
	sout = lixo;
	command += sout;
	command += "\"";
	int result = system(command.c_str());
	return result;
}
void ChangeInternalConnectOrder(TPZCompMesh *mesh, int typeel, int addToOrder) {

	int nEl = mesh->NElements();
	int dim = mesh->Dimension();

	for (int iel = 0; iel < nEl; iel++) {
		TPZCompEl *cel = mesh->ElementVec()[iel];
		if (!cel) continue;
		int ncon = cel->NConnects();
		int corder = 0;
		int nshape = 0;
		int nshape2 = 0;

		if (cel->Dimension() == dim)
		{
			TPZConnect &conel = cel->Connect(ncon - 1);
			corder = conel.Order();
			nshape = conel.NShape();

			int neworder = corder + addToOrder;//Aqui = +1
			int64_t cindex = cel->ConnectIndex(ncon - 1);
			conel.SetOrder(neworder, cindex);

			TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
			intel->SetPreferredOrder(neworder);
			nshape = intel->NConnectShapeF(ncon - 1, neworder);

			if (dim == 2 && addToOrder == 1)
			{
				if (typeel==3) {
					nshape2 = (corder + 2)*(corder + 2) - 1;
				}
				else {//Quadrilateral
					nshape2 = 2 * (corder + 1)*(corder + 2);
				}
				if (nshape2 != nshape)
				{
					DebugStop();
				}
			}

			conel.SetNShape(nshape);
			mesh->Block().Set(conel.SequenceNumber(), nshape);
		}
	}
	mesh->CleanUpUnconnectedNodes();
	mesh->ExpandSolution();
}


// TO GRADIENT RECONSTRUCTION
REAL MeanCell(TPZCompEl *cel, int IntOrder) {
	TPZIntPoints *pointIntRule = ((TPZInterpolatedElement*)cel)->Reference()->CreateSideIntegrationRule((cel->Reference()->NSides()) - 1, IntOrder);
	int it, npoints = pointIntRule->NPoints();
	int dim = cel->Mesh()->Dimension();
	REAL integral = 0.0;
	TPZManVector<REAL> point(3, 0.);
	TPZManVector<REAL> xpoint(3, 0.);
	TPZManVector<STATE> sol(1, 0.);
	REAL weight;
	for (it = 0; it < npoints; it++) {
		pointIntRule->Point(it, point, weight);
		if (dim == cel->Dimension())                                              /// Jorge 2019
			weight /= cel->Reference()->RefElVolume();
		cel->Reference()->X(point, xpoint);
		cel->Solution(point, 1, sol);
		integral += weight * (sol[0]);  // ExactSolution(dim, xpoint);
	}
	//REAL area = cel->Reference()->Volume();
	return integral;
}

void GradientReconstructionByLeastSquares(TPZCompEl *cel, TPZManVector<REAL, 3> &center, TPZVec<REAL> &Grad) {
	TPZFMatrix<REAL> grad;
	int dim;
	dim = cel->Mesh()->Dimension();

	// Nada sera realizado para elementos com dimensao diferente da dimensao do problema
	if (!cel)   // || cel->Dimension() != dim) 
		DebugStop();
	REAL solalfa;
	REAL solbeta;

	int k, side;
	// Integration order to calculate cell mean of solution
	int intOrder = 2;

	TPZStack<TPZCompElSide> neighs;
	int nneighs;

	center.Resize(3, 0.);
	TPZManVector<REAL> centerpsi(3, 0.0), centerbeta(3, 0.0);

	TPZFMatrix<REAL> A(dim, dim);  // Linear System matrix
	grad.Redim(dim, 1);

	// matrizes para aplicar o metodo dos minimos quadrados
	TPZFMatrix<REAL> DeltaH;
	TPZFMatrix<REAL> DeltaHTranspose;
	TPZFMatrix<REAL> DifSol;

	// Encontrando o centro do elemento atual (cel)
	TPZGeoEl* gelalfa = cel->Reference();
	gelalfa->CenterPoint(gelalfa->NSides() - 1, centerpsi);
	center.Fill(0.);
	gelalfa->X(centerpsi, center);

	solalfa = MeanCell(cel, intOrder);

	neighs.Resize(0);

	// Procuramos todos los elementos vecinos a cel (sobre todos los lados) sin duplicados
	for (side = 0; side < cel->Reference()->NSides(); side++)
	{
		TPZCompElSide celside(cel, side);
		celside.ConnectedElementList(neighs, 0, 0);
	}

	// si no hay vecinos continuamos con el siguiente elemento
	nneighs = neighs.NElements();
	if (!nneighs) DebugStop();

	std::set<TPZCompEl *> neighscel;
	for (int i = 0; i < nneighs; i++)
	{
		TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(neighs[i].Element());
		if (!InterpEl || InterpEl->Dimension() != dim) continue;
		neighscel.insert(neighs[i].Element());
	}

	nneighs = neighscel.size();

	// si hay vecinos realizamos el proceso de minimos quadrados para calcular una aproximacion del gradiente
	// Para cada vecino calculamos los deltaH (desde su centro al centro del elemento corriente)
	// y el valor de la solucion en su centro solbeta
	DeltaH.Redim(nneighs, dim);
	DeltaHTranspose.Redim(dim, nneighs);
	DifSol.Redim(nneighs, 1);

	// Montando la matriz de los deltas DeltaH y de las diferencias de las soluciones DifSol
	int ineighs = -1;
	int64_t counter = 0;
	std::set<TPZCompEl *>::iterator it;
	for (it = neighscel.begin(); it != neighscel.end(); ++it)
	{
		//(*it)->Print();
		ineighs++;
		TPZGeoEl * gelbeta = (*it)->Reference();

		if (!gelbeta) DebugStop();

		centerpsi.Fill(0.0);
		centerbeta.Fill(0.0);
		gelbeta->CenterPoint(gelbeta->NSides() - 1, centerpsi);
		gelbeta->X(centerpsi, centerbeta);
		solbeta = MeanCell(gelbeta->Reference(), intOrder);

		for (k = 0; k < dim; k++)
		{
			DeltaH(ineighs, k) = centerbeta[k] - center[k];
		}
		DifSol(ineighs, 0) = solbeta - solalfa;
		counter++;
	}

	// Resolviendo el sistema por los minimos cuadrados: DeltaH_t * DifSol = DeltaH_t * DeltaH * Grad(u)
	A.Zero();
	DeltaH.Transpose(&DeltaHTranspose);
	grad = DeltaHTranspose * DifSol;
	A = DeltaHTranspose * DeltaH;
	if (counter > 0)
		A.SolveDirect(grad, ELU);

	// Return gradient vector
	Grad.Resize(dim);
	for (int j = 0; j < dim; j++)
		Grad[j] = grad(j, 0);
}
// Generate an output of all geomesh to VTK, associating to each one the given data, creates a file with filename given
void PrintDataMeshVTK(TPZCompMesh *cmesh, std::string &filename, TPZFMatrix<REAL> &elData)
{
	std::ofstream file(filename);
#ifdef PZDEBUG
	if (!file.is_open())
		DebugStop();
#endif

	int dim = cmesh->Dimension();
	TPZGeoMesh *gmesh = cmesh->Reference();
	int64_t nelements = elData.Rows();

	std::stringstream connectivity, type, cellval1, cellval2, cellval3;

	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;

	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";

	int64_t t, c, el;
	int64_t actualNode = -1, size = 0, nVALIDelements = 0;
	int64_t counternodes = gmesh->NNodes();
	TPZGeoEl *gel;
	TPZVec<REAL> centerpsi(3);
	TPZManVector<REAL> center(3);
	TPZManVector<REAL> gradient(3);
	int64_t counter = 0;

	for (el = 0; el < nelements; el++)
	{
		gel = cmesh->ElementVec()[elData(el, 2 * dim)]->Reference();
		if (!gel)              /// Jorge 2019     --->  || gel->Reference()->Dimension() != dim)
			continue;

		MElementType elt = gel->Type();
		int elNnodes = MElementType_NNodes(elt);

		size += (1 + elNnodes);
		connectivity << elNnodes;

		for (t = 0; t < elNnodes; t++)
		{
			actualNode = gel->NodeIndex(t);
			if (actualNode < 0)
				DebugStop();

			connectivity << " " << actualNode;
		}
		connectivity << std::endl;

		int elType = TPZVTKGeoMesh::GetVTK_ElType(gel);
		type << elType << std::endl;
		REAL gradN, Norm, tempN = 0.0, temp = 0.0;
		TPZVec<REAL> v(3);
		for (c = 0; c < dim; c++) {
			v[c] = elData(counter, c);
			tempN += v[c] * v[c];
			gradient[c] = elData(counter, dim + c);
			temp += gradient[c] * gradient[c];
		}
		Norm = sqrt(tempN);
		gradN = sqrt(temp);
		for (c = 0; c < dim; c++) {
			if (IsZero(Norm)) v[c] = 0.;
			else v[c] /= Norm;
		}
		TPZVec<REAL> vort(3);
		vort[0] = -v[1];
		vort[1] = v[0];
		vort[2] = 0.0;

		cellval1 << gradN << std::endl;
		if (dim == 2) {
			cellval2 << vort[0] * gradient[0] + vort[1] * gradient[1] << std::endl;
			cellval3 << v[0] * gradient[0] + v[1] * gradient[1] << std::endl;
		}
		else {
			cellval2 << elData(counter, 1) << std::endl;
			cellval3 << elData(counter, 0) << std::endl;
		}
		counter++;
		nVALIDelements++;
	}

	// Printing all nodes of the mesh
	file << counternodes << " float" << std::endl;
	for (t = 0; t < counternodes; t++) {
		TPZGeoNode *node = &(gmesh->NodeVec()[t]);
		for (c = 0; c < 3; c++) {
			double coord = node->Coord(c);
			file << coord << " ";
		}
		file << std::endl;
	}

	file << std::endl << "CELLS " << nVALIDelements << " ";

	file << size << std::endl;
	file << connectivity.str() << std::endl;

	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str() << std::endl;

	file << "CELL_DATA" << " " << nVALIDelements << std::endl;
	file << "SCALARS Magnitude float" << std::endl;
	file << "LOOKUP_TABLE default" << std::endl;

	file << cellval1.str() << std::endl;

	file << "SCALARS DotProductSpecial float" << std::endl;
	file << "LOOKUP_TABLE default" << std::endl;

	file << cellval2.str() << std::endl;

	file << "SCALARS DotProduct float" << std::endl;
	file << "LOOKUP_TABLE default" << std::endl;

	file << cellval3.str();

	file.close();
}
void PosProcessGradientReconstruction(TPZCompMesh *cmesh, TPZFMatrix<REAL> &datagradients) {

	// Redimensionando a matriz dos dados da reconstruca de gradientes
	int dim = cmesh->Dimension();
	int64_t nelem = cmesh->NElements();
	datagradients.Redim(nelem, 2 * dim + 2);
	datagradients.Zero();

	TPZManVector<REAL, 3> center;
	TPZManVector<REAL> Grad(dim);

	int64_t i, k;
	int64_t counter = 0;

	TPZCompEl *cel;
	// Calculando el gradiente por elemento computacional
	for (i = 0; i < nelem; i++) {

		cel = cmesh->ElementVec()[i];

		// Nada sera realizado para elementos con dimension diferente de la dimension del problema
		if (!cel) // || cel->Dimension() != dim)
			continue;
		if (cel->IsInterface())
			continue;
		center.Fill(0.0);
		Grad.Fill(0.0);

		GradientReconstructionByLeastSquares(cel, center, Grad);

		//data of the vector gradiente
		for (k = 0; k < dim; k++) {
			datagradients(counter, k) = center[k];//centro do elemento
			datagradients(counter, dim + k) = Grad[k];//valor do gradiente
		}
		// Increment a last column to store volume of the finite element
		datagradients(counter, 2 * dim) = cel->Index();
		datagradients(counter, 2 * dim + 1) = cel->VolumeOfEl();

		counter++;
	}
	// Redimensionando la matriz de los gradientes
	k = datagradients.Cols();
	datagradients.Resize(counter, k);
}

void GetCommentary(std::istream &input,int nlines) {
  char lixo[256];
  for(int i=0;i<nlines;i++) {
    input >> lixo[0];
    input.getline(lixo,256);
  }
}
void GetDataCommented(std::istream &input,int &value) {
  char lixo[256];
  input >> lixo[0];
  input.getline(lixo,256);
  input >> value;
}
void GetDataCommented(std::istream &input,int const nlines,int &value) {
  char lixo[256];
  input >> lixo[0];
  for(int j=0;j<nlines;j++)
      input.getline(lixo,256);
  input >> value;
}
void GetDataCommented(std::istream *input,REAL &value) {
    char chardata;
    bool getit = false;
  do {
      input->get(chardata);
      if (chardata == '@') {
          getit = true;
          (*input) >> value;
      }
  } while (!(input->eof()) && !getit);
}
void GetDataCommented(std::istream *input,int &value) {
    char chardata;
    bool getit = false;
  do {
      input->get(chardata);
      if (chardata == '@') {
          getit = true;
          (*input) >> value;
      }
  } while (!(input->eof()) && !getit);
}
void GetDataCommented(std::istream &input,REAL &value) {
  char lixo[256];
  input >> lixo[0];
  input.getline(lixo,256);
  input >> value;
}
void GetDataCommented(std::istream &input,TPZVec<REAL> &vector) {
  char lixo[256];
  input >> lixo[0];
  input.getline(lixo,256);
  int i, n = vector.NElements();
  for(i=0;i<n;i++) input >> vector[i];
}
void GetDataCommented(std::istream &input,char *string,int size) {
  char lixo[256];
  input >> lixo[0];
  input.getline(lixo,256);
  input >> string[0];
  string[1] = '\0';
  input.getline(lixo,256);
  lixo[size - 1] = '\0';
  strncat(string,lixo,size-1);
}
void GetDataCommented(std::istream &input, std::string &str) {
    char lixo[256];
    input >> lixo[0];
    str.empty();
    input.getline(lixo+1, 256);
    str = lixo;
}

void TesteNormal(TPZCompMesh *cmesh) {
    int64_t i, nelem = cmesh->NElements();
    for(i=0;i<nelem;i++) {
        TPZCompEl *cel = cmesh->ElementVec()[i];
        if(!cel || cel->IsInterface() || !cel->Dimension())
            continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        TPZManVector<REAL> normal(3,0.);
        TPZManVector<REAL> qsi(3,0.);
        TPZStack<TPZGeoElSide> allneigh;
        for(int j=0;j<gel->NSides();j++) {
            TPZGeoElSide gside(gel,j);
            gside.AllNeighbours(allneigh);
            qsi.Resize(gside.Dimension());
            gside.Normal(qsi,gel,allneigh[0].Element(),normal);
            allneigh.Resize(0);
        }
    }
}

void InvJacob2d(TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &jacinv) {
    REAL det = 1./(jacinv(0,0)*jacinv(1,1));
    if(jacinv(0,0)<0.) det *= -1.;
    REAL coef1 = 1./(fabs(jacinv(1,1))*det), coef2 = -jacinv(0,1);

    jacinv(0,0) = coef1*axes(1,1)+coef2*axes(0,1);
    jacinv(0,1) = -(coef1*axes(1,0)+coef2*axes(0,0));
    jacinv(1,0) = -(jacinv(1,1)*axes(0,1));
    jacinv(1,1) = jacinv(1,1)*axes(0,0);
}

void Multiply(TPZMatrix<STATE> &A,TPZMatrix<STATE> &B,TPZMatrix<STATE> &product) {
  int rows = A.Rows(), cols = B.Cols();
  int n = A.Cols();
#ifndef NOTDEBUG
  if(n!=B.Rows()) {
    PZError << "myheader::Multiply. A.Cols and B.Rows are incompatibles.\n";
    return;
  }
#endif
  int i,j,k;
  for(i=0;i<rows;i++) {
    for(j=0;j<cols;j++) {
      product(i,j)=0.;
      for(k=0;k<n;k++)
        product(i,j) += A(i,k)*B(k,j);
    }
  }
}

/**Imprime os valores da diagonal de matrix para verificacao desde
   matrix(start,start) ateh matrix(end-1,end-1) */
int PrintDiagonal(TPZMatrix<REAL> *matrix, int mask) {
  int valpause = Pause(mask);
  if(valpause) {
    int cols = matrix->Cols();
    int start = 10, end = 50;
//    cout << setprecision(3);
    for(int r=start;r<end;r++) {
      if(cols>r) std::cout << (*matrix)(r,r);
      else std::cout << (*matrix)(r,0) << "*";
      if(r%10) std::cout << "\t";
      else  std::cout << std::endl;
    }
    std::cout << std::endl;
//    std::cout << setprecision(8);
  }
  return valpause;
}

REAL MinimumDeltaX(TPZCompMesh *mesh) {
    int64_t nel = mesh->ElementVec().NElements(), i;
    if(nel == 0) std::cout << "\nTPZCompMesh::MaximumRadiusOfMesh no elements\n";
    REAL mindist = BIGNUMBER, dist=0.0;
    for(i=0;i<nel;i++){
        TPZCompEl *com = mesh->ElementVec()[i];
        if(!com) continue;
        if(com->Type() == EInterface || com->Material()->Id() < 0) continue;
        dist = com->MaximumRadiusOfEl();
        if(dist < mindist) mindist = dist;
    }
    return mindist;
}

