/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Material.h"
#include "Generalapara.h"

using namespace std;

//	Clear an array
template <class type> void clear( type* a, unsigned int N )
{
	for (unsigned int i = 0; i < N; i++)
		a[i] = 0;
}

CDomain* CDomain::_instance = nullptr;

//	Constructor
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0;

	NUMNP = 0;
	NodeList = nullptr;
	
	NUMEG = 0;
	EleGrpList = nullptr;
	
	NLCASE = 0;
	NLOAD = nullptr;
	LoadCases = nullptr;
	
	NEQ = 0;
    totaltime = 0;
	h = 0;

	Force = nullptr;
	StiffnessMatrix = nullptr;
}

//	Desconstructor
CDomain::~CDomain()
{
	delete [] NodeList;

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
	delete StiffnessMatrix;
}

//	Return pointer to the instance of the Domain class
CDomain* CDomain::Instance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile)
{
	Input.open(FileName);

	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	COutputter* Output = COutputter::Instance(OutFile);

//	Read the heading line
	Input.getline(Title, 256);
	Output->OutputHeading();

//	Read the control line
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	Read nodal point data
	if (ReadNodalPoints())
        Output->OutputNodeInfo();
    else
        return false;

//	Update equation number
	CalculateEquationNumber();
	Output->OutputEquationNumber();

//	Read general ¦Á method parameter data
    if(ReadGeneralaparameters())
	    Output->OutputGeneralparaInfo();
    else
         return false;

//	Read load data
	if (ReadLoadCases())
        Output->OutputLoadInfo();
    else
        return false;

//	Read element data
	if (ReadElements())
        Output->OutputElementInfo();
    else
        return false;

	return true;
}

//	Read nodal point data
bool CDomain::ReadNodalPoints()
{

//	Read nodal point data lines
	NodeList = new CNode[NUMNP];

//	Loop over for all nodal points
	for (unsigned int np = 0; np < NUMNP; np++)
		if (!NodeList[np].Read(Input, np))
			return false;

	return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber()
{
	NEQ = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode[dof]) 
				NodeList[np].bcode[dof] = 0;
			else
			{
				NEQ++;
				NodeList[np].bcode[dof] = NEQ;
			}
		}
	}
}

//	Read general ¦Á method parameter data
bool CDomain::ReadGeneralaparameters()
{
//	Read general ¦Á method parameter data lines
//	CGeneralapara* Generalaparas_;	// List all load cases
	Generalaparas = new CGeneralapara[9];
	if (!Generalaparas->Read(Input))
			return false;
	h = Generalaparas->h;
	totaltime = Generalaparas->totaltime;
	return true;
}


//	Read load case data
bool CDomain::ReadLoadCases()
{
//	Read load data lines
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases

//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
		if (!LoadCases[lcase].Read(Input, lcase))
			return false;

	return true;
}

// Read element data
bool CDomain::ReadElements()
{
    EleGrpList = new CElementGroup[NUMEG];

//	Loop over for all element group
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
        if (!EleGrpList[EleGrp].Read(Input))
            return false;
    
    return true;
}

//	Calculate column heights
void CDomain::CalculateColumnHeights()
{
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
        {
            CElement& Element = ElementGrp[Ele];

            // Generate location matrix
            Element.GenerateLocationMatrix();
            
			EffStiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
			
			MassMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());

			DampingMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());

			StiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
        }
    }
    
    EffStiffnessMatrix->CalculateMaximumHalfBandwidth();

    MassMatrix->CalculateMaximumHalfBandwidth();

    DampingMatrix->CalculateMaximumHalfBandwidth();

    StiffnessMatrix->CalculateMaximumHalfBandwidth();
    
#ifdef _DEBUG_
	COutputter* Output = COutputter::Instance();
	Output->PrintColumnHeights();
#endif

}

//	Assemble the banded gloabl stiffness matrix
void CDomain::AssembleStiffnessMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
            Element.ElementStiffness(Matrix);
            StiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
        }

		delete[] Matrix;
		Matrix = nullptr;
	}
}
//	Assemble the banded gloabl mass matrix
void CDomain::AssembleMassMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfMassMatrix();
		double* Matrix = new double[size];

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
            Element.ElementMass(Matrix);
            MassMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
        }

		delete[] Matrix;
		Matrix = nullptr;
	}
}
//	Assemble the banded gloabl damping matrix
void CDomain::AssembleDampingMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfDampingMatrix();
		double* Matrix = new double[size];
		double* Matrix1 = new double[size];
		double* Matrix2 = new double[size];

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
			Element.ElementStiffness(Matrix2);
			Element.ElementMass(Matrix1);
			Element.ElementDamping(Matrix,Matrix1,Matrix2);
            DampingMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
        }

		delete[] Matrix;
		Matrix = nullptr;
	}
}

//	Assemble the banded gloabl effectivestiffness matrix
void CDomain::AssembleEffStiffnessMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];
		double* Matrix1 = new double[size];
		double* Matrix2 = new double[size];
		double* Matrix3 = new double[size];
//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
			Element.ElementStiffness(Matrix2);
			Element.ElementMass(Matrix1);
			Element.ElementDamping(Matrix3,Matrix1,Matrix2);
            Element.ElementEffstiffness(Matrix,Matrix2,Matrix3,Matrix1,Generalaparas);
            EffStiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
        }

		delete[] Matrix;
		Matrix = nullptr;
	}

#ifdef _DEBUG_
	COutputter* Output = COutputter::Instance();
	Output->PrintStiffnessMatrix();
#endif

}

//	Assemble the global nodal force vector for load case LoadCase
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];

//	Loop over for all concentrated loads in load case LoadCase
	for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
        
        if(dof) // The DOF is activated
            Force[dof - 1] += LoadData->load[lnum];
	}

	return true;
}

bool CDomain::AssembleForceGE(unsigned int LoadCase, double ti)
{
	if (LoadCase > NLCASE) 
		return false;

	double Amplitudeinti = 0;//the amplitude in time at ti
	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];
    
	for (unsigned int fnum = 0; fnum < NEQ; fnum++)
	{
		Force[fnum + 2 * NEQ] = Force[ fnum + NEQ];
	}

//	Cauclate the amplitude in time at ti
	for (unsigned int tnum = 0; tnum < LoadData->ntimes; tnum++)
	{	
		if ( ( ti >= LoadData->Time[tnum] ) && ( ti <= LoadData->Time[tnum + 1] ) )
		{
			Amplitudeinti = LoadData->Amplitude[tnum] + ( ti - LoadData->Time[tnum ]) / ( LoadData->Time[tnum + 1] - LoadData->Time[tnum] ) * ( LoadData->Amplitude[tnum + 1] - LoadData->Amplitude[tnum] ) ; 
		    break;
		}
	}

	for (unsigned int i = 0; i < NEQ; i++)
	{
		Force[NEQ + i] = 0;
	}
//	Loop over for all concentrated loads in load case LoadCase
	for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
        
        if(dof) // The DOF is activated
            Force[NEQ + dof - 1] += Amplitudeinti * LoadData->load[lnum];
	}
//	CGeneralapara* generalapara_ = Generalaparas_;	// Pointer to parameter of the general ¦Á method
	double m1,m2,m3,m4,m5,m6,m7,m8,m9,ck;
	m1=Generalaparas->beta1 * (1 - Generalaparas->yita) * pow(Generalaparas->h,2);
	m2=Generalaparas->game * (1 - Generalaparas->delta) * Generalaparas->h;
	m3=1 - Generalaparas->alpha1;
	m4=Generalaparas->beta1 * Generalaparas->yita * pow(Generalaparas->h,2);
	m5=Generalaparas->h*m3;
	m6=(Generalaparas->yipusi * m3 - Generalaparas->alpha1 * Generalaparas->beta1) * pow(Generalaparas->h,2);
	m7=(Generalaparas->game * (1 - Generalaparas->delta) - Generalaparas->beta1) * pow(Generalaparas->h,2);
	m8=(1-Generalaparas->delta) * (Generalaparas->yipusi * Generalaparas->game - Generalaparas->beta1 * Generalaparas->mu) * pow(Generalaparas->h,3);
	m9=-Generalaparas->beta1 * Generalaparas->yita * pow(Generalaparas->h,2);
//**************************************
/*	m1=1;
	m4=0;*/
//**************************************
	for (unsigned int fnum = 0; fnum < NEQ; fnum++)
	{
		Force[fnum] =m1 * Force[ fnum + NEQ] + m4 * Force[fnum + 2 * NEQ] ;
	}
	CSkylineMatrix<double>* K;
	K = GetMassMatrix();
	unsigned int N = K->dim();
    unsigned int* ColumnHeights = K->GetColumnHeights();   // Column Hights
	double* Force1;
	Force1 = new double[NEQ];
/*	for (unsigned int i = 0; i < NEQ; i++)
	{
		Force1[i] = 0;
	}*/
	memset(Force1,0,sizeof(double)*NEQ);
//*****************************************
/*	m3=(1-Generalaparas->alpha1)/Generalaparas->beta1/pow(Generalaparas->h,2);
	m5=Generalaparas->h*m3;
	m6=m5*Generalaparas->h/2-1;*/
//*****************************************
	//cout<<Force1[1];
/*	for (unsigned int fnum = 0; fnum < NEQ; fnum++)
	{*/
	    for (unsigned int i = 0; i < N; i++)
	   {
		   for (unsigned int j = 0; j < N; j++)
		   {
			   if ((i > j) && ((i - j) <= ColumnHeights[i]))
			   {
			   Force1[i] +=  (*K)(j+1,i+1) * (m3 * Force[ j + NEQ * 3 ] + m5 * Force[ j + NEQ * 4] + m6 * Force [ j + NEQ * 5 ]);
			   //cout<<(*K)(j+1,i+1)<<Force1[i];
			   }
			   else if((i <= j) && ((j - i) <= ColumnHeights[j]))
			   {
			   Force1[i] +=  (*K)(i+1,j+1) * (m3 * Force[ j + NEQ * 3 ] + m5 * Force[ j + NEQ * 4] + m6 * Force [ j + NEQ * 5 ]);
			   //cout<<(*K)(i+1,j+1)<<Force1[i];
			   }
		   }
	   }
//	}
	K = GetDampingMatrix();
	double* Force2;
	Force2 = new double[NEQ];
    ck=1-Generalaparas->delta;
//*************************
/*	m2=ck*Generalaparas->game/Generalaparas->beta1/Generalaparas->h;
	m7=ck*Generalaparas->game/Generalaparas->beta1-1;
	m8=ck*(Generalaparas->game/2/Generalaparas->beta1-1)*Generalaparas->h;*/
//*************************
/*	for (unsigned int i = 0; i < NEQ; i++)
	{
		Force2[i] = 0;
	}*/
	memset(Force2,0,sizeof(double)*NEQ);
/*	for (unsigned int fnum = 0; fnum < NEQ; fnum++)
	{*/
	    for (unsigned int i = 0; i < N; i++)
	   {
		   for (unsigned int j = 0; j < N; j++)
		   {
			   if ((i > j) && ((i - j) <= ColumnHeights[i]))
			   {
			   Force2[i] +=  (*K)(j+1,i+1) * (m2 * Force[ j + NEQ * 3 ] + m7 * Force[ j + NEQ * 4] + m8 * Force [ j + NEQ * 5 ]);
			   //cout<<(*K)(j+1,i+1);
			   }
			   else if((i <= j) && ((j - i) <= ColumnHeights[j]))
			   {
			   Force2[i] +=  (*K)(i+1,j+1) * (m2 * Force[ j + NEQ * 3 ] + m7 * Force[ j + NEQ * 4] + m8 * Force [ j + NEQ * 5 ]);
			   //cout<<(*K)(j+1,i+1);
			   }
		   }
	   }
//	}
	K = GetStiffnessMatrix();
	double* Force3;
	Force3 = new double[NEQ];
/*	for (unsigned int i = 0; i < NEQ; i++)
	{
		Force3[i] = 0;
	}*/
    memset(Force3,0,sizeof(double)*NEQ);
//********************************
/*	m9=-Generalaparas->delta;*/
//********************************
/*	for (unsigned int fnum = 0; fnum < NEQ; fnum++)
	{*/
	    for (unsigned int i = 0; i < N; i++)
	   {
		   for (unsigned int j = 0; j < N; j++)
		   {
			   if ((i > j) && ((i - j) <= ColumnHeights[i]))
			   {
			   Force3[i] +=  (*K)(j+1,i+1) * (m9 * Force[ j + NEQ * 3 ]);
			//   cout<<(*K)(j+1,i+1);
			   }
			   else if((i <= j) && ((j - i) <= ColumnHeights[j]))
			   {
			   Force3[i] +=  (*K)(i+1,j+1) * (m9 * Force[ j + NEQ * 3 ]);
			 //  cout<<(*K)(i+1,j+1);
			   }
		   }
	   }
//	}
	for (unsigned int fnum = 0; fnum < NEQ; fnum++)
	{
		Force[fnum] +=(Force1[fnum] + Force2[fnum] + Force3[fnum]) ;
		//cout<<Force[fnum];
		//cout<<Force1[fnum];
	}

	return true;
}


//	Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//	and calculate the column heights and address of diagonal elements
void CDomain::AllocateMatrices()
{
//	Allocate for global force/displacement vector
	Force = new double[NEQ];
    clear(Force, NEQ);

//  Create the banded stiffness matrix
    StiffnessMatrix = new CSkylineMatrix<double>(NEQ);

//	Calculate column heights
	CalculateColumnHeights();

//	Calculate address of diagonal elements in banded matrix
	StiffnessMatrix->CalculateDiagnoalAddress();

//	Allocate for banded global stiffness matrix
    StiffnessMatrix->Allocate();

	COutputter* Output = COutputter::Instance();
	Output->OutputTotalSystemData();
}

void CDomain::AllocateMatricesGE()
{
//	Allocate for global force/displacement vector
	unsigned int NEQGE;
	NEQGE = 6 * NEQ;
	Force = new double[NEQGE];
    clear(Force, NEQGE);
/*	for (unsigned int i = 0; i < NEQGE; i++)
	{
		Force[i] = 0;
	}*/
	memset(Force,0,sizeof(double)*NEQGE);

//  Create the banded stiffness matrix
    StiffnessMatrix = new CSkylineMatrix<double>(NEQ);

    MassMatrix = new CSkylineMatrix<double>(NEQ);

    DampingMatrix = new CSkylineMatrix<double>(NEQ);

	EffStiffnessMatrix = new CSkylineMatrix<double>(NEQ);
//	Calculate column heights
	CalculateColumnHeights();

	MassMatrix->CalculateDiagnoalAddress();

	MassMatrix->Allocate();
//	Calculate address of diagonal elements in banded matrix
	StiffnessMatrix->CalculateDiagnoalAddress();

//	Allocate for banded global stiffness matrix
    StiffnessMatrix->Allocate();


	DampingMatrix->CalculateDiagnoalAddress();

	DampingMatrix->Allocate();

	EffStiffnessMatrix->CalculateDiagnoalAddress();

	EffStiffnessMatrix->Allocate();

	COutputter* Output = COutputter::Instance();
	Output->OutputTotalSystemData();
}