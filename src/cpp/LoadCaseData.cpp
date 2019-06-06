/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "LoadCaseData.h"

#include <iomanip>
#include <iostream>

using namespace std;

CLoadCaseData :: ~CLoadCaseData()
{
	delete [] node;
	delete [] dof;
	delete [] load;
}

void CLoadCaseData :: Allocate(unsigned int num, unsigned int num2)
{
	nloads = num;
	node = new unsigned int[nloads];
	dof = new unsigned int[nloads];
	load = new double[nloads];
	ntimes = num2;
	Time = new double[ntimes];
	Amplitude = new double[ntimes];
}; 

//	Read load case data from stream Input
bool CLoadCaseData :: Read(ifstream& Input, unsigned int lcase)
{
//	Load case number (LL) and number of concentrated loads in this load case(NL) and number of Tables(NT)
	
	unsigned int LL, NL, NT;

	Input >> LL >> NL >> NT;	

	if (LL != lcase + 1) 
	{
		cerr << "*** Error *** Load case must be inputted in order !" << endl 
			 << "   Expected load case : " << lcase + 1 << endl
			 << "   Provided load case : " << LL << endl;

		return false;
	}

	Allocate(NL ,NT);

	for (unsigned int i = 0; i < NL; i++)
	{	
		Input >> node[i] >> dof[i] >> load[i];
	}
	for (unsigned int j = 0; j < NT; j++)
	{
		Input >> Time[j] >> Amplitude[j];
	}
	return true;
}

//	Write load case data to stream
void CLoadCaseData::Write(COutputter& output, unsigned int lcase)
{
	for (unsigned int i = 0; i < nloads; i++)
	{
		output << setw(7) << node[i] << setw(13) << dof[i]  << setw(19) << load[i] << endl;
	}
	for (unsigned int j = 0; j < ntimes; j++)
	{
		output << setw(7) << Time[j] << setw(13) << Amplitude[j]  << setw(19) << endl;
	}
}
