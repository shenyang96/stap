/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Generalapara.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CGeneralapara::Read(ifstream& Input)
{

	Input >> alpha1 >> delta >> yita >> yipusi >> beta1 >> mu >> game >> h >> totaltime;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream
void CGeneralapara::Write(COutputter& output)
{
	output << setw(5) << setw(16) << alpha1 << setw(16) << delta << setw(16) << yita << setw(16) << yipusi << setw(16) << beta1 << setw(16) << mu << setw(16) << game << setw(16) << h <<setw(16) << totaltime << endl;
}