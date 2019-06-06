/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Outputter.h"

using namespace std;

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CGeneralapara
{
public:

	double alpha1;  //¦Á
	double delta;  //¦Ä
	double yita;  //¦Ç
	double yipusi;  //¦Å
	double beta1;  //¦Â
	double mu;  //¦Ì
	double game;  //¦Ã
	double h;  // Solution interval
	double totaltime;//the total time of solution

public:

//! Virtual deconstructor
    virtual ~CGeneralapara() {};

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
    virtual void Write(COutputter& output);

};