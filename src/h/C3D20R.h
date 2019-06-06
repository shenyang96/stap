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

#include "Element.h"

using namespace std;

// C3D20R element class
class CC3D20R : public CElement
{
public:

// Constructor
	CC3D20R();

// Desconstructor
	~CC3D20R();

// Read element data from stream Input
	virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);

//	Write element data to stream
	virtual void Write(COutputter& output, unsigned int Ele);

// Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1
	virtual void GenerateLocationMatrix();

//	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//	Calculate element mass matrix
	virtual void ElementMass(double* Matrix);

//	Calculate element damping matrix
	virtual void ElementDamping(double* Matrix, double* Matrix1, double* Matrix2);

//	Calculate element effectivestiffness matrix
	virtual void ElementEffstiffness(double* Matrix, double* Matrix1, double* Matrix2, double* Matrix3, CGeneralapara* Generalaparas);

//	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

//	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix();

//	Return the size of the element mass matrix (stored as an array column by column)
	virtual unsigned int SizeOfMassMatrix();

//	Return the size of the element damping matrix (stored as an array column by column)
	virtual unsigned int SizeOfDampingMatrix();

};
