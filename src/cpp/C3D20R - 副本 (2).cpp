/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "C3D20R.h"
#include "Domain.h"
#include "Generalapara.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//  Constructor
CC3D20R::CC3D20R()
{
	NEN_ = 20;	// Each element has 20 nodes
	nodes_ = new CNode*[NEN_];

	ND_ = 60;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//  Desconstructor
CC3D20R::~CC3D20R()
{
}

//	Read element data from stream Input
bool CC3D20R::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;

	Input >> N;	// element number

	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl
			<< "    Expected element : " << Ele + 1 << endl
			<< "    Provided element : " << N << endl;

		return false;
	}

	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, N14, N15, N16, N17, N18, N19, N20;	// node number

	Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> N9 >> N10 
		>> N11 >> N12 >> N13 >> N14 >> N15 >> N16 >> N17 >> N18 >> N19 >> N20 >> MSet;
	ElementMaterial_ = dynamic_cast<C3D20RMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	nodes_[4] = &NodeList[N5 - 1];
	nodes_[5] = &NodeList[N6 - 1];
	nodes_[6] = &NodeList[N7 - 1];
	nodes_[7] = &NodeList[N8 - 1];
	nodes_[8] = &NodeList[N9 - 1];
	nodes_[9] = &NodeList[N10 - 1];
	nodes_[10] = &NodeList[N11 - 1];
	nodes_[11] = &NodeList[N12 - 1];
	nodes_[12] = &NodeList[N13 - 1];
	nodes_[13] = &NodeList[N14 - 1];
	nodes_[14] = &NodeList[N15 - 1];
	nodes_[15] = &NodeList[N16 - 1];
	nodes_[16] = &NodeList[N17 - 1];
	nodes_[17] = &NodeList[N18 - 1];
	nodes_[18] = &NodeList[N19 - 1];
	nodes_[19] = &NodeList[N20 - 1];

	return true;
}

//	Write element data to stream
void CC3D20R::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele + 1 << setw(11) << nodes_[0]->NodeNumber
		<< setw(9) << nodes_[1]->NodeNumber << setw(9) << nodes_[2]->NodeNumber 
		<< setw(9) << nodes_[3]->NodeNumber << setw(9) << nodes_[4]->NodeNumber 
		<< setw(9) << nodes_[5]->NodeNumber << setw(9) << nodes_[6]->NodeNumber
		<< setw(9) << nodes_[7]->NodeNumber << setw(9) << nodes_[8]->NodeNumber
		<< setw(9) << nodes_[9]->NodeNumber << setw(9) << nodes_[10]->NodeNumber
		<< setw(9) << nodes_[11]->NodeNumber << setw(9) << nodes_[12]->NodeNumber
		<< setw(9) << nodes_[13]->NodeNumber << setw(9) << nodes_[14]->NodeNumber
		<< setw(9) << nodes_[15]->NodeNumber << setw(9) << nodes_[16]->NodeNumber
		<< setw(9) << nodes_[17]->NodeNumber << setw(9) << nodes_[18]->NodeNumber
		<< setw(9) << nodes_[19]->NodeNumber << setw(12) << ElementMaterial_->nset << endl; 
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CC3D20R::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
		for (unsigned int D = 0; D < 3; D++)
			LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 20 nodes cube element, element stiffness is a 60x60 matrix, whose upper triangular part
//	has 1830 elements
unsigned int CC3D20R::SizeOfStiffnessMatrix() { return 1830; }

//	Return the size of the element mass matrix (stored as an array column by column)
//	For 20 nodes cube element, element mass is a 60x60 matrix, whose upper triangular part
//	has 1830 elements
unsigned int CC3D20R::SizeOfMassMatrix() { return 1830; }

//	Return the size of the element damping matrix (stored as an array column by column)
//	For 20 nodes cube element, element Rayleighdamping is a 60x60 matrix, whose upper triangular part
//	has 1830 elements
unsigned int CC3D20R::SizeOfDampingMatrix() { return 1830; }


//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CC3D20R::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	int SIGN[20][3] =		// element points sign
	{
		{-1, 1, 1},
		{-1, -1, 1},
		{-1, -1, -1},
		{-1, 1, -1},
		{1, 1, 1},
		{1, -1, 1},
		{1, -1, -1},
		{1, 1, -1},
		{-1, 0, 1},
		{-1, -1, 0},
		{-1, 0, -1},
		{-1, 1, 0},
		{1, 0, 1},
		{1, -1, 0},
		{1, 0, -1},
		{1, 1, 0},
		{0, 1, 1},
		{0, -1, 1},
		{0, -1, -1},
		{0, 1, -1}
	};

	double e = 0.577350269189626;

	double GC[8][3] =
	{
		{ e, e, e },
		{ -e, e, e },
		{ -e, -e, e },
		{ e, -e, e },
		{ e, e, -e },
		{ -e, e, -e },
		{ -e, -e, -e },
		{ e, -e, -e },
	};

	double GW1 = 1;

	double GW[8] = { GW1, GW1, GW1, GW1, GW1, GW1, GW1, GW1 };

	C3D20RMaterial* material_ = dynamic_cast<C3D20RMaterial*>(ElementMaterial_);	// Pointer to material of the element
	double G = material_->E / (2 * (1 + material_->nu));
	double lambda = material_->nu*material_->E / ((1 + material_->nu)*(1 - 2 * material_->nu));

	double D[6][6] =						// Elastic matrix
	{
		{2 * G + lambda,lambda,lambda,0,0,0},
		{lambda,2 * G + lambda,lambda,0,0,0},
		{lambda,lambda,2 * G + lambda,0,0,0},
		{0,0,0,lambda,0,0},
		{0,0,0,0,lambda,0},
		{0,0,0,0,0,lambda},
	};
	
	double D_SHAPE[3][20];		//  Derivative of shape function
	double JACOBI[3][3];
	double JACOBI_DET;
	double JACOBI_INV[3][3];
	double D_SHAPE_XY[3][20];
	double B[6][60];

	double GX;
	double GY;
	double GZ;

	for (unsigned int I = 0; I < 8; I++)
	{
		clear(*D_SHAPE, 60);
		for (unsigned int J = 0; J < 8; J++)		// Caculate corner points derivative of shape function
		{
			GX = SIGN[J][0] * GC[I][0];
			GY = SIGN[J][1] * GC[I][1];
			GZ = SIGN[J][2] * GC[I][2];
			D_SHAPE[0][J] = 1. / 8 * (1 + GY) * (1 + GZ) * (2 * GX + GY + GZ - 1) * SIGN[J][0];
			D_SHAPE[1][J] = 1. / 8 * (1 + GX) * (1 + GZ) * (2 * GY + GX + GZ - 1) * SIGN[J][1];
			D_SHAPE[2][J] = 1. / 8 * (1 + GX) * (1 + GY) * (2 * GZ + GX + GY - 1) * SIGN[J][2];
		}

		for (unsigned int J = 8; J < 20; J++)			// Caculate mid points derivative of shape function
		{
			if (SIGN[J][0] == 0)
			{
				GX = GC[I][0];
				GY = SIGN[J][1] * GC[I][1];
				GZ = SIGN[J][2] * GC[I][2];
				D_SHAPE[0][J] = -1. / 2 * GX * (1 + GY) * (1 + GZ);
				D_SHAPE[1][J] = 1. / 4 * SIGN[J][1] * (1 - GX * GX)*(1 + GZ);
				D_SHAPE[2][J] = 1. / 4 * SIGN[J][2] * (1 - GX * GX)*(1 + GY);
			}
			if (SIGN[J][1] == 0)
			{
				GX = SIGN[J][0] * GC[I][0];
				GY = GC[I][1];
				GZ = SIGN[J][2] * GC[I][2];
				D_SHAPE[0][J] = 1. / 4 * SIGN[J][0] * (1 - GY * GY)*(1 + GZ);
				D_SHAPE[1][J] = -1. / 2 * GY * (1 + GX) * (1 + GZ);
				D_SHAPE[2][J] = 1. / 4 * SIGN[J][2] * (1 - GY * GY)*(1 + GX);
			}
			if (SIGN[J][2] == 0)
			{
				GX = SIGN[J][0] * GC[I][0];
				GY = SIGN[J][1] * GC[I][1];
				GZ = GC[I][2];
				D_SHAPE[0][J] = 1. / 4 * SIGN[J][0] * (1 - GZ * GZ)*(1 + GY);
				D_SHAPE[1][J] = 1. / 4 * SIGN[J][1] * (1 - GZ * GZ)*(1 + GX);
				D_SHAPE[2][J] = -1. / 2 * GZ * (1 + GX) * (1 + GY);
			}
		}

		clear(*JACOBI, 9);

		for (unsigned int J = 0; J < 20; J++)			// caculate jacobi matrix
		{
			JACOBI[0][0] += D_SHAPE[0][J] * nodes_[J]->XYZ[0];
			JACOBI[0][1] += D_SHAPE[0][J] * nodes_[J]->XYZ[1];
			JACOBI[0][2] += D_SHAPE[0][J] * nodes_[J]->XYZ[2];
			JACOBI[1][0] += D_SHAPE[1][J] * nodes_[J]->XYZ[0];
			JACOBI[1][1] += D_SHAPE[1][J] * nodes_[J]->XYZ[1];
			JACOBI[1][2] += D_SHAPE[1][J] * nodes_[J]->XYZ[2];
			JACOBI[2][0] += D_SHAPE[2][J] * nodes_[J]->XYZ[0];
			JACOBI[2][1] += D_SHAPE[2][J] * nodes_[J]->XYZ[1];
			JACOBI[2][2] += D_SHAPE[2][J] * nodes_[J]->XYZ[2];
		}

		JACOBI_DET = JACOBI[0][0] * JACOBI[1][1] * JACOBI[2][2] + JACOBI[0][1] * JACOBI[1][2] * JACOBI[2][0]	// det(JACOBI)
			+ JACOBI[0][2] * JACOBI[1][0] * JACOBI[2][1] - JACOBI[0][0] * JACOBI[1][2] * JACOBI[2][1]
			- JACOBI[0][1] * JACOBI[1][0] * JACOBI[2][2] - JACOBI[0][2] * JACOBI[1][1] * JACOBI[2][0];

		JACOBI_INV[0][0] = (JACOBI[1][1] * JACOBI[2][2] - JACOBI[1][2] * JACOBI[2][1]) / JACOBI_DET;		// inv(JACOBI)
		JACOBI_INV[0][1] = -(JACOBI[0][1] * JACOBI[2][2] - JACOBI[0][2] * JACOBI[2][1]) / JACOBI_DET;
		JACOBI_INV[0][2] = (JACOBI[0][1] * JACOBI[1][2] - JACOBI[0][2] * JACOBI[1][1]) / JACOBI_DET;
		JACOBI_INV[1][0] = -(JACOBI[1][0] * JACOBI[2][2] - JACOBI[1][2] * JACOBI[2][0]) / JACOBI_DET;
		JACOBI_INV[1][1] = (JACOBI[0][0] * JACOBI[2][2] - JACOBI[0][2] * JACOBI[2][0]) / JACOBI_DET;
		JACOBI_INV[1][2] = -(JACOBI[0][0] * JACOBI[1][2] - JACOBI[0][2] * JACOBI[1][0]) / JACOBI_DET;
		JACOBI_INV[2][0] = (JACOBI[1][0] * JACOBI[2][1] - JACOBI[1][1] * JACOBI[2][0]) / JACOBI_DET;
		JACOBI_INV[2][1] = -(JACOBI[0][0] * JACOBI[2][1] - JACOBI[0][1] * JACOBI[2][0]) / JACOBI_DET;
		JACOBI_INV[2][2] = (JACOBI[0][0] * JACOBI[1][1] - JACOBI[0][1] * JACOBI[1][0]) / JACOBI_DET;

		clear(*D_SHAPE_XY, 60);

		for (unsigned int J1 = 0; J1 < 3; J1++)			// caculate derivative of shape function in actual element
		{
			for (unsigned int J2 = 0; J2 < 20; J2++)
			{
				for (unsigned int J3 = 0; J3 < 3; J3++)
					D_SHAPE_XY[J1][J2] += JACOBI_INV[J1][J3] * D_SHAPE[J3][J2];
			}			
		}

		clear(*B, 360);

		for (unsigned int J = 0; J < 20; J++)		//  strain matrix
		{
			B[0][3 * J]     = D_SHAPE_XY[0][J];
			B[1][3 * J + 1] = D_SHAPE_XY[1][J];
			B[2][3 * J + 2] = D_SHAPE_XY[2][J];
			B[3][3 * J + 1] = D_SHAPE_XY[2][J];
			B[3][3 * J + 2] = D_SHAPE_XY[1][J];
			B[4][3 * J]     = D_SHAPE_XY[2][J];
			B[4][3 * J + 2] = D_SHAPE_XY[0][J];
			B[5][3 * J]     = D_SHAPE_XY[1][J];
			B[5][3 * J + 1] = D_SHAPE_XY[0][J];
		}

		for (unsigned int J1 = 0; J1 < 60; J1++)		// element stifness matrix
		{
			for (unsigned int J2 = 0; J2 < J1 + 1; J2++)
			{
				for (unsigned int J3 = 0; J3 < 6; J3++)
				{
					for (unsigned int J4 = 0; J4 < 6; J4++)
					{
						Matrix[J1*(J1 + 1) / 2 + J1 - J2] += B[J3][J2] * D[J3][J4] * B[J4][J1] * JACOBI_DET * GW[I];
					}
				}
			}
		}
	}
}
//	Calculate element mass matrix 
void CC3D20R::ElementMass(double* Matrix)
{
		clear(Matrix, SizeOfMassMatrix());
	int SIGN[20][3] =		// element points sign
	{
		{-1, 1, 1},
		{-1, -1, 1},
		{-1, -1, -1},
		{-1, 1, -1},
		{1, 1, 1},
		{1, -1, 1},
		{1, -1, -1},
		{1, 1, -1},
		{-1, 0, 1},
		{-1, -1, 0},
		{-1, 0, -1},
		{-1, 1, 0},
		{1, 0, 1},
		{1, -1, 0},
		{1, 0, -1},
		{1, 1, 0},
		{0, 1, 1},
		{0, -1, 1},
		{0, -1, -1},
		{0, 1, -1}
	};

	double e = 0.577350269189626;

	double GC[8][3] =
	{
		{ e, e, e },
		{ -e, e, e },
		{ -e, -e, e },
		{ e, -e, e },
		{ e, e, -e },
		{ -e, e, -e },
		{ -e, -e, -e },
		{ e, -e, -e },
	};

	double GW1 = 1;

	double GW[8] = { GW1, GW1, GW1, GW1, GW1, GW1, GW1, GW1 };

	C3D20RMaterial* material_ = dynamic_cast<C3D20RMaterial*>(ElementMaterial_);	// Pointer to material of the element
	double DENSITY = material_->den;   // The density of material

	double SHAPE[3][60] = { 0 };       // Shape function
	double D_SHAPE[3][20];	 	//  Derivative of shape function
	double JACOBI[3][3];
	double JACOBI_DET;
	double JACOBI_INV[3][3];
	double D_SHAPE_XY[3][20];
	double B[6][60];
	double Matrix1[1830] = { 0 };

	double GX;
	double GY;
	double GZ;

	for (unsigned int I = 0; I < 8; I++)
	{
		clear(Matrix1, 1830);
		clear(*SHAPE, 180);
		for (unsigned int J = 0; J < 8; J++)		// Caculate corner points shape function
		{
			GX = SIGN[J][0] * GC[I][0];
			GY = SIGN[J][1] * GC[I][1];
			GZ = SIGN[J][2] * GC[I][2];
			SHAPE[0][3 * J] = 1. / 8 * (1 + GX) * (1 + GY) * (1 + GZ) * (GX + GY + GZ - 2);
			SHAPE[1][3 * J + 1 ] = SHAPE[0][3 * J];
			SHAPE[2][3 * J + 2 ] = SHAPE[1][3 * J + 1 ];
		}
		for (unsigned int J = 8; J < 20; J++)			// Caculate mid points shape function
		{
			if (SIGN[J][0] == 0)
			{
				GX = GC[I][0];
				GY = SIGN[J][1] * GC[I][1];
				GZ = SIGN[J][2] * GC[I][2];
				SHAPE[0][3 * J] = 1. / 4 * (1 - GX * GX ) * (1 + GY) * (1 + GZ);
				SHAPE[1][3 * J + 1] = SHAPE[0][3 * J];
				SHAPE[2][3 * J + 2] = SHAPE[1][3 * J + 1];
			}
			if (SIGN[J][1] == 0)
			{
				GX = SIGN[J][0] * GC[I][0];
				GY = GC[I][1];
				GZ = SIGN[J][2] * GC[I][2];
				SHAPE[0][3 * J] = 1. / 4 * (1 + GX) * (1 - GY * GY)*(1 + GZ);
				SHAPE[1][3 * J + 1] = SHAPE[0][3 * J];
				SHAPE[2][3 * J + 2] = SHAPE[1][3 * J + 1];
			}
			if (SIGN[J][2] == 0)
			{
				GX = SIGN[J][0] * GC[I][0];
				GY = SIGN[J][1] * GC[I][1];
				GZ = GC[I][2];
				SHAPE[0][3 * J] = 1. / 4 * (1 + GX) * (1 + GY) * (1 - GZ * GZ );
				SHAPE[1][3 * J + 1] = SHAPE[0][3 * J];
				SHAPE[2][3 * J + 2] = SHAPE[1][3 * J + 1];
			}
		}
		clear(*D_SHAPE, 60);
		for (unsigned int J = 0; J < 8; J++)		// Caculate corner points derivative of shape function
		{
			GX = SIGN[J][0] * GC[I][0];
			GY = SIGN[J][1] * GC[I][1];
			GZ = SIGN[J][2] * GC[I][2];
			D_SHAPE[0][J] = 1. / 8 * (1 + GY) * (1 + GZ) * (2 * GX + GY + GZ - 1) * SIGN[J][0];
			D_SHAPE[1][J] = 1. / 8 * (1 + GX) * (1 + GZ) * (2 * GY + GX + GZ - 1) * SIGN[J][1];
			D_SHAPE[2][J] = 1. / 8 * (1 + GX) * (1 + GY) * (2 * GZ + GX + GY - 1) * SIGN[J][2];
		}

		for (unsigned int J = 8; J < 20; J++)			// Caculate mid points derivative of shape function
		{
			if (SIGN[J][0] == 0)
			{
				GX = GC[I][0];
				GY = SIGN[J][1] * GC[I][1];
				GZ = SIGN[J][2] * GC[I][2];
				D_SHAPE[0][J] = -1. / 2 * GX * (1 + GY) * (1 + GZ);
				D_SHAPE[1][J] = 1. / 4 * SIGN[J][1] * (1 - GX * GX)*(1 + GZ);
				D_SHAPE[2][J] = 1. / 4 * SIGN[J][2] * (1 - GX * GX)*(1 + GY);
			}
			if (SIGN[J][1] == 0)
			{
				GX = SIGN[J][0] * GC[I][0];
				GY = GC[I][1];
				GZ = SIGN[J][2] * GC[I][2];
				D_SHAPE[0][J] = 1. / 4 * SIGN[J][0] * (1 - GY * GY)*(1 + GZ);
				D_SHAPE[1][J] = -1. / 2 * GY * (1 + GX) * (1 + GZ);
				D_SHAPE[2][J] = 1. / 4 * SIGN[J][2] * (1 - GY * GY)*(1 + GX);
			}
			if (SIGN[J][2] == 0)
			{
				GX = SIGN[J][0] * GC[I][0];
				GY = SIGN[J][1] * GC[I][1];
				GZ = GC[I][2];
				D_SHAPE[0][J] = 1. / 4 * SIGN[J][0] * (1 - GZ * GZ)*(1 + GY);
				D_SHAPE[1][J] = 1. / 4 * SIGN[J][1] * (1 - GZ * GZ)*(1 + GX);
				D_SHAPE[2][J] = -1. / 2 * GZ * (1 + GX) * (1 + GY);
			}
		}

		clear(*JACOBI, 9);

		for (unsigned int J = 0; J < 20; J++)			// caculate jacobi matrix
		{
			JACOBI[0][0] += D_SHAPE[0][J] * nodes_[J]->XYZ[0];
			JACOBI[0][1] += D_SHAPE[0][J] * nodes_[J]->XYZ[1];
			JACOBI[0][2] += D_SHAPE[0][J] * nodes_[J]->XYZ[2];
			JACOBI[1][0] += D_SHAPE[1][J] * nodes_[J]->XYZ[0];
			JACOBI[1][1] += D_SHAPE[1][J] * nodes_[J]->XYZ[1];
			JACOBI[1][2] += D_SHAPE[1][J] * nodes_[J]->XYZ[2];
			JACOBI[2][0] += D_SHAPE[2][J] * nodes_[J]->XYZ[0];
			JACOBI[2][1] += D_SHAPE[2][J] * nodes_[J]->XYZ[1];
			JACOBI[2][2] += D_SHAPE[2][J] * nodes_[J]->XYZ[2];
		}

		JACOBI_DET = JACOBI[0][0] * JACOBI[1][1] * JACOBI[2][2] + JACOBI[0][1] * JACOBI[1][2] * JACOBI[2][0]	// det(JACOBI)
			+ JACOBI[0][2] * JACOBI[1][0] * JACOBI[2][1] - JACOBI[0][0] * JACOBI[1][2] * JACOBI[2][1]
			- JACOBI[0][1] * JACOBI[1][0] * JACOBI[2][2] - JACOBI[0][2] * JACOBI[1][1] * JACOBI[2][0];

		for (unsigned int J1 = 0; J1 < 60; J1++)		// element mass matrix
		{
			for (unsigned int J2 = 0; J2 < J1 + 1; J2++)
			{
				for (unsigned int J3 = 0; J3 < 2; J3++)
				{

						Matrix1[J1*(J1 + 1) / 2 + J1 - J2] += SHAPE[J3][J2] * SHAPE[J3][J1];

				}
				Matrix[J1*(J1 + 1) / 2 + J1 - J2] +=  DENSITY * Matrix1[J1*(J1 + 1) / 2 + J1 - J2] * JACOBI_DET*GW[I];
				//cout<<Matrix1[J1*(J1 + 1) / 2 + J1 - J2];
			}
		}
	}
}
//	Calculate element damping matrix 
void CC3D20R::ElementDamping(double* Matrix, double* Matrix1, double* Matrix2)
{
	clear(Matrix, SizeOfDampingMatrix());
	C3D20RMaterial* material_ = dynamic_cast<C3D20RMaterial*>(ElementMaterial_);	// Pointer to material of the element
	double ALPHA = material_->alpha;   // The cofficient of Rayleigh
	double BETA = material_->beta;
	for (unsigned int J1=0; J1 < SizeOfDampingMatrix(); J1++)
	{
		Matrix[J1] = Matrix1[J1] * ALPHA + Matrix2[J1] * BETA;
	}  


}

//	Calculate element efficiency stiffiness matrix 
void CC3D20R::ElementEffstiffness(double* Matrix, double* Matrix1, double* Matrix2, double* Matrix3, CGeneralapara* Generalaparas)
{
	clear(Matrix, SizeOfDampingMatrix());
	double m1,m2,m3;
	m1=Generalaparas->beta1*(1-Generalaparas->yita) * pow(Generalaparas->h,2);
	m2=Generalaparas->game*(1-Generalaparas->delta) * Generalaparas->h;
	m3=1-Generalaparas->alpha1;
	for (unsigned int J1=0; J1 < SizeOfDampingMatrix(); J1++)
	{
		Matrix[J1] = Matrix1[J1] * m1 + Matrix2[J1] * m2 + Matrix3[J1] * m3;
	//cout<<Matrix3[J1];
	}  
}

//	Calculate element stress 
void CC3D20R::ElementStress(double* stress, double* Displacement)
{
	C3D20RMaterial* material_ = dynamic_cast<C3D20RMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double G = material_->E / (2 * (1 + material_->nu));
	double lambda = material_->nu*material_->E / ((1 + material_->nu)*(1 - 2 * material_->nu));

	double D[6][6] =
	{
		{ 2 * G + lambda,lambda,lambda,0,0,0 },
		{ lambda,2 * G + lambda,lambda,0,0,0 },
		{ lambda,lambda,2 * G + lambda,0,0,0 },
		{ 0,0,0,lambda,0,0 },
		{ 0,0,0,0,lambda,0 },
		{ 0,0,0,0,0,lambda },
	};

	int SIGN[20][3] =
	{
		{ -1, 1, 1 },
		{ -1, -1, 1 },
		{ -1, -1, -1 },
		{ -1, 1, -1 },
		{ 1, 1, 1 },
		{ 1, -1, 1 },
		{ 1, -1, -1 },
		{ 1, 1, -1 },
		{ -1, 0, 1 },
		{ -1, -1, 0 },
		{ -1, 0, -1 },
		{ -1, 1, 0 },
		{ 1, 0, 1 },
		{ 1, -1, 0 },
		{ 1, 0, -1 },
		{ 1, 1, 0 },
		{ 0, 1, 1 },
		{ 0, -1, 1 },
		{ 0, -1, -1 },
		{ 0, 1, -1 }
	};

	double e = 0.577350269189626;

	double GC[8][3] =
	{
		{ -e, e, e },
		{ -e, -e, e },
		{ -e, -e, -e },
		{ -e, e, -e },
		{ e, e, e },
		{ e, -e, e },
		{ e, -e, -e },
		{ e, e, -e },
	};

	double GW1 = 1;

	double GW[8] = { GW1, GW1, GW1, GW1, GW1, GW1, GW1, GW1 };

	double S[8][6] = { 0 };
	double D_SHAPE[3][20];		//  Derivative of shape function
	double JACOBI[3][3];
	double JACOBI_DET;
	double JACOBI_INV[3][3];
	double D_SHAPE_XY[3][20];
	double B[6][60];

	double GX;
	double GY;
	double GZ;

	for (unsigned int I = 0; I < 8; I++)
	{
		clear(*D_SHAPE, 60);
		for (unsigned int J = 0; J < 8; J++)
		{
			GX = SIGN[J][0] * GC[I][0];
			GY = SIGN[J][1] * GC[I][1];
			GZ = SIGN[J][2] * GC[I][2];
			D_SHAPE[0][J] = 1. / 8 * (1 + GY) * (1 + GZ) * (2 * GX + GY + GZ - 1) * SIGN[J][0];
			D_SHAPE[1][J] = 1. / 8 * (1 + GX) * (1 + GZ) * (2 * GY + GX + GZ - 1) * SIGN[J][1];
			D_SHAPE[2][J] = 1. / 8 * (1 + GX) * (1 + GY) * (2 * GZ + GX + GY - 1) * SIGN[J][2];
		}

		for (unsigned int J = 8; J < 20; J++)
		{
			if (SIGN[J][0] == 0)
			{
				GX = GC[I][0];
				GY = SIGN[J][1] * GC[I][1];
				GZ = SIGN[J][2] * GC[I][2];
				D_SHAPE[0][J] = -1. / 2 * GX * (1 + GY) * (1 + GZ);
				D_SHAPE[1][J] = 1. / 4 * SIGN[J][1] * (1 - GX * GX)*(1 + GZ);
				D_SHAPE[2][J] = 1. / 4 * SIGN[J][2] * (1 - GX * GX)*(1 + GY);
			}
			if (SIGN[J][1] == 0)
			{
				GX = SIGN[J][0] * GC[I][0];
				GY = GC[I][1];
				GZ = SIGN[J][2] * GC[I][2];
				D_SHAPE[0][J] = 1. / 4 * SIGN[J][0] * (1 - GY * GY)*(1 + GZ);
				D_SHAPE[1][J] = -1. / 2 * GY * (1 + GX) * (1 + GZ);
				D_SHAPE[2][J] = 1. / 4 * SIGN[J][2] * (1 - GY * GY)*(1 + GX);
			}
			if (SIGN[J][2] == 0)
			{
				GX = SIGN[J][0] * GC[I][0];
				GY = SIGN[J][1] * GC[I][1];
				GZ = GC[I][2];
				D_SHAPE[0][J] = 1. / 4 * SIGN[J][0] * (1 - GZ * GZ)*(1 + GY);
				D_SHAPE[1][J] = 1. / 4 * SIGN[J][1] * (1 - GZ * GZ)*(1 + GX);
				D_SHAPE[2][J] = -1. / 2 * GZ * (1 + GX) * (1 + GY);
			}
		}

		clear(*JACOBI, 0);

		for (unsigned int J = 0; J < 20; J++)
		{
			JACOBI[0][0] += D_SHAPE[0][J] * nodes_[J]->XYZ[0];
			JACOBI[0][1] += D_SHAPE[0][J] * nodes_[J]->XYZ[1];
			JACOBI[0][2] += D_SHAPE[0][J] * nodes_[J]->XYZ[2];
			JACOBI[1][0] += D_SHAPE[1][J] * nodes_[J]->XYZ[0];
			JACOBI[1][1] += D_SHAPE[1][J] * nodes_[J]->XYZ[1];
			JACOBI[1][2] += D_SHAPE[1][J] * nodes_[J]->XYZ[2];
			JACOBI[2][0] += D_SHAPE[2][J] * nodes_[J]->XYZ[0];
			JACOBI[2][1] += D_SHAPE[2][J] * nodes_[J]->XYZ[1];
			JACOBI[2][2] += D_SHAPE[2][J] * nodes_[J]->XYZ[2];
		}

		JACOBI_DET = JACOBI[0][0] * JACOBI[1][1] * JACOBI[2][2] + JACOBI[0][1] * JACOBI[1][2] * JACOBI[2][0]
			+ JACOBI[0][2] * JACOBI[1][0] * JACOBI[2][1] - JACOBI[0][0] * JACOBI[1][2] * JACOBI[2][1]
			- JACOBI[0][1] * JACOBI[1][0] * JACOBI[2][2] - JACOBI[0][2] * JACOBI[1][1] * JACOBI[2][0];

		JACOBI_INV[0][0] = (JACOBI[1][1] * JACOBI[2][2] - JACOBI[1][2] * JACOBI[2][1]) / JACOBI_DET;
		JACOBI_INV[0][1] = -(JACOBI[0][1] * JACOBI[2][2] - JACOBI[0][2] * JACOBI[2][1]) / JACOBI_DET;
		JACOBI_INV[0][2] = (JACOBI[0][1] * JACOBI[1][2] - JACOBI[0][2] * JACOBI[1][1]) / JACOBI_DET;
		JACOBI_INV[1][0] = -(JACOBI[1][0] * JACOBI[2][2] - JACOBI[1][2] * JACOBI[2][0]) / JACOBI_DET;
		JACOBI_INV[1][1] = (JACOBI[0][0] * JACOBI[2][2] - JACOBI[0][2] * JACOBI[2][0]) / JACOBI_DET;
		JACOBI_INV[1][2] = -(JACOBI[0][0] * JACOBI[1][2] - JACOBI[0][2] * JACOBI[1][0]) / JACOBI_DET;
		JACOBI_INV[2][0] = (JACOBI[1][0] * JACOBI[2][1] - JACOBI[1][1] * JACOBI[2][0]) / JACOBI_DET;
		JACOBI_INV[2][1] = -(JACOBI[0][0] * JACOBI[2][1] - JACOBI[0][1] * JACOBI[2][0]) / JACOBI_DET;
		JACOBI_INV[2][2] = (JACOBI[0][0] * JACOBI[1][1] - JACOBI[0][1] * JACOBI[1][0]) / JACOBI_DET;

		clear(*D_SHAPE_XY, 60);

		for (unsigned int J1 = 0; J1 < 3; J1++)
		{
			for (unsigned int J2 = 0; J2 < 20; J2++)
			{
				for (unsigned int J3 = 0; J3 < 3; J3++)
					D_SHAPE_XY[J1][J2] += JACOBI_INV[J1][J3] * D_SHAPE[J3][J2];
			}
		}

		clear(*B, 360);

		for (unsigned int J = 0; J < 20; J++)
		{
			B[0][3 * J]     = D_SHAPE_XY[0][J];
			B[1][3 * J + 1] = D_SHAPE_XY[1][J];
			B[2][3 * J + 2] = D_SHAPE_XY[2][J];
			B[3][3 * J + 1] = D_SHAPE_XY[2][J];
			B[3][3 * J + 2] = D_SHAPE_XY[1][J];
			B[4][3 * J]     = D_SHAPE_XY[2][J];
			B[4][3 * J + 2] = D_SHAPE_XY[0][J];
			B[5][3 * J]     = D_SHAPE_XY[1][J];
			B[5][3 * J + 1] = D_SHAPE_XY[0][J];
		}

		for (unsigned int i = 0; i < 6; i++)
		{
			for (unsigned int j = 0; j < 60; j++)
			{
				if (LocationMatrix_[j])
				{
					for (unsigned int k = 0; k < 6; k++)
					{
						S[I][i] += D[i][k] * B[k][j] * Displacement[LocationMatrix_[j] - 1];
					}
				}
			}
		}

	}

	double a = (5 + 3 * 1.732) / 4;
	double b = -(1.732 + 1) / 4;
	double c = (1.732 - 1) / 4;
	double d = (5 - 3 * 1.732) / 4;

	double T[20][8] = 
	{
		{a,b,c,b,b,c,d,c},
		{b,a,b,c,c,b,c,d},
		{c,b,a,b,d,c,b,c},
		{b,c,b,a,c,d,c,b},
		{b,c,b,a,c,d,c,b},
		{c,b,c,d,b,a,b,c},
		{d,c,b,c,c,b,a,b},
		{c,d,c,b,b,c,b,a},
		{ (a + b) / 2,(a + b) / 2,(b + c) / 2,(b + c) / 2,(b + c) / 2,(b + c) / 2,(c + d) / 2,(c + d) / 2 },
		{ (b + c) / 2,(a + b) / 2,(a + b) / 2,(b + c) / 2,(c + d) / 2,(b + c) / 2,(b + c) / 2,(c + d) / 2 },
		{ (b + c) / 2,(b + c) / 2,(a + b) / 2,(a + b) / 2,(c + d) / 2,(c + d) / 2,(b + c) / 2,(b + c) / 2 },
		{ (a + b) / 2,(b + c) / 2,(b + c) / 2,(a + b) / 2,(b + c) / 2,(c + d) / 2,(c + d) / 2,(b + c) / 2 },
		{ (b + c) / 2,(b + c) / 2,(b + c) / 2,(a + d) / 2,(b + c) / 2,(a + d) / 2,(b + c) / 2,(b + c) / 2 },
		{ (c + d) / 2,(b + c) / 2,(b + c) / 2,(c + d) / 2,(b + c) / 2,(a + b) / 2,(a + b) / 2,(b + c) / 2 },
		{ (c + d) / 2,(c + d) / 2,(b + c) / 2,(b + c) / 2,(b + c) / 2,(b + c) / 2,(a + b) / 2,(a + b) / 2 },
		{ (b + c) / 2,(c + d) / 2,(b + c) / 2,(a + b) / 2,(b + c) / 2,(c + d) / 2,(b + c) / 2,(a + b) / 2 },
		{ (a + b) / 2,(b + c) / 2,(b + c) / 2,(a + b) / 2,(b + c) / 2,(c + d) / 2,(c + d) / 2,(b + c) / 2 },
		{ (b + c) / 2,(a + b) / 2,(b + c) / 2,(c + d) / 2,(b + c) / 2,(a + b) / 2,(b + c) / 2,(c + d) / 2 },
		{ (c + d) / 2,(b + c) / 2,(a + b) / 2,(b + c) / 2,(c + d) / 2,(b + c) / 2,(a + b) / 2,(b + c) / 2 },
		{ (b + c) / 2,(c + d) / 2,(b + c) / 2,(a + b) / 2,(b + c) / 2,(c + d) / 2,(b + c) / 2,(a + b) / 2 },
	};

	clear(stress, 120);

	for (unsigned int i = 0; i < 20; i++)
	{
		for (unsigned int j = 0; j < 6; j++)
		{
			for (unsigned int k = 0; k < 8; k++)
			{
				stress[6 * i + j] += T[i][k] * S[k][j];
			}
		}
	}


}