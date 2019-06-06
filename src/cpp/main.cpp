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
#include "Bar.h"
#include "Outputter.h"
#include "Clock.h"
#include <math.h>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 2) //  Print help message
	{
	    cout << "Usage: stap++ InputFileName\n";
		exit(1);
	}

	string filename(argv[1]);
    size_t found = filename.find_last_of('.');
	unsigned int realmode = 0;
    // If the input file name is provided with an extension
    if (found != std::string::npos) {
        if (filename.substr(found) == ".in")
            filename = filename.substr(0, found);
		else if (filename.substr(found) == ".ge")
		{
				filename = filename.substr(0, found);
				realmode = 1;
		}
        else {
            // The input file name must has an extension of 'dat'
            cout << "*** Error *** Invalid file extension: "
                 << filename.substr(found+1) << endl;
            exit(1);
        }
    }
	string InFile = filename + ".ge";
	if (realmode)
    {
		string InFile = filename + ".ge";
	}
	else
	{
		string InFile = filename + ".in";
	}
	string OutFile = filename + ".out";

	CDomain* FEMData = CDomain::Instance();

    Clock timer;
    timer.Start();
if (realmode)
{
//  Read data and define the problem domain
	if (!FEMData->ReadData(InFile, OutFile))
	{
		cerr << "*** Error *** Data input failed!" << endl;
		exit(1);
	}
    
    double time_input = timer.ElapsedTime();

//  Allocate global vectors and matrices, such as the Force, ColumnHeights,
//  DiagonalAddress and StiffnessMatrix, and calculate the column heights
//  and address of diagonal elements
	FEMData->AllocateMatricesGE();
    
//  Assemble the banded gloabl effective stiffness matrix
	FEMData->AssembleEffStiffnessMatrix();
	
	FEMData->AssembleMassMatrix();

    FEMData->AssembleStiffnessMatrix();

	FEMData->AssembleDampingMatrix();
    
    double time_assemble = timer.ElapsedTime();

//  Solve the linear equilibrium equations for displacements
	CLDLTSolver* Solver = new CLDLTSolver(FEMData->GetEffStiffnessMatrix());
    
//  Perform L*D*L(T) factorization of stiffness matrix
    Solver->LDLT();

    COutputter* Output = COutputter::Instance();

#ifdef _DEBUG_
    Output->PrintStiffnessMatrix();

    Output->PrintMassMatrix();

    Output->PrintDampingMatrix();
#endif
        
//  Loop over for all load cases
    for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
    {
		Output->OutputNsteps(ceil(FEMData->Gettotaltime()/FEMData->Geth()));
//  Loop over for all time cases
		for (unsigned int tcase = 0; tcase < ceil(FEMData->Gettotaltime()/FEMData->Geth()); tcase++)
		{

//      Assemble righ-hand-side vector (force vector)
        FEMData->AssembleForceGE(lcase + 1, tcase * FEMData->Geth());
            
//      Reduce right-hand-side force vector and back substitute
        Solver->BackSubstitutionGE(FEMData->GetForce(), FEMData->GetGeneralapara());
            
#ifdef _DEBUG_
        Output->PrintDisplacement(lcase);
#endif
            
        Output->OutputNodalDisplacement(lcase,tcase * FEMData->Geth());

//  Calculate and output stresses of all elements    
		Output->OutputElementStress();
	    }

    double time_solution = timer.ElapsedTime();


    
    double time_stress = timer.ElapsedTime();
    
    timer.Stop();
    
    *Output << "\n S O L U T I O N   T I M E   L O G   I N   S E C \n\n"
            << "     TIME FOR INPUT PHASE = " << time_input << endl
            << "     TIME FOR CALCULATION OF STIFFNESS MATRIX = " << time_assemble - time_input << endl
            << "     TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS = " << time_solution - time_assemble << endl << endl
            << "     T O T A L   S O L U T I O N   T I M E = " << time_stress << endl;

	
   }
}
return 0;
}