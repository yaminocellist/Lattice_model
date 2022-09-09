#include "settings.h"
#include "functions.h"
#include "main.h"
#include <eigen3/Eigen/StdVector>

/*******************************************************************
 *   1st argument: (mg) binding site length of 1st type of protein;
 *   2nd argument: (Vg) max interaction distance of 1st type of protein;
 *   3rd argument: (K_n_g) binding constant for 1st type of protein;
 *   4th argument: (C_0_g) concentration of 1st type of protein;
 *   ... ...
 *   and so on:
********************************************************************/
int main(int argc, char *argv[]) {
    Timer global;
    Parameters p;
    set_param(p, argc, argv);

    std::cout << "The results:" << std::endl;
    
    MatrixXd M = p.Q_n[0];
    for (int i = 1; i < p.Q_n.size(); i++) {
        M = M*p.Q_n[i];
    }
    
    std::cout << "The Qn matrix is: " << std::endl << M << std::endl;
    /***********************************************************
    *   Summing all elements of Qn to get partition function Z:
    ***********************************************************/
    double Z = 0.;
    for (int i = 0; i < M.rows(); i++) {
        for (int j = 0; j < M.cols(); j++) {
            Z += M(i, j); 
        }
    }
    std::cout << "The partition function result is: " << Z << std::endl;

    /********************************************
    *   Claiming input variables:
    ********************************************/
    std::cout << "There are " << p.f << " types of proteins involved in competition:" << std::endl << std::endl;
    for (int i = 0; i< p.f; i++) {
        std::cout << "Type " << i+1 << ", g = " << i << ", mg[" << i+1 << "] = " << p.mg[i] << ", Vg[" << i+1 << "] = " << p.Vg[i] << ";" << std::endl;
    }
    double runtime = global.elapsed();
    
    std::cout << "Number of DNA segments is " << p.N << ", total elapsed time is " << runtime << " seconds.\n";
}
