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
    int target_n = std::stoi(argv[argc - 2]);
    int target_g = std::stoi(argv[argc - 1]);
    /************************************************************************************************************
    *                                       Claiming input variables:
    *************************************************************************************************************/
    std::cout << "There are " << p.f << " types of proteins involved in competition:" << std::endl;
    for (int i = 0; i< p.f; i++) {
        std::cout << "Type " << i+1 << ", g = " << i << ", mg[" << i+1 << "] = " << p.mg[i] << ", Vg[" << i+1 << "] = " << p.Vg[i] << ";" << std::endl;
    }

    std::cout << std::endl;
    std::cout << "------------    Results    ------------" << std::endl << std::endl;

    MatrixXd M = p.Q_n[0];
    for (int i = 1; i < p.Q_n.size(); i++) {
        M = M*p.Q_n[i];
    }
    
    // std::cout << "The Qn matrix is: " << std::endl << M << std::endl;
    /***********************************************************
    *   Summing all elements of Qn to get partition function Z:
    ***********************************************************/
    double Z = 0.;
    for (int i = 0; i < M.rows(); i++) {
        for (int j = 0; j < M.cols(); j++) {
            Z += M(i, j); 
        }
    }
    std::cout << "The original Z is: " << Z << std::endl;

    /****************************************************
    *                   To get ∂Z:
    ****************************************************/
    vector<MatrixXd> Qn_delta = set_Qn_delta(p, target_n, target_g);  // newly set Qn matrix using K_n_g*(1 + delta);
    //  Calculating the multiplication over all Qn matrix through each DNA segments, from n ~ [0, p.N);
    MatrixXd Qn0 = Qn_delta[0];
    for (int i = 1; i < Qn_delta.size(); i++) {
        Qn0 *= Qn_delta[i];
    }

    double Z_delta = 0.;
    for (int i = 0; i < Qn0.rows(); i++) {
        for (int j = 0; j < Qn0.cols(); j++) {
            Z_delta += Qn0(i, j); 
        }
    }

    double delta_z = Z_delta - Z;
    double delta_Kng = p.K_n_g[target_n][target_g] * p.delta;
    std::cout << "The ∂Z is:         " << delta_z << std::endl;
    std::cout << "The ∂K(" << target_n << ", " << target_g << ") is:   " << delta_Kng << std::endl;
    double c_n_g = (delta_z * p.K_n_g[target_n][target_g])/(delta_Kng * Z);
    std::cout << "The probability for binding type " << target_g + 1 << " protein on the " << target_n+1 << "th segment is: " << std::endl;

    std::cout << "\033[1;31m" << "                " << c_n_g*100. << "%" << "\033[0m" << std::endl;
    double runtime = global.elapsed();
    
    std::cout << "Number of DNA segments is " << p.N << ", total elapsed time is " << runtime << " seconds.\n";
}
