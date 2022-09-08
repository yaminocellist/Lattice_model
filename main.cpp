#include "settings.h"
#include "functions.h"

/*******************************************************************
 *   1st argument: (mg) binding site length of 1st type of protein;
 *   2nd argument: (Vg) max interaction distance of 1st type of protein;
 *   3rd argument: (K_n_g) binding constant for 1st type of protein;
 *   4th argument: (C_0_g) concentration of 1st type of protein;
 *   ... ...
 *   and so on:
********************************************************************/



int main(int argc, char *argv[]) {
    Parameters p;
    set_param(p, argc, argv);

    /********************************************
    *   Claiming input variables:
    ********************************************/
    std::cout << "There are " << p.f << " types of proteins involved in competition:" << std::endl << std::endl;
    for (int i = 0; i< p.f; i++) {
        std::cout << "Type " << i+1 << ", g = " << i << ", mg[" << i+1 << "] = " << p.mg[i] << ", Vg[" << i+1 << "] = " << p.Vg[i] << ";" << std::endl;
    }
    
    std::cout << "The results:" << std::endl;
    // std::cout << p.sum_mgVg(0) << "  " << p.sum_mgVg(1) << "  " << p.sum_mgVg(2) << "  " << p.sum_mgVg(3) << std::endl;
    // std::cout << "It should be:" << std::endl << "20  26  39  40" << std::endl; 
    std::cout << p.sum_mgVg(2)+1 << std::endl;
}