#include "settings.h"
#include <iostream>

double pi = boost::math::constants::pi<double>();

void set_param (Parameters &p, int argc, char *argv[]) {
    int d = floor(argc/2);
    p.mg.resize(d);       // decide the size;
    p.Vg.resize(d);
    p.K_n_g.resize(d);
    p.C_0_g.resize(d);
    p.N = std::floor(p.DNA_length/p.b_B);
    for (int i = 0; i < argc - 1; i++) {
        p.dimensions += std::atoi(argv[i+1]);
    }

    /*  Resize the Q_n matrix:  */
    p.Q_n.resize(p.dimensions);
    for (int i = 0; i < p.dimensions; i++) {
        p.Q_n[i].resize(p.dimensions);
    }

}

/*******************************************************************
 *   1st argument: binding site length of 1st type of protein;
 *   2nd argument: max interaction distance of 1st type of protein;
 *   ... ...
 *   and so on:
********************************************************************/
int main(int argc, char *argv[]) {
    Parameters p;
    
    set_param(p, argc, argv);
    std::cout << p.dimensions << std::endl;
}