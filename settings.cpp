#include "settings.h"
#include <iostream>
#include <eigen3/Eigen/Dense>

double pi = boost::math::constants::pi<double>();

void set_param (Parameters &p, int argc, char *argv[]) {
    // number of type of DNA:
    p.f = floor(argc/4);                 
    p.N = 147;              // 146 or 147 for nucleosome;
    // mg of each type g:                  
    for (int i = 0; i < argc - 1; i++) {
        if (((i + 1) % 4) == 1) {
            p.mg.push_back(std::stoi(argv[i + 1]));
        }
    }

    // Vg of each type g:
    for (int i = 0; i < argc - 1; i++) {
        if (((i + 1) % 4) == 2) {
            p.Vg.push_back(std::stoi(argv[i + 1]));
        }
    }

    /*  binding constant of each type g,
        the binding sites is set to start from n = 3:  */
    vector<double> temp;    // temporally stores input K values;
    for (int i = 0; i < argc - 1; i++) {
        if (((i + 1) % 4) == 3) {
            temp.push_back(std::stod(argv[i + 1]));
        }
    }
    p.K_n_g.resize(p.N);
    for (int i = 0; i < p.N; i++) {
        p.K_n_g[i].resize(p.f);
    }
    for (int n = 0; n < p.N; n++) {
        for (int g = 0; g < p.f; g++) {
            if (n >= 3 && n < 3 + p.mg[g]) {
                p.K_n_g[n][g] = std::pow(temp[g], 1./p.mg[g]);
            }
            else {
                p.K_n_g[n][g] = 0;
            }
        }
    }
    std::vector<double>().swap(temp);   // clear up temporal vector variable;

    // C_g of each type g:
    for (int i = 0; i < argc - 1; i++) {
        if (((i + 1) % 4) == 0) {
            p.C_0_g.push_back(std::stod(argv[i + 1]));
        }
    }

    /*  Unwrap of each type g:  */
    for (int i = 0; i < p.f; i++) {
        p.Unwrap.push_back(1.);
    }

    p.w_g1_g2.resize(p.f);
    for (int i = 0; i < p.f; i++) {
        p.w_g1_g2[i].resize(p.f, 1);
    }

    // Calculate the size of Q_n:
    for (int i = 0; i < argc - 1; i++) {
        if (((i + 1) % 4) == 1 || ((i + 1) % 4) == 2)
            p.dimensions += std::stoi(argv[i + 1]);
    }
    auto it = std::max_element(p.Vg.begin(), p.Vg.end());
    p.dimensions += p.Vg[it - p.Vg.begin()] + 3;
    /*  Resize the Q_n matrix:  */
    p.Q_n.resize(p.N);
    for (int i = 0; i < p.N; i++) {
        p.Q_n[i].resize(p.dimensions);
        for (int j = 0; j < p.dimensions; j++) {
            p.Q_n[i][j].resize(p.dimensions, 0.);
        }
    }
    
    /*  Set value for each elements of Q_n:  */
    // set_Q(p);
}

/***********************************************************
 *  index 'n' indicates the nth segment of 1-d DNA molecule;
 *  index 'i' and 'j' indicates one of the STATEs,
 *      'i' indicates state of current segment,
 *      'j' indicates state of its next segment,
***********************************************************/
/*****************************************************************************************
 *  Look-up table for non-zero Q_n values: 
 *  value 0 means impossible or not considered;
 *  value 1 means energy is zero;
******************************************************************************************/
void set_Q (Parameters &p) {
    for (int n = 0; n < p.Q_n.size(); n++) {
        for (int i = 0; i < p.Q_n[n].size(); i++) {
            for (int j = 0; j < p.Q_n[n][i].size(); j++) {
                /*  Table 1.:  */
                for (int g = 0; g < p.f; g++) {
                    if (i == sum_mg(p.mg, g)) {
                        if (n == 0) {
                            p.Q_n[n][i][i + 1] = p.K_n_g[n][g] * p.C_0_g[g];
                        }
                        else {
                            p.Q_n[n][i][i + 1] = p.K_n_g[n][g];
                        }
                    }
                }

                /*  Table 2.:  */
                for (int g = 0; g < p.f; g++) {
                    if (i > sum_mg(p.mg, g) && i < sum_mg(p.mg, g + 1) - 1) {
                        p.Q_n[n][i][i + 1] = p.K_n_g[n][g];
                    }
                }
            }
        }
    }
}
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
    std::cout << "The results:" << std::endl;

    std::cout << p.f << std::endl;
    std::cout << sum_mg(p.mg, 0) << ' ' << sum_mg(p.mg, 1) << std::endl;
    std::cout << sum_mg(p.mg, 1) << ' ' << sum_mg(p.mg, 2) << std::endl;
    std::cout << sum_mg(p.mg, 2) << ' ' << sum_mg(p.mg, 3) << std::endl;
    // std::cout << p.Q_n[144][99][95] << std::endl;
    // std::cout << first_non_zero(p.K_n_g[7]) << std::endl;
}
