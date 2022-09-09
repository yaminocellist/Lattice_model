#include "settings.h"
#include <iostream>


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
        p.Unwrap.push_back(.5);
    }

    p.w_g1_g2.resize(p.f);
    for (int i = 0; i < p.f; i++) {
        p.w_g1_g2[i].resize(p.f, 1);
    }

    /*  Upper limit of bound states:  */
    p.bound_cap = p.sum_mg(p.f); 

    /*  Upper limit of g1-g2 gap states:  */
    p.gap_cap = p.sum_mgVg(p.f);

    // Calculate the size of Q_n:
    for (int i = 0; i < argc - 1; i++) {
        if (((i + 1) % 4) == 1 || ((i + 1) % 4) == 2)
            p.dimensions += std::stoi(argv[i + 1]);
    }
    auto it = std::max_element(p.Vg.begin(), p.Vg.end());
    p.dimensions += p.Vg[it - p.Vg.begin()] + 3;
    /*  Resize the Q_n matrix:  */
    for (int i = 0; i < p.N; i++) {
        p.Q_n.push_back(MatrixXd::Zero(p.dimensions, p.dimensions));
    }

    // p.Q_n.resize(p.N);
    // for (int i = 0; i < p.N; i++) {
    //     p.Q_n[i].resize(p.dimensions);
    //     for (int j = 0; j < p.dimensions; j++) {
    //         p.Q_n[i][j].resize(p.dimensions, 0.);
    //     }
    // }
    
    /*  Set value for each elements of Q_n:  */
    set_Q_n(p);
}


// void set_Q (Parameters &p) {
//     for (int n = 0; n < p.Q_n.size(); n++) {
//         for (int i = 0; i < p.Q_n[n].size(); i++) {
//             for (int j = 0; j < p.Q_n[n][i].size(); j++) {
//                 /*  Table 1.:  */
//                 for (int g = 0; g < p.f; g++) {
//                     if (i == p.sum_mg(g)) {
//                         if (n == 0) {
//                             p.Q_n[n][i][i + 1] = p.K_n_g[n][g] * p.C_0_g[g];
//                         }
//                         else {
//                             p.Q_n[n][i][i + 1] = p.K_n_g[n][g];
//                         }
//                     }
//                 }

//                 /*  Table 2.:  */
//                 for (int g = 0; g < p.f; g++) {
//                     if (i > p.sum_mg(g) && i < p.sum_mg(g + 1) - 1) {
//                         p.Q_n[n][i][i + 1] = p.K_n_g[n][g];
//                     }
//                 }

//                 /*  Table 3.:  */
//                 if (i < p.bound_cap) {
//                     int temp_g = p.find_g(i);
//                     if (n == 0) {
//                         p.Q_n[n][i][p.bound_cap + 1] = p.K_n_g[n][temp_g] * p.C_0_g[temp_g] * p.Unwrap[temp_g];
//                     }
//                     else {
//                         p.Q_n[n][i][p.bound_cap + 1] = p.K_n_g[n][temp_g] * p.Unwrap[temp_g];
//                     }
//                 }

//                 /*  Table 4.:  */
//                 if (i < p.bound_cap) {
//                     int temp_g1 = p.find_g(i);
//                     int temp_g2 = p.find_g(j);
//                     if (j == p.sum_mg(temp_g2)) {
//                         if (n == 0) {
//                             p.Q_n[n][i][j] = p.K_n_g[n][temp_g1]*p.w_g1_g2[temp_g1][temp_g2]*p.Unwrap[temp_g1]*p.C_0_g[temp_g1]*p.C_0_g[temp_g2];
//                         }
//                         else {
//                             p.Q_n[n][i][j] = p.K_n_g[n][temp_g1]*p.w_g1_g2[temp_g1][temp_g2]*p.Unwrap[temp_g1]*p.C_0_g[temp_g1];
//                         }
//                     }
//                 }

//                 /*  Table 5.:  */
//                 if (i < p.bound_cap && j < p.bound_cap) {
//                     int temp_g1 = p.find_g(i);
//                     int temp_g2 = p.find_g(j);
//                     if (i == p.sum_mg(temp_g1 + 1) - 1 && (j > p.sum_mg(temp_g2) && j < p.sum_mg(temp_g2 + 1))) {
//                         p.Q_n[n][i][j] = p.K_n_g[n][temp_g1]*p.w_g1_g2[temp_g1][temp_g2]*p.Unwrap[temp_g2]*p.C_0_g[temp_g2];
//                     }
//                 }

//                 /*  Table 6.:  */
//                 if (i == p.bound_cap && j == p.bound_cap) {
//                     p.Q_n[n][i][j] = 1.;
//                 }

//                 /*  Table 7.:  */
//                 if (i == p.bound_cap + 1 && j == p.bound_cap + 1) {
//                     p.Q_n[n][i][j] = 1.;
//                 }

//                 /*  Table 8.:  */
//                 if (i == p.bound_cap && j < p.bound_cap) {
//                     int g2 = p.find_g(j);
//                     p.Q_n[n][i][j] = p.w_g1_g2[0][g2]*p.Unwrap[g2]*p.C_0_g[g2];
//                 }

//                 /*  Table 9.:  */
//                 if (i < p.bound_cap) {
//                     int g1 = p.find_g(i);
//                     if (j == p.gap_cap + 1) {
//                         if (n == 1) {
//                             p.Q_n[n][i][j] = p.w_g1_g2[g1][0]*p.K_n_g[0][g1]*p.Unwrap[g1]*p.C_0_g[g1];
//                         }
//                         else {
//                             p.Q_n[n][i][j] = p.w_g1_g2[g1][0]*p.K_n_g[0][g1];
//                         }
//                     }
//                 }

//                 /*  Table 10.:  */
//                 if (i > p.gap_cap && i <= p.dimensions - 2) {
//                     p.Q_n[n][i][i + 1] = 1.;
//                 }

//                 /*  Table 11.:  */
//                 if (i == p.dimensions - 1) {
//                     p.Q_n[n][i][i] = 1.;
//                 }

//                 /*  Table 12.:  */
//                 if (i == p.dimensions - 1 && j < p.bound_cap) {
//                     int g2 = p.find_g(j);
//                     p.Q_n[n][i][j] = p.Unwrap[g2]*p.C_0_g[g2];
//                 }

//                 /*  Table 13.:  */
//                 if (i < p.bound_cap && (j > p.bound_cap + 1 && j <= p.gap_cap)) {
//                     int g1 = p.find_g(i);
//                     int g2 = p.find_g(j);
//                     if (n == 0) {
//                         p.Q_n[n][i][j] = p.w_g1_g2[g1][g2]*p.K_n_g[n][g1]*p.Unwrap[g1]*p.C_0_g[g1];
//                     }
//                     else {
//                         p.Q_n[n][i][j] = p.w_g1_g2[g1][g2]*p.K_n_g[n][g1]*p.Unwrap[g1];
//                     }
//                 }

//                 /*  Table 14.:  */
//                 if (i > p.bound_cap + 2 && i <= p.gap_cap) {
//                     p.Q_n[n][i][i - 1] = 1;
//                 }

//                 /*  Table 15.:  */
//                 if ((i > p.bound_cap + 1 && i <= p.gap_cap) && j < p.bound_cap) {
//                     int g1 = p.find_g(i);
//                     int g2 = p.find_g(j);
//                     if (i == p.sum_mgVg(g1) + 1) {
//                         p.Q_n[n][i][j] = p.Unwrap[g2]*p.C_0_g[g2];
//                     }
//                 }
//             }
//         }
//     }
// }

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
void set_Q_n (Parameters &p) {
    for (int n = 0; n < p.Q_n.size(); n++) {
        /*  Algorithm 6:  */
        p.Q_n[n](p.bound_cap, p.bound_cap) = 1.;

        /*  Algorithm 7:  */
        p.Q_n[n](p.bound_cap + 1, p.bound_cap + 1) = 1.;
            
        /*  Algorithm 11:  */
        p.Q_n[n](p.dimensions - 1, p.dimensions - 1) = 1.;

        for (int i = 0; i < p.Q_n[n].rows(); i++) {
            int g1 = p.find_g(i);
            /*  Algorithm 1:  */
            if (i == p.sum_mg(g1)) {
                if (n == 3) {
                    p.Q_n[n](i, i+1) = p.K_n_g[n][g1]*p.C_0_g[g1];
                }
                else {
                    p.Q_n[n](i, i+1) = p.K_n_g[n][g1];
                }
            }

            /*  Algorithm 2:  */
            if (i < p.sum_mg(g1 + 1) - 1 && i > p.sum_mg(g1)) {
                p.Q_n[n](i, i+1) = p.K_n_g[n][g1];
            }

            /*  Algorithm 3:  */
            if (i < p.bound_cap) {
                if (n == 3) {
                    p.Q_n[n](i, p.bound_cap + 1) = p.K_n_g[n][g1]*p.Unwrap[g1]*p.C_0_g[g1];
                }
                else {
                    p.Q_n[n](i, p.bound_cap + 1) = p.K_n_g[n][g1]*p.Unwrap[g1];
                }
            }

            /*  Algorithm 9:  */
            if (i < p.bound_cap) {
                if (n == 3) {
                    p.Q_n[n](i, p.gap_cap + 1) = p.w_g1_g2[g1][0]*p.K_n_g[n][g1]*p.Unwrap[g1]*p.C_0_g[g1];
                }
                else {
                    p.Q_n[n](i, p.gap_cap + 1) = p.w_g1_g2[g1][0]*p.K_n_g[n][g1]*p.Unwrap[g1];
                }
            }

            /*  Algorithm 10:  */
            if (i >= p.gap_cap + 1 && i < p.dimensions - 1) {
                p.Q_n[n](i, i + 1) = 1.;
            }

            /*  Algorithm 14:  */
            if (i >= p.sum_mgVg(g1) + 2 && i <= p.sum_mgVg(g1 + 1)) {
                p.Q_n[n](i, i - 1) = 1.;
            }

            for (int j = 0; j < p.Q_n[n].cols(); j++) {
                int g2 = p.find_g(j);
                /*  Algorithm 4:  */
                if (i < p.bound_cap && j == p.sum_mg(g2)) {
                    if (n == 3) {
                        p.Q_n[n](i, j) = p.K_n_g[n][g1]*p.w_g1_g2[g1][g2]*p.Unwrap[g1]*p.C_0_g[g1]*p.C_0_g[g2];
                    }
                    else {
                        p.Q_n[n](i, j) = p.K_n_g[n][g1]*p.w_g1_g2[g1][g2]*p.Unwrap[g1]*p.C_0_g[g1];
                    }
                }

                /*  Algorithm 5:  */
                if (i == p.sum_mg(g1 + 1) - 1 && j < p.bound_cap && j != p.sum_mg(g2)) {
                    p.Q_n[n](i, j) = p.K_n_g[n][g1]*p.w_g1_g2[g1][g2]*p.Unwrap[g2]*p.C_0_g[g2];
                }

                /*  Algorithm 8:  */
                if (j < p.bound_cap) {
                    p.Q_n[n](p.bound_cap, j) = p.w_g1_g2[0][g2]*p.Unwrap[g2]*p.C_0_g[g2];
                }

                /*  Algorithm 12:  */
                if (j < p.bound_cap) {
                    p.Q_n[n](p.dimensions - 1, j) = p.Unwrap[g2]*p.C_0_g[g2];
                }

                /*  Algorithm 13:  */
                if (i < p.bound_cap && (j >= p.bound_cap + 2 && j <= p.gap_cap)) {
                    if (n == 3) {
                        p.Q_n[n](i, j) = p.w_g1_g2[g1][g2]*p.K_n_g[n][g1]*p.Unwrap[g1]*p.C_0_g[g1];
                    }
                    else {
                        p.Q_n[n](i, j) = p.w_g1_g2[g1][g2]*p.K_n_g[n][g1]*p.Unwrap[g1];
                    }
                }

                /*  Algorithm 15:  */
                if (i == p.sum_mgVg(g1) + 1 && j < p.bound_cap) {
                    p.Q_n[n](i, j) = p.Unwrap[g2]*p.C_0_g[g2];
                }
            }
        }
    }
}
