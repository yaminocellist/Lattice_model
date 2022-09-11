#include "functions.h"
using vector2d = vector<vector <double>>;
const double epsilon = 1e-5;

void printvd (vector<double> vd) {
    for (int i = 0; i < vd.size(); i++) {
        std::cout << vd[i] << std::endl;
    }
}

void printvd (vector<int> vd) {
    for (int i = 0; i < vd.size(); i++) {
        std::cout << vd[i] << std::endl;
    }
}

void print2dv (vector2d v) {
    for (int i = 0; i < v.size(); i++) {
        for (int j = 0; j < v[i].size(); j++) {
            std::cout << v[i][j] << "  ";
        }
        std::cout << std::endl;
    }
}

void print3dv (vector3d v) {
    for (int i = 0; i < v.size(); i++) {
        for (int j = 0; j < v[i].size(); j++) {
            for (int k = 0; k < v[i][j].size(); k++) {
                std::cout << v[i][j][k] << "   ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

/*  Returning ∑mg:  */
// int sum_mg (vector<int> mg, int g) {
//     int res = 0; int index = 0;
//     while (index < g) {
//         res += mg[index];
//         index++;
//     }
//     return res;
// }

/*  Returns the value of current g (For Algorithm Table 1, 2, 3):  */
// int find_g_1 (int i, vector<int> mg, int f) {
//     int g = 0;
//     for (int c = 0; c < f; c++) {
//         if (i >= sum_mg(mg, c) && i <= sum_mg(mg, c + 1) - 1) {
//             g = c;
//             break;
//         }
//     }
//     return g;
// }

/*  Find the first non-zero element of the vector, return its index:  */
int first_non_zero (vector<double> v) {
    int ind;
    for (int i = 0; i < v.size(); i++) {
        if (std::abs(v[i] - 0.) > epsilon) {
            ind = i;
            break;
        }
    }
    return ind;
}

/*  Return rest elements of 'mg' vector other than an appointed value:  */
vector<int> rest_mg (vector<int> mg, int g) {
    mg.erase(std::remove(mg.begin(), mg.end(), g), mg.end());
    return mg;
}

vector<MatrixXd> set_Qn_delta(Parameters p, int n, int g) {
    vector2d Kng_temp = p.K_n_g;
    Kng_temp[n][g] += Kng_temp[n][g] * p.delta;
    std::vector<MatrixXd> Qn_delta;
    for (int i = 0; i < p.N; i++) {
        Qn_delta.push_back(MatrixXd::Zero(p.dimensions, p.dimensions));
    }
    for (int n = 0; n < Qn_delta.size(); n++) {
        /*  Algorithm 6:  */
        Qn_delta[n](p.bound_cap, p.bound_cap) = 1.;

        /*  Algorithm 7:  */
        Qn_delta[n](p.bound_cap + 1, p.bound_cap + 1) = 1.;
            
        /*  Algorithm 11:  */
        Qn_delta[n](p.dimensions - 1, p.dimensions - 1) = 1.;

        for (int i = 0; i < Qn_delta[n].rows(); i++) {
            int g1 = p.find_g(i);
            /*  Algorithm 1:  */
            if (i == p.sum_mg(g1)) {
                if (n == 3) {
                    Qn_delta[n](i, i+1) = Kng_temp[n][g1]*p.C_0_g[g1];
                }
                else {
                    Qn_delta[n](i, i+1) = Kng_temp[n][g1];
                }
            }

            /*  Algorithm 2:  */
            if (i < p.sum_mg(g1 + 1) - 1 && i > p.sum_mg(g1)) {
                Qn_delta[n](i, i+1) = Kng_temp[n][g1];
            }

            /*  Algorithm 3:  */
            if (i < p.bound_cap) {
                if (n == 3) {
                    Qn_delta[n](i, p.bound_cap + 1) = Kng_temp[n][g1]*p.Unwrap[g1]*p.C_0_g[g1];
                }
                else {
                    Qn_delta[n](i, p.bound_cap + 1) = Kng_temp[n][g1]*p.Unwrap[g1];
                }
            }

            /*  Algorithm 9:  */
            if (i < p.bound_cap) {
                if (n == 3) {
                    Qn_delta[n](i, p.gap_cap + 1) = p.w_g1_g2[g1][0]*Kng_temp[n][g1]*p.Unwrap[g1]*p.C_0_g[g1];
                }
                else {
                    Qn_delta[n](i, p.gap_cap + 1) = p.w_g1_g2[g1][0]*Kng_temp[n][g1]*p.Unwrap[g1];
                }
            }

            /*  Algorithm 10:  */
            if (i >= p.gap_cap + 1 && i < p.dimensions - 1) {
                Qn_delta[n](i, i + 1) = 1.;
            }

            /*  Algorithm 14:  */
            if (i >= p.sum_mgVg(g1) + 2 && i <= p.sum_mgVg(g1 + 1)) {
                Qn_delta[n](i, i - 1) = 1.;
            }

            for (int j = 0; j < Qn_delta[n].cols(); j++) {
                int g2 = p.find_g(j);
                /*  Algorithm 4:  */
                if (i < p.bound_cap && j == p.sum_mg(g2)) {
                    if (n == 3) {
                        Qn_delta[n](i, j) = Kng_temp[n][g1]*p.w_g1_g2[g1][g2]*p.Unwrap[g1]*p.C_0_g[g1]*p.C_0_g[g2];
                    }
                    else {
                        Qn_delta[n](i, j) = Kng_temp[n][g1]*p.w_g1_g2[g1][g2]*p.Unwrap[g1]*p.C_0_g[g1];
                    }
                }

                /*  Algorithm 5:  */
                if (i == p.sum_mg(g1 + 1) - 1 && j < p.bound_cap && j != p.sum_mg(g2)) {
                    Qn_delta[n](i, j) = Kng_temp[n][g1]*p.w_g1_g2[g1][g2]*p.Unwrap[g2]*p.C_0_g[g2];
                }

                /*  Algorithm 8:  */
                if (j < p.bound_cap) {
                    Qn_delta[n](p.bound_cap, j) = p.w_g1_g2[0][g2]*p.Unwrap[g2]*p.C_0_g[g2];
                }

                /*  Algorithm 12:  */
                if (j < p.bound_cap) {
                    Qn_delta[n](p.dimensions - 1, j) = p.Unwrap[g2]*p.C_0_g[g2];
                }

                /*  Algorithm 13:  */
                if (i < p.bound_cap && (j >= p.bound_cap + 2 && j <= p.gap_cap)) {
                    if (n == 3) {
                        Qn_delta[n](i, j) = p.w_g1_g2[g1][g2]*Kng_temp[n][g1]*p.Unwrap[g1]*p.C_0_g[g1];
                    }
                    else {
                        Qn_delta[n](i, j) = p.w_g1_g2[g1][g2]*Kng_temp[n][g1]*p.Unwrap[g1];
                    }
                }

                /*  Algorithm 15:  */
                if (i == p.sum_mgVg(g1) + 1 && j < p.bound_cap) {
                    Qn_delta[n](i, j) = p.Unwrap[g2]*p.C_0_g[g2];
                }
            }
        }
    }
    return Qn_delta;
}