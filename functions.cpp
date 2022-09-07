#include "functions.h"

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

/*  Returning ∑mg:  */
int sum_mg (vector<int> mg, int g) {
    int res = 0; int index = 0;
    while (index < g) {
        res += mg[index];
        index++;
    }
    return res;
}

/*  Returns the value of current g (For Algorithm Table 1, 2, 3):  */
int find_g_1 (int i, vector<int> mg, int f) {
    int g = -1;
    for (int c = 0; c < f; c++) {
        if (i > sum_mg(mg, c - 1) && i <= sum_mg(mg, c)) {
            g = c;
            break;
        }
    }
    return g;
}

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
