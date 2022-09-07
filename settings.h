#pragma once

#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <cmath>
#include <valarray>
#include <fstream>
#include <string>
#include <cctype>
#include <vector>
#include <sstream>
#include <algorithm>
#include "functions.h"
#include <boost/math/constants/constants.hpp>

extern double pi;
using std::valarray;
using std::vector;
using vector2d = vector<vector <double>>;
using vector3d = vector<vector<vector <double>>>;

struct Parameters {
    /**************************
    *   Default values:
    ***************************/
    /*  these values are referred to 10.1103/PhysRevE.94.032404 */
    // B-DNA:
    double b_B = 0.5;       // in [nm] unit;
    double A_B = 50;        // in [nm] unit;
    double C_B = 95;        // in [nm] unit;
    double lambda_B = 4.3;  // λ for B_DNA;
    // L-DNA:
    double b_L = 1.35*b_B;
    double A_L = 7.;
    double C_L = 15.;
    double lambda_L = 4.3;
    double mu_L = 4.9;      // base-pair energy difference μL, in unit of k_B*T;
    double lk_L_0 = -1.0 / 16 - 1 / 10.4;   // Estimated to be - 1/h_L - 1/h_B, nearly -0.16;
    double coop = 9.0;      // energy penalty due to B - L - DNA, B - P - DNA and P - L - DNA boundaries;
    // P-DNA:
    double b_P = 1.7*b_B;
    double A_P = 15.;
    double C_P = 25.;
    double lambda_P = -0.5;
    double mu_P = 17.8;
    double lk_P_0 = 1.0 / 3 - 1 / 10.4;     // Estimated to be 1/h_P - 1/h_B, nearly 0.24;
    // Others:
    double cutoff = pi;
    double plus_coeff = 1.0001;
    double minus_coeff = 0.9999;

	/*	for the row vector on the left-handed side, and the column vector on the right-handed side:  */
	vector<double> left; vector<double> right;

    /**************************
    *   Input variables:
    ***************************/
    /*  For DNA:  */
    double force = 0.1;         // in [pN] unit;
    double torque = 0.;         // in [pN*nm] unit;
    double DNA_length = 3400;   // in [nm] unit;    
    int cutoff_harmonics = 14;  // partition function angular momentum cutoff (10 ~ 15);
    int threads = 1;            // number of concurrent threads;

    /*  For proteins:  */
    int f = 2;                  // total number of types of proteins, default to be 2 (TF & histone octamer);
    vector<int> mg;          	// storing binding length mg of type-g protein;
	vector<int> Vg;				// storing max interaction distance of type-g protein;
	vector2d K_n_g;				// binding constant of type-g protein with binding sites information n;
	vector<double> C_0_g;		// concentration of type-g protein;
	vector<double> Unwrap;		// weight for unwrapping branches of type-g protein;
	vector2d w_g1_g2;			// interaction between type-g1 and type-g2 proteins; 
	

    /**************************
    *   Calculated variables:
    ***************************/
   	int dimensions = 0;				// Dimension for each Q_n(i, j);
	vector3d Q_n;				// Conditional probability matrix for each DNA segment;




    double a_B;                 // B - DNA bending elasticity;
	double c_B;                 // B - DNA twisting elasticity;
	double q;                   // number of DNA base - pairs in individual DNA segment;
	double a_L;                 // L - DNA bending elasticity;
	double c_L;                 // L - DNA twisting elasticity;
	double a_P;                 // P - DNA bending elasticity;
	double c_P;                 // P - DNA twisting elasticity;
	double mu_L_plus;
	double mu_L_minus;
	double mu_P_plus;
	double mu_P_minus;
    int N;                      // total number of segments of DNA;
    double tau;                 // torque / (k_B*T);
    double eff;                 // force / (k_B*T);    
    double chi_B;               // χ = τλ/2πc;
	double chi_L;
	double chi_P;
	double omega_B;             // ω = arctan(χ);
	double omega_L;
	double omega_P;
	double B_coeff;
	double L_coeff;
	double P_coeff;
	double f_plus;
	double f_minus;
	double B_coeff_plus;
	double B_coeff_minus;
	double L_coeff_plus;
	double L_coeff_minus;
	double P_coeff_plus;
	double P_coeff_minus;
	valarray<double> besseli_coeff_B;
	valarray<double> besseli_coeff_plus_B;
	valarray<double> besseli_coeff_minus_B;
	valarray<double> besseli_coeff_L;
	valarray<double> besseli_coeff_plus_L;
	valarray<double> besseli_coeff_minus_L;
	valarray<double> besseli_coeff_P;
	valarray<double> besseli_coeff_plus_P;
	valarray<double> besseli_coeff_minus_P;
	valarray<double> besseli_c_B_torque;
	valarray<double> besseli_c_L_torque;
	valarray<double> besseli_c_P_torque;
	valarray<double> coeff_omega_B_all;
	valarray<double> coeff_omega_L_all;
	valarray<double> coeff_omega_P_all;
	valarray<double> coeff_r_total_B;
	valarray<double> coeff_r_total_B_plus;
	valarray<double> coeff_r_total_B_minus;
	valarray<double> coeff_r_total_L;
	valarray<double> coeff_r_total_L_plus;
	valarray<double> coeff_r_total_L_minus;
	valarray<double> coeff_r_total_P;
	valarray<double> coeff_r_total_P_plus;
	valarray<double> coeff_r_total_P_minus;
};



/*	Part for functions:  */
void set_param (Parameters &p, int argc, char *argv[]);
void set_Q (Parameters &p);
#endif
