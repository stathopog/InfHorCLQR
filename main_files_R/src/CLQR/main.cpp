//
//  main.cpp
//  Constrained Linear Quadratic Regulator
//
//  Created by Georgios on 11/10/14.
//  Copyright (c) 2014 Georgios. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <chrono>
#include "ImportDataFromFile.hpp"
#include "SolveCLQR.hpp"
#include "SaveDataToFile.hpp"
#include "OptimalValue.hpp"
#include </opt/local/include/armadillo>

# define PERTURB 0.0 // percentage of initial state perturbation

using namespace std;
using namespace arma;

int main(int argc, const char * argv[])
{
    // Read in data (scalars)
    string filePath = "/Users/georgios/Documents/Research/CLQR_via_FBS/Cpp_Code/Matlab_files/";
    data_dims* p_data_dims = new data_dims;
    ReadInData(p_data_dims, filePath);

    // Read in data (matrices and vectors)
    data* p_data = new data;
    ReadInMatVecs(p_data_dims, filePath, p_data);

    // Initialize vars @ 0
    vars* p_vars = new vars;
    InitVarsCold(p_data_dims, p_vars);

    // Keep track of horizon lengths, factorizations and durations of execution
    uvec HORZS(p_data_dims->N,fill::zeros); // horizon lengths
    vec TIMES(p_data_dims->N,fill::zeros); // durations in ms
    uvec ITERS(p_data_dims->N,fill::zeros); // iterations
    vec OPT_VAL(1,fill::zeros);; // optimal value of cost

    // Closed-loop solution initialize
    mat X(p_vars->x);
    mat U(p_vars->u);
    X.col(0) = p_data->xinit;

    // Main execution
    int k = 1;
    do
    {
        cout << "\n***Solving problem " << k << "..." << "\n\n";

        // main
        SolveCLQR(p_data_dims, p_data, p_vars, HORZS, TIMES, ITERS, k);

        // perturb initial state
        cout << "perturbation is: " << p_data->perturb(k-1) << "\n";
        p_data->xinit = p_vars->x.col(1) + p_data->perturb(k-1)*p_vars->x.col(1);
        p_data->xinit.print("xinit is");

        // Check if initial state in the interior of maximal postitively invariant LQR set and stop
        if (all(p_data->Hf*p_data->xinit-p_data->hf <= 0)) // H*xinit < hf, t=0
        {
            cout << "xinit inside invariant set; use LQ controller \n";
            U.col(k-1) = p_vars->u.col(0);
            X.col(k) = p_data->xinit;
            HORZS(k-1) = 0;
            break;
        }
        else
        {
            U.col(k-1) = p_vars->u.col(0);
            X.col(k) = p_data->xinit;

            // warm-start
            InitVarsWarm(p_data_dims, p_vars, HORZS(k));
        }

        k++;

    } while(HORZS(k-1)>1); // Until x0 in LQR invariant set

    // Closed-loop solution
    X.cols(0,k).print("x =");
    U.cols(0,k-1).print("u =");

    // Compute optimal value
    OPT_VAL = OptimalValue(p_data_dims, p_data, X, U, k);

    // reports
    cout << "\n\nAverage time: " << mean(TIMES.rows(1,k-1)) << " micro-seconds over " << k-1 << " solves.\n";
    cout << "Max time: " << max(TIMES.rows(1,k-1)) << " micro-seconds over " << k-1 << " solves.\n";
    cout << "Average No. of iterations: " << mean(ITERS.rows(1,k-1)) << "\n";
    cout << "Max No. of iterations: " << max(ITERS.rows(1,k-1)) << "\n";
    cout << "The CLQR optimal value is: " << OPT_VAL << "\n";

    // save to file
    WriteData(X, U, filePath, HORZS, TIMES, ITERS, OPT_VAL, k-1);

    return(0);
}
