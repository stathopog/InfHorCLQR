//
//  ImportDataFromFile.h
//  Constrained-LQR
//
//  Created by Georgios on 09/10/15.
//
//

#ifndef __Constrained_LQR__ImportDataFromFile__
#define __Constrained_LQR__ImportDataFromFile__

#include <stdio.h>
#include </opt/local/include/armadillo>

using namespace arma;
using namespace std;

//! Header for the ImportDataFromFile.cpp
/*!
 Contains two structures, one for the scalar data, one for the matrix/vector data and two function prototypes
 */

//! struct containing standard problem data dimensions
typedef struct problem_data_dims {
    int nx,nu,px,pu,pf,N;
    double beta, w;
}data_dims;

//! struct containing standard problem data
typedef struct problem_data {
    mat A, B, Q, R, S, K, Hf, C, D, M, L;
    vec xinit, hf, c, d, perturb;
}data;

//! function prototypes:
void ReadInData(data_dims* p_data_dims, string filePath); // reads scalar data provided in .csv files

void ReadInMatVecs(data_dims* p_data_dims, string filePath, data* p_data); // reads matrices and vectors provided in .csv files

#endif /* defined(__Constrained_LQR__ImportDataFromFile__) */
