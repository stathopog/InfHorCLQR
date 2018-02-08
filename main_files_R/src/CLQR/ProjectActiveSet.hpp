//
//  ProjectActiveSet.h
//  Constrained-LQR
//
//  Created by Georgios on 04/10/15.
//
//

#ifndef __Constrained_LQR__ProjectActiveSet__
#define __Constrained_LQR__ProjectActiveSet__

#include <stdio.h>
#include </opt/local/include/armadillo>
#include "InitVars.hpp"

using namespace arma;
using namespace std;

//! Header for the ProjectActiveSet.cpp
/*!
 Contains the functioin that projects onto the active set; happens once before exiting with the solution
 */


//! function prototypes:

void ProjectActiveSet(data_dims* p_data_dims, data* p_data, vars* p_vars, int& newHor); // Performs projection onto active set

mat BlkDiag(mat& m1, mat& m2); // creates block diagonal sparse matrix from two other matrices

#endif /* defined(__Constrained_LQR__ProjectActiveSet__) */
