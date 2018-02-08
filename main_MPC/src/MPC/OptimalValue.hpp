//
//  OptimalValue.h
//  Constrained-LQR
//
//  Created by Georgios on 10/10/15.
//
//

#ifndef __Constrained_LQR__OptimalValue__
#define __Constrained_LQR__OptimalValue__

#include <stdio.h>
#include </opt/local/include/armadillo>
#include "InitVars.hpp"

using namespace arma;
using namespace std;

//! Header for the OptimalValue.cpp
/*!
 Contains the function that computes the cost optimal value
 */

//! function prototypes:
double OptimalValue(data_dims* p_data_dims, data* p_data, mat& X, mat& U, int k); // compute optimal value

#endif /* defined(__Constrained_LQR__OptimalValue__) */
