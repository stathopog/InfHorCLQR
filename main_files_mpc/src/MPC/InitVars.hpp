//
//  InitVars.h
//  Constrained-LQR
//
//  Created by Georgios on 09/10/15.
//
//

#ifndef __Constrained_LQR__InitVars__
#define __Constrained_LQR__InitVars__

#include <stdio.h>
#include </opt/local/include/armadillo>
#include "ImportDataFromFile.hpp"

using namespace arma;
using namespace std;

//! Header for the InitVars.cpp
/*!
 Contains a structure for the variables (u,x,lambda,oldlambda,hatlambda) of the algorithm and two functions to cold or warm-start
 */

//! variables' structure
typedef struct problem_vars {
    mat u,x; // primals
    mat lambda, oldlambda, hatlambda; // duals
    vec lambdaf, oldlambdaf, hatlambdaf; // terminal dual
}vars;

//! function prototypes:
void InitVarsCold(data_dims* p_data_dims, vars* p_vars); // cold-start

void InitVarsWarm(data_dims* p_data_dims, vars* p_vars); // warm-start

#endif /* defined(__Constrained_LQR__InitVars__) */
