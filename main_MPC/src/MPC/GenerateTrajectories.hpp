//
//  GenerateTrajectories.h
//  Constrained-LQR
//
//  Created by Georgios on 09/10/15.
//
//

#ifndef __Constrained_LQR__GenerateTrajectories__
#define __Constrained_LQR__GenerateTrajectories__

#include <stdio.h>
#include </opt/local/include/armadillo>
#include "InitVars.hpp"

using namespace arma;
using namespace std;

//! Header for the GenerateTrajectories.cpp
/*!
 Contains four functions, GeneratePrimalTrajectories, SimDynSys, TV_Riccati and GenerateDualTrajectories
 */


//! function prototypes:

void GeneratePrimalTrajectories(data_dims* p_data_dims, data* p_data, vars* p_vars); // Generate the trajectories ny Riccati recursion

void TV_Riccati(data_dims* p_data_dims, data* p_data, mat& x, mat& u, mat& tempHatLambda, vec& tempHatLambdaf); // (Time-varying) Riccati recursion for state-input computation

void GenerateDualTrajectories(data_dims* p_data_dims, vars* p_vars, mat& grad, mat& gradf, double& rho); // Dual update

#endif /* defined(__Constrained_LQR__GenerateTrajectories__) */
