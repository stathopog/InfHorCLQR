//
//  SolveCLQR.h
//  Constrained-LQR
//
//  Created by Georgios on 13/09/15.
//  Copyright (c) 2015 Georgios. All rights reserved.
//

#ifndef Constrained_LQR_SolveCLQR_h
#define Constrained_LQR_SolveCLQR_h

#include </opt/local/include/armadillo>
//#include "ImportDataFromFile.hpp"
#include "Backtrack.hpp"
#include "GenerateTrajectories.hpp"
#include "ProjectActiveSet.hpp"

using namespace arma;
using namespace std;
using namespace std::chrono;

//! Header for the SolveCLQR.c
/*!
 Contains one function prototype; solution of CLQR problem
 */
 int SolveCLQR(data_dims* p_data_dims, data* p_data, vars* p_vars, uvec& HORZS, vec& TIMES, uvec& ITERS, int& noProb);


#endif /* defined(__Constrained_LQR__SolveCLQR__) */
