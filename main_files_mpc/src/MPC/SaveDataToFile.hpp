//
//  SaveDataToFile.h
//  Constrained-LQR
//
//  Created by Georgios on 10/10/15.
//
//

#ifndef __Constrained_LQR__SaveDataToFile__
#define __Constrained_LQR__SaveDataToFile__

#include <stdio.h>
#include </opt/local/include/armadillo>
#include "InitVars.hpp"

using namespace arma;
using namespace std;

//! Header for the SaveDataToFile.cpp
/*!
 Contains the function that saves the results to a .csv file
 */

//! function prototypes:
void WriteData(mat& X, mat& U, string filePath, vec& TIMES, uvec& ITERS, vec& OPTVAL, int noProb); // saves data to file

#endif /* defined(__Constrained_LQR__SaveDataToFile__) */
