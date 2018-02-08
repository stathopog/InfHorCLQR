//
//  SaveDataToFile.cpp
//  Constrained-LQR
//
//  Created by Georgios on 10/10/15.
//
//

#include <fstream>
#include <cassert>
#include "SaveDataToFile.hpp"

using namespace arma;
using namespace std;


void WriteData(mat& X, mat& U, string filePath, vec& TIMES, uvec& ITERS, vec& OPTVAL, int noProb)
{
	  vec times = TIMES.rows(1,noProb);
    times.save(filePath+"MPC_TIMES.dat",csv_ascii);

    uvec iters = ITERS.rows(1,noProb);
    iters.save(filePath+"MPC_ITERS.dat",csv_ascii);

    mat tempX = X.cols(0,noProb);
    tempX.save(filePath+"MPC_x_opt.dat",csv_ascii);

    mat tempU = U.cols(0,noProb-1);
    tempU.save(filePath+"MPC_u_opt.dat",csv_ascii);

    OPTVAL.save(filePath+"MPC_optval.dat",csv_ascii);
}

