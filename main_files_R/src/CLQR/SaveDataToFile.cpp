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


void WriteData(mat& X, mat& U, string filePath, uvec& HORZS, vec& TIMES, uvec& ITERS, vec& OPTVAL, int noProb)
{
    vec times = TIMES.rows(1,noProb);
    times.save(filePath+"CLQR_TIMES.dat",csv_ascii);

    uvec horzs = HORZS.rows(1,noProb)+1;
    horzs.save(filePath+"CLQR_HORZS.dat",csv_ascii);

    uvec iters = ITERS.rows(1,noProb);
    iters.save(filePath+"CLQR_ITERS.dat",csv_ascii);

    mat tempX = X.cols(0,noProb);
    tempX.save(filePath+"CLQR_x_opt.dat",csv_ascii);

    mat tempU = U.cols(0,noProb-1);
    tempU.save(filePath+"CLQR_u_opt.dat",csv_ascii);

    OPTVAL.save(filePath+"CLQR_optval.dat",csv_ascii);
}

