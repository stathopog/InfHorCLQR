//
//  OptimalValue.cpp
//  Constrained-LQR
//
//  Created by Georgios on 14/05/16.
//
//

#include <fstream>
#include <cassert>
#include "OptimalValue.hpp"

using namespace arma;
using namespace std;


double OptimalValue(data_dims* p_data_dims, data* p_data, mat& X, mat& U, int k)
{
  double value = 0.5*trace(X.cols(0,k-2).t()*p_data->Q*X.cols(0,k-2)) + 0.5*trace(U.t()*p_data->R*U);
  value += 0.5*trace(X.col(k-1).t()*p_data->S*X.col(k-1));
  return value;
}

