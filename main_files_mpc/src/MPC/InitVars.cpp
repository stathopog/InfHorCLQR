//
//  InitVars.cpp
//  Constrained-LQR
//
//  Created by Georgios on 09/10/15.
//
//

#include <fstream>
#include <cassert>
#include "InitVars.hpp"

using namespace arma;
using namespace std;

void InitVarsCold(data_dims* p_data_dims, vars* p_vars)
{
    p_vars->u = zeros(p_data_dims->nu,p_data_dims->N); // primal variable: input
    p_vars->x = zeros(p_data_dims->nx,p_data_dims->N+1); // primal variable: state
    p_vars->lambda = zeros(p_data_dims->px+p_data_dims->pu,p_data_dims->N); // dual variable
    p_vars->oldlambda = p_vars->lambda; // copy of previous value of the dual
    p_vars->hatlambda = p_vars->lambda;
    p_vars->lambdaf = zeros(p_data_dims->pf); // dual variable for terminal set
    p_vars->oldlambdaf = p_vars->lambdaf;
    p_vars->hatlambdaf = p_vars->lambdaf;
}

void InitVarsWarm(data_dims* p_data_dims, vars* p_vars)
{
    mat initLAM = p_data_dims->w*join_horiz(p_vars->lambda.cols(1,p_data_dims->N-1), zeros(p_data_dims->px+p_data_dims->pu,1));
    p_vars->lambda = initLAM;
    p_vars->oldlambda = p_vars->lambda;
    p_vars->hatlambda = p_vars->lambda;
    p_vars->oldlambdaf = p_vars->lambdaf;
    p_vars->hatlambdaf = p_vars->lambdaf;
}



