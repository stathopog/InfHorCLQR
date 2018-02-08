//
//  GenerateTrajectories.cpp
//  Constrained-LQR
//
//  Created by Georgios on 09/10/15.
//
//

#include <fstream>
#include <cassert>
#include <string>
#include "GenerateTrajectories.hpp"

using namespace arma;
using namespace std;

// Perform the backtracking and return a local Lipschitz
void GeneratePrimalTrajectories(data_dims* p_data_dims, data* p_data, vars* p_vars)
{

    mat tempHatLambda = p_vars->hatlambda.cols(0, p_data_dims->N-1);
    vec tempHatLambdaf = p_vars->hatlambdaf;
    // Run time-varying Riccati for nonzero multipliers
    TV_Riccati(p_data_dims, p_data, p_vars->x, p_vars->u, tempHatLambda, tempHatLambdaf);


}


/*****************************************************************************
 *                  Riccati Recursion                                        *
 *****************************************************************************/
void TV_Riccati(data_dims* p_data_dims, data* p_data, mat& x, mat& u, mat& tempHatLambda, vec& tempHatLambdaf)
{
    mat k = zeros(p_data_dims->nu, p_data_dims->N+1);
    mat p = join_horiz(zeros(p_data_dims->nx, p_data_dims->N), -p_data->Hf.t()*tempHatLambdaf);

    // Perform Riccati recursion for the time-varying part
    for (int t=p_data_dims->N-1; t>-1; t--)
    {
        k.col(t) = solve(trimatu(p_data->L.t()), solve(trimatl(p_data->L), -pow(p_data_dims->w,t)*p_data->D.t()*tempHatLambda(span(p_data_dims->px,p_data_dims->px+p_data_dims->pu-1), t) + p_data->B.t()*p.col(t+1) ));
        p.col(t) = -pow(p_data_dims->w,t)*p_data->C.t()*tempHatLambda(span(0,p_data_dims->px-1), t) + p_data->A.t()*p.col(t+1) - p_data->K.t()*p_data->M*k.col(t);
    }

    // Update state and input
    x.col(0) = p_data->xinit;
    for (int t=0; t<p_data_dims->N; t++)
    {
        u.col(t) = -p_data->K*x.col(t) - k.col(t);
        x.col(t+1) = p_data->A*x.col(t) + p_data->B*u.col(t);
    }
}

/*****************************************************************************
 *                  Dual update                                              *
 *****************************************************************************/

void GenerateDualTrajectories(data_dims* p_data_dims, vars* p_vars, mat& grad, mat& gradf, double& rho)
{
    p_vars->oldlambda.cols(0,p_data_dims->N-1) = p_vars->lambda.cols(0,p_data_dims->N-1);
    p_vars->oldlambdaf = p_vars->lambdaf;
    mat tempHatLambda = p_vars->hatlambda.cols(0,p_data_dims->N-1);
    vec tempHatLambdaf = p_vars->hatlambdaf;
    mat tempgrad = tempHatLambda - grad*rho;
    vec tempgradf = tempHatLambdaf - gradf*rho;
    p_vars->lambda.cols(0, p_data_dims->N-1) = arma::min(tempgrad, zeros(p_data_dims->px+p_data_dims->pu,p_data_dims->N));
    p_vars->lambdaf = arma::min(tempgradf, zeros(p_data_dims->pf));
    //cout << "lambda is \n" << p_vars->lambda.cols(0, newHor-1) << "\n";
}

