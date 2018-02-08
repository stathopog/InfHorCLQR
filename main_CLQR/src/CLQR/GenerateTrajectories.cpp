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
void GeneratePrimalTrajectories(data_dims* p_data_dims, data* p_data, vars* p_vars, int& oldHor, int& newHor, int lastHor, int& T)
{
    // CHECK ASSERTIONS
    if (T == 0)
    {
        if (any(vectorise(p_vars->lambda))) // If warm-starting is on
        {
            newHor = lastHor;
            mat tempHatLambda = p_vars->hatlambda.cols(0, newHor-1); // temporary dual variable
            // Run time-varying Riccati for nonzero multipliers
            //cout << "hatlambda is \n" << tempHatLambda << "\n";
            TV_Riccati(p_data_dims, p_data, p_vars->x, p_vars->u, tempHatLambda, newHor);
            // Simulate the rest
            mat tempu = p_vars->u;
            vec tempxinit = p_vars->x.col(newHor);
            SimDynSys(p_data_dims, p_vars, p_data, tempu, tempxinit, oldHor, newHor, T);
        }
        else // cold start (zero multipliers)
        {
            // Find the first hitting time by forward simulation
            mat tempu = p_vars->u;
            SimDynSys(p_data_dims, p_vars, p_data, tempu, p_data->xinit, oldHor, newHor, T);
        }
    }
    else
    {
        mat tempHatLambda = p_vars->hatlambda.cols(0, newHor-1);
        // Run time-varying Riccati for nonzero multipliers
        TV_Riccati(p_data_dims, p_data, p_vars->x, p_vars->u, tempHatLambda, newHor);
        // Find hitting time
        mat tempu = p_vars->u;
        vec tempxinit = p_vars->x.col(newHor);
        SimDynSys(p_data_dims, p_vars, p_data, tempu, tempxinit, oldHor, newHor, T);
    }
}

/*****************************************************************************
 * Forward simulation of the dynamical system.                               *
 * The state trajectory evolves until it enters a positively invariant set   *
 * H*x(t)<=hf. The function outputs the hitting time 't' and updates         *
 * the global variables x and u.                                             *
 *****************************************************************************/
int SimDynSys(data_dims* p_data_dims,  vars* p_vars, data* p_data, mat& tempu, vec& x0, int& oldHor, int& newHor, int& T)
{
    oldHor = newHor;
    int t = oldHor; // time index
    p_vars->x.col(oldHor) = x0; // initialize state x(0) at xinit
    vec tempcheck = p_data->Hf*x0-p_data->hf; // checking set inclusion

    // Case x0 already in the positively invariant set
    /*
    if (all(p_data->Hf*p_data->xinit-p_data->hf <= 0) && (T==0)) // H*xinit < hf, t=0
    {
        p_vars->u.col(t) = -p_data->K*p_vars->x.col(t);
        p_vars->x.col(t+1) = p_data->A*p_vars->x.col(t) + p_data->B*p_vars->u.col(t);
        newHor = oldHor;
        cout << "xinit inside invariant set; use LQ controller \n";
        exit(1);
    }
    */

    if (all(tempcheck <= 0)) // Hf*xinit < hf
    {
        newHor = oldHor;
    }

    // Other cases
    while (any(tempcheck>0))
    {
        p_vars->u.col(t) = -p_data->K*p_vars->x.col(t);
        p_vars->x.col(t+1) = p_data->A*p_vars->x.col(t) + p_data->B*p_vars->u.col(t);
        t += 1;
        tempcheck = p_data->Hf*p_vars->x.col(t)-p_data->hf;
    }
    int Tstop = t;
    if (Tstop <= oldHor)
    {
        newHor = oldHor;
    }
    else
    {
        oldHor = newHor;
        newHor = Tstop;
    }
    return 0; // return new hitting time
}


/*****************************************************************************
 *                  Riccati Recursion                                        *
 *****************************************************************************/
void TV_Riccati(data_dims* p_data_dims, data* p_data, mat& x, mat& u, mat& tempHatLambda, int& oldHor)
{
    mat k = zeros(p_data_dims->nu, oldHor+1);
    mat p = zeros(p_data_dims->nx, oldHor+1);

    // Perform Riccati recursion for the time-varying part
    for (int t=oldHor-1; t>-1; t--)
    {
        k.col(t) = solve(trimatu(p_data->L.t()), solve(trimatl(p_data->L), -pow(p_data_dims->w,t)*p_data->D.t()*tempHatLambda(span(p_data_dims->px,p_data_dims->px+p_data_dims->pu-1), t) + p_data->B.t()*p.col(t+1) ));
        p.col(t) = -pow(p_data_dims->w,t)*p_data->C.t()*tempHatLambda(span(0,p_data_dims->px-1), t) + p_data->A.t()*p.col(t+1) - p_data->K.t()*p_data->M*k.col(t);
    }

    // Update state and input
    x.col(0) = p_data->xinit;
    for (int t=0; t<oldHor; t++)
    {
        u.col(t) = -p_data->K*x.col(t) - k.col(t);
        x.col(t+1) = p_data->A*x.col(t) + p_data->B*u.col(t);
    }
}

/*****************************************************************************
 *                  Dual update                                              *
 *****************************************************************************/

void GenerateDualTrajectories(data_dims* p_data_dims, vars* p_vars, mat& grad, int& newHor, double& rho)
{
    p_vars->oldlambda.cols(0,newHor-1) = p_vars->lambda.cols(0,newHor-1);
    mat tempHatLambda = p_vars->hatlambda.cols(0,newHor-1);
    mat tempgrad = tempHatLambda - grad*rho;
    p_vars->lambda.cols(0, newHor-1) = arma::min(tempgrad, zeros(p_data_dims->px+p_data_dims->pu,newHor));
    //cout << "lambda is \n" << p_vars->lambda.cols(0, newHor-1) << "\n";
}

