//
//  SolveMPC.cpp
//  Constrained-LQR
//
//  Created by Georgios on 10/10/15.
//
//

#include <fstream>
#include <cassert>
#include "SolveMPC.hpp"

#define MAXITER 2e4 // maximum number of iterations MPC will perform
#define EPS 1e-12 // tolerance (for testing convergence)
#define a 5 // relaxation parameter for acceleration (alpha = (k-1)/(k+a))

int SolveMPC(data_dims* p_data_dims, data* p_data, vars* p_vars, vec& TIMES, uvec& ITERS, int& noProb)
{
    // Algorithmic initializations
    double rho = 1/p_data_dims->beta;
    double Lip = 1e-2;
    double alpha = 1; // over-relaxation parameter (see Chambolle FISTA)

    // Execution time begin
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

    for (int t=0; t<MAXITER; t++)
    {
        // printf("Iteration %d\t\n", t);

        /*********************************************************************
         *****            Step 1: Generate the trajectories            *******
         *********************************************************************/
        GeneratePrimalTrajectories(p_data_dims, p_data, p_vars);

        /**************************************************************
         *****           Step 2: Backtrack for stepsize         *******
         **************************************************************/
        Backtrack(p_data_dims, p_data, p_vars, Lip);
        rho = 1/Lip;
        mat g_L = zeros(p_data_dims->nu+p_data_dims->nx,p_data_dims->N);
        vec g_Lf = zeros(p_data_dims->pf);
        GradEval(p_data_dims, p_data, p_vars->x, p_vars->u, g_L, g_Lf);

        /**************************************************************
         *****            Step 3: Dual update                   *******
         **************************************************************/
        GenerateDualTrajectories(p_data_dims, p_vars, g_L, g_Lf, rho);

         /**************************************************************
         *****            Step 4: Acceleration                  *******
         **************************************************************/
        alpha = ((double)(t+1) / (double)(t+a+2));
        mat tempLambda = p_vars->lambda.cols(0,p_data_dims->N-1);
        mat tempOldLambda = p_vars->oldlambda.cols(0,p_data_dims->N-1);
        p_vars->hatlambda.cols(0,p_data_dims->N-1) = tempLambda + alpha*(tempLambda-tempOldLambda);
        vec tempLambdaf = p_vars->lambdaf;
        vec tempOldLambdaf = p_vars->oldlambdaf;
        p_vars->hatlambdaf = tempLambdaf + alpha*(tempLambdaf-tempOldLambdaf);

        /**************************************************************
         *****            Step 5: Checking convergence          *******
         **************************************************************/
        if ( (norm(p_vars->lambda - p_vars->oldlambda,"fro") / norm(p_vars->lambda,"fro") <= EPS) || (t == MAXITER-1) || (norm(p_vars->lambda,"fro") < 1e-6) ) // Extra condition; if state inside invariant avoid zero division
        {
            // Projection on the set of active constraints
            //ProjectActiveSet(p_data_dims, p_data, p_vars, newHor);

            // Execution time end
            chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>( t2 - t1 ).count();
            TIMES(noProb) = duration;
            cout <<  "The MPC algorithm needed " << TIMES(noProb) << " micros to converge" "\n";

            // Printing
            cout <<  "Vector u optimal is: " << "\n" << p_vars->u.cols(0, p_data_dims->N-1) << "\n";
            cout <<  "Vector x optimal is: " << "\n" << p_vars->x.cols(0, p_data_dims->N) << "\n";
            ITERS(noProb) = t;
            cout <<  "Total number of iterations: " << ITERS(noProb) << "\n";
            return 0;
        }
    }
    return 0;
}
