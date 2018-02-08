//
//  SolveCLQR.cpp
//  Constrained-LQR
//
//  Created by Georgios on 10/10/15.
//
//

#include <fstream>
#include <cassert>
#include "SolveCLQR.hpp"

#define MAXITER 2e4 // maximum number of iterations CLQR will perform
#define EPS 1e-12 // tolerance (for testing convergence)
#define a 5 // relaxation parameter for acceleration (alpha = (k-1)/(k+a))

using namespace arma;
using namespace std;

int SolveCLQR(data_dims* p_data_dims, data* p_data, vars* p_vars, uvec& HORZS, vec& TIMES, uvec& ITERS, int& noProb)
{
    // Algorithmic initializations
    double rho = 1/p_data_dims->beta;
    double Lip = 1e-2;
    double alpha = 1; // over-relaxation parameter (see Chambolle FISTA)

    // Initialize hitting time
    int oldHor = HORZS(noProb); int newHor = oldHor;

    // Execution time begin
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    for (int t=0; t<MAXITER; t++)
    {
        // printf("Iteration %d\t\n", t);

        /*********************************************************************
         *****            Step 1: Generate the trajectories            *******
         *********************************************************************/
        GeneratePrimalTrajectories(p_data_dims, p_data, p_vars, oldHor, newHor, HORZS(noProb-1)-1, t);
        int dT = newHor-oldHor;
        // cout << "dT is \n" << dT << "\n";

        /**************************************************************
         *****           Step 2: Backtrack for stepsize         *******
         **************************************************************/
        vec tempw = zeros(newHor);
        // cout << "newHor is \n" << newHor << "\n";
        for (int ttt=0; ttt<newHor; ttt++) {tempw(ttt) = pow(p_data_dims->w,ttt);} // temporary weight matrix
        mat tempW = diagmat(tempw);
        Backtrack(p_data_dims, p_data, p_vars, tempW, newHor, Lip);
        rho = 1/Lip;
        // cout << "rho is: " << rho << "\n";
        mat grad = GradEval(p_data, p_vars->x, p_vars->u, tempW, newHor);

        /**************************************************************
         *****            Step 3: Dual update                   *******
         **************************************************************/
        GenerateDualTrajectories(p_data_dims, p_vars, grad, newHor, rho);


         /**************************************************************
         *****            Step 4: Acceleration                  *******
         **************************************************************/
        alpha = ((double)(t+1) / (double)(t+a+2));
        mat tempLambda = p_vars->lambda.cols(0,newHor-1);
        mat tempOldLambda = p_vars->oldlambda.cols(0,newHor-1);
        p_vars->hatlambda.cols(0,newHor-1) = tempLambda + alpha*(tempLambda-tempOldLambda);

        /**************************************************************
         *****            Step 5: Checking convergence          *******
         **************************************************************/
        if ( (norm(p_vars->lambda - p_vars->oldlambda,"fro") / norm(p_vars->lambda,"fro") <= EPS) || (t == MAXITER-1) )
        {
            // Projection on the set of active constraints
            // ProjectActiveSet(p_data_dims, p_data, p_vars, newHor);

            // Recover optimal horizon length
            rowvec index = arma::min(p_vars->lambda.cols(0, newHor-1),0);
            //index.print("index is:");
            if (index(newHor-1)<0)
            {
                HORZS(noProb) = newHor;
            }
            else
            {
                uvec temp1 = find(index!=0);
                HORZS(noProb) = max(temp1)+1;
            }

            // Execution time end
            high_resolution_clock::time_point t2 = high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
            TIMES(noProb) = duration;
            cout <<  "The CLQR algorithm needed " << TIMES(noProb) << " micros to converge" "\n";

            // Printing
            std::cout <<  "Vector u optimal is: " << "\n" << p_vars->u.cols(0, HORZS(noProb)-1) << "\n";
            std::cout <<  "Vector x optimal is: " << "\n" << p_vars->x.cols(0, HORZS(noProb)) << "\n";
            ITERS(noProb) = t;
            std::cout <<  "Identified horizon length T is: " << HORZS(noProb) << "\n";
            std::cout <<  "Total number of iterations: " << ITERS(noProb) << "\n";
            return 0;
        }
    }
    return 0;
}
