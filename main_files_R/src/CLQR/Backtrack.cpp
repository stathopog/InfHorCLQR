//
//  Backtrack.cpp
//  Constrained-LQR
//
//  Created by Georgios on 09/10/15.
//
//

#include <fstream>
#include <cassert>
#include "Backtrack.hpp"

using namespace arma;
using namespace std;

// Perform the backtracking and return a local Lipschitz
void Backtrack(data_dims* p_data_dims, data* p_data, vars* p_vars, mat& tempW, int& newHor, double& Lip)
{
    // CHECK ASSERTIONS
    double eta = 1.2; // Lipschitz contraction constant
    mat tempHatLambda = p_vars->hatlambda.cols(0, newHor-1); // temporary dual variable
    mat tempx = p_vars->x.cols(0, newHor); // temporary primals
    mat tempu = p_vars->u.cols(0, newHor-1);
       // cout << "hatlambda is \n" << tempHatLambda;
    mat g = GradEval(p_data, tempx, tempu, tempW, newHor); // @ hat.lambda
    //cout << "grad is \n" << g;
    double h = FunEval(p_data, tempx, tempu, tempW, tempHatLambda, newHor); // @ hat.lambda
    //cout << "obj is \n" << h << "\n";
    mat tempzero(tempHatLambda.n_rows, tempHatLambda.n_cols, fill::zeros);
    mat p_L = arma::min(tempHatLambda-(1/Lip)*g, tempzero);
    //cout << "p_L is \n" << p_L << "\n";
    TV_Riccati(p_data_dims, p_data, tempx, tempu, p_L, newHor);
     //   cout << "x is \n" << p_vars->x.cols(0,newHor) << "\n";
     //       cout << "u is \n" << p_vars->u.cols(0,newHor-1) << "\n";
    double F = FunEval(p_data, tempx, tempu, tempW, p_L, newHor); // @ p_L
    double Q_L = h + trace((p_L-tempHatLambda).t()*(g*tempW)) + 0.5*Lip * trace(((p_L-tempHatLambda).t()*(p_L-tempHatLambda))*tempW);
    int kkk = 0;
    double diff = F - Q_L - 1e-8;
    while ( (diff>=0) && (Lip <= p_data_dims->beta) && (kkk <= 100) )
    {
        Lip = eta*Lip;
        mat p_L = arma::min(tempHatLambda-(1/Lip)*g, tempzero);
        TV_Riccati(p_data_dims, p_data, tempx, tempu, p_L, newHor);
        double F = FunEval(p_data, tempx, tempu, tempW, p_L, newHor); // @ p_L
        //cout << "F is: " << F << "\n";
        double Q_L = h + trace((p_L-tempHatLambda).t()*(g*tempW)) + 0.5 * Lip * trace(((p_L-tempHatLambda).t()*(p_L-tempHatLambda))*tempW);
        //cout << "Q_L is: " << Q_L << "\n";
        if ( Lip >= p_data_dims->beta )
        {
            Lip = p_data_dims->beta;
        }
        kkk +=1;
        diff = F - Q_L;
    }
}

// Function evaluation
double FunEval(data* p_data, mat& tempx, mat& tempu, mat& tempW, mat& templambda, int& newHor)
{
    double h1 = 0.5*trace(tempx.cols(0,newHor-1).t()*p_data->Q*tempx.cols(0,newHor-1)) + 0.5*trace(tempu.t()*p_data->R*tempu);
    double h2 = -trace(join_vert(p_data->C*tempx.cols(0,newHor-1)-repmat(p_data->c,1,newHor), p_data->D*tempu-repmat(p_data->d,1,newHor)).t() * (templambda*tempW));
    double h = - (h1 + h2 + 0.5*trace(tempx.col(newHor).t()*p_data->S*tempx.col(newHor)) );
    return h;
}

// Gradient evaluation
mat GradEval(data* p_data, mat& tempx, mat& tempu, mat& tempW, int& newHor)
{
   //cout << "tempx is: \n" << tempx.cols(0,newHor) << "\n";
   //cout << "tempu is: \n" << tempu.cols(0,newHor-1) << "\n";
   //cout << "W is: \n" << tempW << "\n";
   mat g = join_vert( p_data->C*tempx.cols(0,newHor-1)-repmat(p_data->c,1,newHor), p_data->D*tempu.cols(0,newHor-1)-repmat(p_data->d,1,newHor) ) * tempW;
   return g;
}

