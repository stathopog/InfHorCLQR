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
void Backtrack(data_dims* p_data_dims, data* p_data, vars* p_vars, double& Lip)
{
    // CHECK ASSERTIONS
    double eta = 1.2; // Lipschitz contraction constant
    mat tempHatLambda = p_vars->hatlambda.cols(0, p_data_dims->N-1); // temporary dual variable
    vec tempHatLambdaf = p_vars->hatlambdaf; // temporary dual for terminal constraint
    mat tempx = p_vars->x; // temporary primals
    mat tempu = p_vars->u;
       // cout << "hatlambda is \n" << tempHatLambda;
    mat g_L = zeros(p_data_dims->nu+p_data_dims->nx,p_data_dims->N);
    vec g_Lf = zeros(p_data_dims->pf);
    GradEval(p_data_dims, p_data, tempx, tempu, g_L, g_Lf); // @ hat.lambda
    // cout << "grad is \n" << g_Lf;
    double h = FunEval(p_data_dims, p_data, tempx, tempu, tempHatLambda, tempHatLambdaf); // @ hat.lambda
    //cout << "obj is \n" << h << "\n";
    mat tempzero(tempHatLambda.n_rows, tempHatLambda.n_cols, fill::zeros);
    mat p_L = arma::min(tempHatLambda-(1/Lip)*g_L, tempzero);
    vec p_Lf = arma::min(tempHatLambdaf-(1/Lip)*g_Lf, zeros(p_data_dims->pf));
    // cout << "p_Lf is \n" << p_Lf << "\n";
    TV_Riccati(p_data_dims, p_data, tempx, tempu, p_L, p_Lf);
     //   cout << "x is \n" << p_vars->x.cols(0,newHor) << "\n";
     //       cout << "u is \n" << p_vars->u.cols(0,newHor-1) << "\n";
    double F = FunEval(p_data_dims, p_data, tempx, tempu, p_L, p_Lf); // @ p_L
    double Q_L = h + trace((p_L-tempHatLambda).t()*g_L)  + 0.5*Lip * trace((p_L-tempHatLambda).t()*(p_L-tempHatLambda)) + dot(p_Lf-tempHatLambdaf,g_Lf) + 0.5*Lip*dot(p_Lf-tempHatLambdaf,p_Lf-tempHatLambdaf);
    int kkk = 0;
    double diff = F - Q_L - 1e-8;
    while ( (diff>=0) && (Lip <= p_data_dims->beta) && (kkk <= 100) )
    {
        Lip = eta*Lip;
        mat p_L = arma::min(tempHatLambda-(1/Lip)*g_L, tempzero);
        vec p_Lf = arma::min(tempHatLambdaf-(1/Lip)*g_Lf, zeros(p_data_dims->pf));
        TV_Riccati(p_data_dims, p_data, tempx, tempu, p_L, p_Lf);
        double F = FunEval(p_data_dims, p_data, tempx, tempu, p_L, p_Lf); // @ p_L
        //cout << "F is: " << F << "\n";
        double Q_L = h + trace((p_L-tempHatLambda).t()*g_L)  + 0.5*Lip * trace((p_L-tempHatLambda).t()*(p_L-tempHatLambda)) + dot(p_Lf-tempHatLambdaf,g_Lf) + 0.5*Lip*dot(p_Lf-tempHatLambdaf,p_Lf-tempHatLambdaf);
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
double FunEval(data_dims* p_data_dims, data* p_data, mat& tempx, mat& tempu, mat& templambda, vec& templambdaf)
{
    double h1 = 0.5*trace(tempx.cols(0,p_data_dims->N-1).t()*p_data->Q*tempx.cols(0,p_data_dims->N-1)) + 0.5*trace(tempu.t()*p_data->R*tempu);
    double h2 = -trace(join_vert(p_data->C*tempx.cols(0,p_data_dims->N-1)-repmat(p_data->c,1,p_data_dims->N), p_data->D*tempu-repmat(p_data->d,1,p_data_dims->N)).t() * templambda);
    double h = - (h1 + h2 + 0.5*trace(tempx.col(p_data_dims->N).t()*p_data->S*tempx.col(p_data_dims->N)) - dot(templambdaf, p_data->Hf*tempx.col(p_data_dims->N) - p_data->hf) );
    return h;
}

// Gradient evaluation
void GradEval(data_dims* p_data_dims, data* p_data, mat& tempx, mat& tempu, mat& gradL, vec& gradLf)
{
   //cout << "tempx is: \n" << tempx.cols(0,newHor) << "\n";
   //cout << "tempu is: \n" << tempu.cols(0,newHor-1) << "\n";
   //cout << "W is: \n" << tempW << "\n";
   gradL  = join_vert( p_data->C*tempx.cols(0,p_data_dims->N-1)-repmat(p_data->c,1,p_data_dims->N), p_data->D*tempu.cols(0,p_data_dims->N-1)-repmat(p_data->d,1,p_data_dims->N) );
   gradLf = p_data->Hf*tempx.col(p_data_dims->N) - p_data->hf;
}

