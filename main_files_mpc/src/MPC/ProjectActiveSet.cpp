//
//  ProjectActiveSet.cpp
//  Constrained-LQR
//
//  Created by Georgios on 04/10/15.
//
//

#include <fstream>
#include <cassert>
#include "ProjectActiveSet.hpp"

using namespace arma;
using namespace std;

// Perform projection on the active set
void ProjectActiveSet(data_dims* p_data_dims, data* p_data, vars* p_vars, int& newHor)
{
    // Find active set
    uvec activeMatRow((p_data_dims->px+p_data_dims->pu)*newHor,fill::zeros);
    uvec activeMatCol((p_data_dims->px+p_data_dims->pu)*newHor,fill::zeros);
    int count1 = 0;
    int count2 = 0;
    for (int iii=0; iii<p_data_dims->px+p_data_dims->pu; iii++)
    {
        for (int jjj=0; jjj<newHor; jjj++)
        {
            if (p_vars->lambda(iii,jjj)!=0)
            {
                activeMatRow(count1) = iii;
                activeMatCol(count2) = jjj;
                count1++;   count2++;   
            }
        }
    }
    uint dim = max(find(activeMatRow));
    uvec activeRows = activeMatRow.rows(0,dim);
    uvec activeCols = activeMatCol.rows(0,dim);
    activeRows.print("Nonzero row indices:\n");
    activeCols.print("Nonzero column indices:\n");

    mat temp1 = BlkDiag(p_data->Q, p_data->R);
    mat temp2 = kron(eye<mat>(newHor,newHor),temp1);
    mat E = BlkDiag(temp2,p_data->S);
    E.print("Enormous E matrix\n");

    mat temp3 = join_horiz(-p_data->A, -p_data->B);
    mat temp4 = join_horiz(eye<mat>(p_data_dims->nx,p_data_dims->nx), zeros(p_data_dims->nx,p_data_dims->nu));
    /*
    uvec cols = linspace<uvec>(0,p_temp_data->C.n_cols-1,p_temp_data->C.n_cols);
    mat activeC  = p_temp_data->C.submat(activeVec, cols);

    // Form KKT matrix
    mat tempKKT1 = join_horiz(p_temp_data->H, activeC.t());
    mat tempKKT2 = join_horiz(activeC, zeros(activeVec.n_rows,activeVec.n_rows));
    mat KKT = join_vert(tempKKT1,tempKKT2);
    vec rhsKKT = join_vert(-p_temp_data->h,p_temp_data->d(activeVec));

    // Solve linear system
    vec solKKT = solve(KKT,rhsKKT);
    mat tempu = solKKT(span(0, p_data_dims->m*(newHor-1)));  
    tempu.reshape(p_data_dims->m, newHor); // reshape to matrix form
    p_vars->u(span(0, p_data_dims->m-1), span(0, newHor-1)) = tempu;
    p_vars->x(span::all,span(1, newHor)) = p_data->A*p_vars->x(span::all,span(0, newHor-1)) + p_data->B*p_vars->u(span::all,span(0,newHor-1));
    */
}

// Method for constructing block diagonal matrix
mat BlkDiag(mat& m1, mat& m2)
{
    mat B(m1.n_rows+m2.n_rows,m1.n_cols+m2.n_cols,fill::zeros);
    B(span(0,m1.n_rows-1),span(0,m1.n_cols-1)) = m1;
    B(span(m1.n_rows,m1.n_rows+m2.n_rows-1),span(m1.n_cols,m1.n_cols+m2.n_cols-1)) = m2;
    return B;
}



