//
//  ImportDataFromFile.cpp
//  Constrained-LQR
//
//  Created by Georgios on 09/10/15.
//
//

#include <fstream>
#include <cassert>
#include "ImportDataFromFile.hpp"

using namespace arma;
using namespace std;

void ReadInData(data_dims* p_data_dims, string filePath)
{
    ifstream read_file;
    string temp = filePath + "sizes_data";
    cout << temp << "\n";
    read_file.open(temp);
    assert(read_file.is_open());
    read_file >> p_data_dims->nx >> p_data_dims->nu >> p_data_dims->px >> p_data_dims->pu >> p_data_dims->pf >> p_data_dims->N >> p_data_dims->beta >> p_data_dims->w;
    read_file.close();
    // p_data_dims->w = 1;
    // p_data_dims->N = 17;
}


void ReadInMatVecs(data_dims* p_data_dims, string filePath, data* p_data)
{
    // Matrices
    p_data->A.load(filePath+"A_data.dat",csv_ascii);
    //p_data->A.print("A =");

    p_data->B.load(filePath+"B_data.dat",csv_ascii);

    p_data->Q.load(filePath+"Q_data.dat",csv_ascii);

    p_data->R.load(filePath+"R_data.dat",csv_ascii);

    p_data->S.load(filePath+"S_data.dat",csv_ascii);

    p_data->M.load(filePath+"M_data.dat",csv_ascii); // The matrix R+B'SB

    p_data->L.load(filePath+"L_data.dat",csv_ascii); // Lower Cholesky factor of inv(R+B'SB)

    p_data->K.load(filePath+"K_data.dat",csv_ascii);

    p_data->Hf.load(filePath+"H_data.dat",csv_ascii);

    p_data->C.load(filePath+"C_data.dat",csv_ascii); // State constraints Cx<c

    p_data->D.load(filePath+"D_data.dat",csv_ascii); // Input constraints Du<d

    // Vectors
    p_data->xinit.load(filePath+"xinit_data.dat",csv_ascii);
    //p_data->xinit.print("xinit =");

    p_data->hf.load(filePath+"hf_data.dat",csv_ascii);

    p_data->c.load(filePath+"cx_data.dat",csv_ascii);

    p_data->d.load(filePath+"cu_data.dat",csv_ascii);

    p_data->perturb.load(filePath+"perturb.dat",csv_ascii);
}

