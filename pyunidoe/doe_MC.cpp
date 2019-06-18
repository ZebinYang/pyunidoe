//
// Created by Administrator on 2019/4/7.
//

#include "doe_MC.h"

void MC::init_design(vector<vector<double> > init)
{
    int i,j;
    A = nsamp*(nlevel*nlevel-1.0)/12.0;
    CORR.assign(static_cast<unsigned long long int>(nv), vector<double>(nv, 0));
    if(!init.empty())
    {
        for(i=0;i<nsamp;i++) for(j=0;j<nv;j++)
            {
                x[i][j]=(2*init[i][j] - (nlevel+1))/(2*sqrt(A));
            }
    }
}

void MC::update_design(vector<vector<double> > init)
{
    int i,j;
    if(!init.empty())
    {
        for(i=0;i<nsamp;i++) for(j=0;j<nv;j++)
            {
                x[i][j]=init[i][j];
            }
    }
    a2->update_design(init);
    a2->evaluate_criteria();
}

vector<vector<double> > MC::get_design()
{
    int i,j;
    for(j=0;j<nv;j++)
    {
        for(i=0;i<nsamp;i++)
            invx[i][j]=(x[i][j]*(2*sqrt(A))+(nlevel+1))/2;
    }
    return(invx);
}

void MC::evaluate_criteria()
{
    int i,j,k,count = 0;
    double a2_surrogate_obj;
    obj = surrogate_obj = 0;
    CORR.assign(static_cast<unsigned long long int>(nv), vector<double>(nv, 0));
    for(i=0;i<nv;i++)
    {
        for(j=i+1;j<nv;j++)
        {
            for (k=0;k<nsamp;k++){
                CORR[i][j]+=x[k][i]*x[k][j];
            }
            if ( (abs(CORR[i][j]) - obj ) > EPS)
            {
                obj = abs(CORR[i][j]);
            }
            CORR[j][i] = CORR[i][j];
        }
    }
    a2->evaluate_criteria();
    a2_surrogate_obj = a2->get_surrogate_criteria();
    surrogate_obj = obj * M + a2_surrogate_obj;
}

double MC::columnwise_exchange(int ncol, int ncp, vector<int> idx1,vector<int> idx2)
{
    int i1,i2,i,j;
    double temp_obj,diff, a2_diff;
    temp_obj = 0;
    tempx = x;
    tempCORR = CORR;
    for(i=0;i<nv;i++)
    {
        for(j=i+1;j<nv;j++)
        {
            if ((i!=ncol)&(j!=ncol)) {
                if ( (abs(tempCORR[i][j]) - temp_obj ) > EPS) {
                    temp_obj = abs(tempCORR[i][j]);
                }
            }
        }
    }
    for(j=0;j<ncp;j++)
    {
        i1=idx1[j]; i2=idx2[j];
        for (i=0;i<nv;i++)
        {
            if (ncol>i) {
                tempCORR[i][ncol] += x[i1][ncol] * x[i2][i] +  x[i2][ncol] * x[i1][i] -  x[i1][ncol] * x[i1][i] -  x[i2][ncol] * x[i2][i];
                if ( (abs(tempCORR[i][ncol]) - temp_obj ) > EPS) {
                    temp_obj = abs(tempCORR[i][ncol]);
                }
            } else if (ncol<i)
            {
                tempCORR[ncol][i] += x[i1][ncol] * x[i2][i] +  x[i2][ncol] * x[i1][i] -  x[i1][ncol] * x[i1][i] -  x[i2][ncol] * x[i2][i];
                if ( (abs(tempCORR[ncol][i]) - temp_obj ) > EPS) {
                    temp_obj = abs(tempCORR[ncol][i]);
                }
            }
        }
        swap(tempx[i1][ncol], tempx[i2][ncol]);
    }
    a2_diff = a2->columnwise_exchange(ncol, ncp, idx1, idx2);
    diff = temp_obj * M + a2_diff - surrogate_obj;
    return(diff);
}

vector<double> MC::get_criteria_matrix()
{
    int i,j;
    vector<double> col_weight(nv, 1);
    for(i=0;i<nv;i++)
    {
        for(j=i+1;j<nv;j++)
        {
            if (abs(abs(CORR[i][j]) - obj) < EPS)
            {
                col_weight[i] +=1;
                col_weight[j] +=1;
            }
        }
    }
    return (col_weight);
}