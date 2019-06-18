//
// Created by Administrator on 2019/4/7.
//

#include "doe_A2.h"

void A2::init_design(vector<vector<double> > init)
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

void A2::update_design(vector<vector<double> > init)
{
    int i,j;
    if(!init.empty())
    {
        for(i=0;i<nsamp;i++) for(j=0;j<nv;j++)
            {
                x[i][j]=init[i][j];
            }
    }
}

vector<vector<double> > A2::get_design()
{
    int i,j;
    for(j=0;j<nv;j++)
    {
        for(i=0;i<nsamp;i++)
            invx[i][j]=(x[i][j]*(2*sqrt(A))+(nlevel+1))/2;
    }
    return(invx);
}

void A2::evaluate_criteria()
{
    int i,j,k;
    obj = 0;
    CORR.assign(static_cast<unsigned long long int>(nv), vector<double>(nv, 0));
    for(i=0;i<nv;i++)
    {
        for(j=i+1;j<nv;j++)
        {
            for (k=0;k<nsamp;k++){
                CORR[i][j]+=x[k][i]*x[k][j];
            }
            obj += pow(CORR[i][j], 2);
        }
    }
    obj /=(nv*(nv-1)/2);
    surrogate_obj = obj;
}

double A2::columnwise_exchange(int ncol, int ncp, vector<int> idx1,vector<int> idx2)
{
    int i1,i2,i,j;
    double temp_obj,diff;
    temp_obj = surrogate_obj;
    tempx = x; tempCORR = CORR;
    for(j=0;j<ncp;j++)
    {
        i1=idx1[j]; i2=idx2[j];
        for (i=0;i<nv;i++)
        {
            if (ncol>i) {
                tempCORR[i][ncol] += x[i1][ncol] * x[i2][i] +  x[i2][ncol] * x[i1][i] -  x[i1][ncol] * x[i1][i] -  x[i2][ncol] * x[i2][i];
            } else if (ncol<i)
            {
                tempCORR[ncol][i] += x[i1][ncol] * x[i2][i] +  x[i2][ncol] * x[i1][i] -  x[i1][ncol] * x[i1][i] -  x[i2][ncol] * x[i2][i];
            }
            temp_obj += (pow(tempCORR[i][ncol], 2) - pow(CORR[i][ncol], 2))/(nv*(nv-1)/2);
        }
        swap(tempx[i1][ncol], tempx[i2][ncol]);
    }
    diff = temp_obj - surrogate_obj;
    return(diff);
}