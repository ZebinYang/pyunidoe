//
// Created by Administrator on 2019/4/8.
//

#include "doe_WD2.h"

void WD2::init_design(vector<vector<double> > init)
{
    int i,j;
    xc.assign(static_cast<unsigned long long int>(nsamp), vector<double>(nv, 0));
    if(!init.empty())
    {
        for(i=0;i<nsamp;i++) for(j=0;j<nv;j++)
            {
                x[i][j]=(2*init[i][j]-1)/(2*nlevel);
                xc[i][j]=abs(x[i][j]-0.5);
            }
    }
}

void WD2::update_design(vector<vector<double> > init)
{
    int i,j;
    if(!init.empty())
    {
        for(i=0;i<nsamp;i++) for(j=0;j<nv;j++)
            {
                x[i][j]=init[i][j];
                xc[i][j]=abs(x[i][j]-0.5);
            }
    }
}

vector<vector<double> > WD2::get_design()
{
    int i,j;
    for(j=0;j<nv;j++)
    {
        for(i=0;i<nsamp;i++)
            invx[i][j]=x[i][j]*nlevel+0.5;
    }
    return(invx);
}

void WD2::evaluate_criteria()
{
    int i,j,k;
    double part2 = 0;
    D2.assign(static_cast<unsigned long long int>(nsamp), vector<double>(nsamp, 1));
    for(i=0;i<nsamp;i++)
    {
        for(j=i;j<nsamp;j++) {
            for(k=0;k<nv;k++) {
                D2[i][j] *= 3.0/2.0 - abs(x[i][k] - x[j][k]) + pow((x[i][k] - x[j][k]),2);
            }
            D2[j][i] = D2[i][j];
        }
    }
    for(i=0;i<nsamp;i++) {
        for(j=0;j<nsamp;j++) {
            part2 += D2[i][j];
        }
    }
    surrogate_obj = obj = -pow(4.0/3.0, nv) + 1.0/(nsamp*nsamp)*part2;
}

double WD2::columnwise_exchange(int ncol,int ncp, vector<int> idx1,vector<int> idx2)
{
    int i1,i2,i,j,k;
    double diff, diff2 = 0;
    tempx = x; tempD2 = D2;
    for(k=0;k<ncp;k++)
    {
        i1=idx1[k]; i2=idx2[k];
        if (i1>i2) swap(i1,i2);
        for (i=0;i<nsamp;i++)
        {
            for (j=i+1;j<nsamp;j++) {
                if (j == i1) {
                    tempD2[i][j] *= (3.0/2.0 - abs(x[i2][k] - x[i][k]) + pow((x[i2][k] - x[i][k]),2))/
                                    (3.0/2.0 - abs(x[i1][k] - x[i][k]) + pow((x[i1][k] - x[i][k]),2));
                } else if ((j == i2) & (i != i1)) {
                    tempD2[i][j] *= (3.0/2.0 - abs(x[i1][k] - x[i][k]) + pow((x[i1][k] - x[i][k]),2))/
                                    (3.0/2.0 - abs(x[i2][k] - x[i][k]) + pow((x[i2][k] - x[i][k]),2));
                } else if ((i == i1) & (j != i2)) {
                    tempD2[i][j] *= (3.0/2.0 - abs(x[i2][k] - x[j][k]) + pow((x[i2][k] - x[j][k]),2))/
                                    (3.0/2.0 - abs(x[i1][k] - x[j][k]) + pow((x[i1][k] - x[j][k]),2));
                } else if (i == i2) {
                    tempD2[i][j] *= (3.0/2.0 - abs(x[i1][k] - x[j][k]) + pow((x[i1][k] - x[j][k]),2))/
                                    (3.0/2.0 - abs(x[i2][k] - x[j][k]) + pow((x[i2][k] - x[j][k]),2));
                }
                tempD2[j][i] = tempD2[i][j];
                diff2 += 2*(tempD2[i][j] - D2[i][j]);
            }
        }
        swap(tempx[i1][ncol],tempx[i2][ncol]);
    }
    diff = 1.0/(nsamp*nsamp)*diff2;
    return(diff);
}