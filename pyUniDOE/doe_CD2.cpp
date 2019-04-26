//
// Created by Administrator on 2019/4/8.
//
#include "doe_CD2.h"

void CD2::init_design(vector<vector<double> > init)
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

void CD2::update_design(vector<vector<double> > init)
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

vector<vector<double> > CD2::get_design()
{
    int i,j;
    for(j=0;j<nv;j++)
    {
        for(i=0;i<nsamp;i++)
        {
            invx[i][j]=x[i][j]*nlevel+0.5;
        }
    }
    return(invx);
}

void CD2::evaluate_criteria()
{
    int i,j,k;
    double part1 = 0, part2 = 0;
    D1.assign(nsamp, 1);
    D2.assign(static_cast<unsigned long long int>(nsamp), vector<double>(nsamp, 1));
    for(i=0;i<nsamp;i++)
    {
        for(k=0;k<nv;k++) {
            D1[i] *= 1 + 0.5 * xc[i][k] - 0.5 * xc[i][k] * xc[i][k];
        }
        for(j=i;j<nsamp;j++) {
            for(k=0;k<nv;k++) {
                D2[i][j] *= (1 + 0.5 * xc[i][k] + 0.5 * xc[j][k] - 0.5 * abs(x[i][k] - x[j][k]));
            }
            D2[j][i] = D2[i][j];
        }
    }
    for(i=0;i<nsamp;i++) {
        part1 += D1[i];
        for(j=0;j<nsamp;j++) {
            part2 += D2[i][j];
        }
    }
    surrogate_obj = obj = pow(13.0/12.0, nv) - 2.0/nsamp*part1 + 1.0/(nsamp*nsamp)*part2;
}

double CD2::columnwise_exchange(int ncol,int ncp, vector<int> idx1,vector<int> idx2)
{
    int i1,i2,i,j,k;
    double diff, diff1 = 0, diff2 = 0;
    tempx = x; tempD1 = D1; tempD2 = D2;
    for(k=0;k<ncp;k++)
    {
        i1=idx1[k]; i2=idx2[k];
        if (i1>i2) swap(i1,i2);
        for (i=0;i<nsamp;i++)
        {
            for (j=i+1;j<nsamp;j++) {
                if (j == i1) {
                    tempD2[i][j] *= (1+0.5*xc[i2][ncol]+0.5*xc[i][ncol]-0.5*abs(x[i][ncol]-x[i2][ncol]))/(1+0.5*xc[i1][ncol]+0.5*xc[i][ncol]-0.5*abs(x[i][ncol]-x[i1][ncol]));
                } else if ((j == i2) & (i != i1)) {
                    tempD2[i][j] *= (1+0.5*xc[i1][ncol]+0.5*xc[i][ncol]-0.5*abs(x[i][ncol]-x[i1][ncol]))/(1+0.5*xc[i2][ncol]+0.5*xc[i][ncol]-0.5*abs(x[i][ncol]-x[i2][ncol]));
                } else if ((i == i1) & (j != i2)) {
                    tempD2[i][j] *= (1+0.5*xc[i2][ncol]+0.5*xc[j][ncol]-0.5*abs(x[j][ncol]-x[i2][ncol]))/(1+0.5*xc[i1][ncol]+0.5*xc[j][ncol]-0.5*abs(x[j][ncol]-x[i1][ncol]));
                } else if (i == i2) {
                    tempD2[i][j] *= (1+0.5*xc[i1][ncol]+0.5*xc[j][ncol]-0.5*abs(x[j][ncol]-x[i1][ncol]))/(1+0.5*xc[i2][ncol]+0.5*xc[j][ncol]-0.5*abs(x[j][ncol]-x[i2][ncol]));
                }
                tempD2[j][i] = tempD2[i][j];
                diff2 += 2*(tempD2[i][j] - D2[i][j]);
            }
        }
        tempD1[i1] *= (1+0.5*xc[i2][ncol]-0.5*xc[i2][ncol]*xc[i2][ncol])/(1+0.5*xc[i1][ncol]-0.5*xc[i1][ncol]*xc[i1][ncol]);
        tempD1[i2] *= (1+0.5*xc[i1][ncol]-0.5*xc[i1][ncol]*xc[i1][ncol])/(1+0.5*xc[i2][ncol]-0.5*xc[i2][ncol]*xc[i2][ncol]);
        diff1 += tempD1[i1]- D1[i1] + tempD1[i2]- D1[i2];
        tempD2[i1][i1] *= (1+0.5*xc[i2][ncol]+0.5*xc[i2][ncol])/(1+0.5*xc[i1][ncol]+0.5*xc[i1][ncol]);
        tempD2[i2][i2] *= (1+0.5*xc[i1][ncol]+0.5*xc[i1][ncol])/(1+0.5*xc[i2][ncol]+0.5*xc[i2][ncol]);
        diff2 += tempD2[i1][i1] - D2[i1][i1] + tempD2[i2][i2] - D2[i2][i2];
        swap(tempx[i1][ncol],tempx[i2][ncol]);
    }
    diff = -2.0/nsamp*diff1+1.0/(nsamp*nsamp)*diff2;
    return(diff);
}