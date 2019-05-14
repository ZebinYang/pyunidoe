//
// Created by Administrator on 2019/4/8.
//

#include "doe_MD2.h"

void MD2::init_design(vector<vector<double> > init)
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

void MD2::update_design(vector<vector<double> > init)
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

vector<vector<double> > MD2::get_design()
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

void MD2::evaluate_criteria()
{
    int i,j,k;
    double part1 = 0, part2 = 0;
    D1.assign(nsamp, 1);
    D2.assign(static_cast<unsigned long long int>(nsamp), vector<double>(nsamp, 1));
    for(i=0;i<nsamp;i++)
    {
        for(k=0;k<nv;k++) {
            D1[i] *= 5.0/3.0 - 0.25 * xc[i][k] - 0.25 * xc[i][k] * xc[i][k];
        }
        for(j=i;j<nsamp;j++) {
            for(k=0;k<nv;k++) {
                D2[i][j] *= 15.0/8.0 - 0.25 * xc[i][k] - 0.25 * xc[j][k] -
                             0.75 * abs(x[i][k] - x[j][k]) + 0.5*pow((x[i][k] - x[j][k]),2);
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
    surrogate_obj = obj = pow(19.0/12.0, nv) - 2.0/nsamp*part1 + 1.0/(nsamp*nsamp)*part2;
}

double MD2::columnwise_exchange(int ncol,int ncp, vector<int> idx1,vector<int> idx2)
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
                    tempD2[i][j] *= (15.0/8.0 - 0.25 * xc[i2][ncol] - 0.25 * xc[i][ncol] - 0.75 * abs(x[i2][ncol] - x[i][ncol]) + 0.5*pow((x[i2][k] - x[i][k]),2))/
                                    (15.0/8.0 - 0.25 * xc[i1][ncol] - 0.25 * xc[i][ncol] - 0.75 * abs(x[i1][ncol] - x[i][ncol]) + 0.5*pow((x[i1][k] - x[i][k]),2));
                } else if ((j == i2) & (i != i1)) {
                    tempD2[i][j] *= (15.0/8.0 - 0.25 * xc[i1][ncol] - 0.25 * xc[i][ncol] - 0.75 * abs(x[i1][ncol] - x[i][ncol]) + 0.5*pow((x[i1][k] - x[i][k]),2))/
                                    (15.0/8.0 - 0.25 * xc[i2][ncol] - 0.25 * xc[i][ncol] - 0.75 * abs(x[i2][ncol] - x[i][ncol]) + 0.5*pow((x[i2][k] - x[i][k]),2));
                } else if ((i == i1) & (j != i2)) {
                    tempD2[i][j] *= (15.0/8.0 - 0.25 * xc[i2][ncol] - 0.25 * xc[j][ncol] - 0.75 * abs(x[i2][ncol] - x[j][ncol]) + 0.5*pow((x[i2][k] - x[j][k]),2))/
                                    (15.0/8.0 - 0.25 * xc[i1][ncol] - 0.25 * xc[j][ncol] - 0.75 * abs(x[i1][ncol] - x[j][ncol]) + 0.5*pow((x[i1][k] - x[j][k]),2));
                } else if (i == i2) {
                    tempD2[i][j] *= (15.0/8.0 - 0.25 * xc[i1][ncol] - 0.25 * xc[j][ncol] - 0.75 * abs(x[i1][ncol] - x[j][ncol]) + 0.5*pow((x[i1][k] - x[j][k]),2))/
                                    (15.0/8.0 - 0.25 * xc[i2][ncol] - 0.25 * xc[j][ncol] - 0.75 * abs(x[i2][ncol] - x[j][ncol]) + 0.5*pow((x[i2][k] - x[j][k]),2));
                }
                tempD2[j][i] = tempD2[i][j];
                diff2 += 2*(tempD2[i][j] - D2[i][j]);
            }
        }
        tempD1[i1] *= (5.0/3.0 - 0.25 * xc[i2][ncol] - 0.25 * xc[i2][ncol] * xc[i2][ncol])/(5.0/3.0 - 0.25 * xc[i1][ncol] - 0.25 * xc[i1][ncol] * xc[i1][ncol]);
        tempD1[i2] *= (5.0/3.0 - 0.25 * xc[i1][ncol] - 0.25 * xc[i1][ncol] * xc[i1][ncol])/(5.0/3.0 - 0.25 * xc[i2][ncol] - 0.25 * xc[i2][ncol] * xc[i2][ncol]);
        diff1 += tempD1[i1]- D1[i1] + tempD1[i2]- D1[i2];
        tempD2[i1][i1] *= (15.0/8.0 - 0.25 * xc[i2][ncol] - 0.25 * xc[i2][ncol])/(15.0/8.0 - 0.25 * xc[i1][ncol] - 0.25 * xc[i1][ncol]);
        tempD2[i2][i2] *= (15.0/8.0 - 0.25 * xc[i1][ncol] - 0.25 * xc[i1][ncol])/(15.0/8.0 - 0.25 * xc[i2][ncol] - 0.25 * xc[i2][ncol]);
        diff2 += tempD2[i1][i1] - D2[i1][i1] + tempD2[i2][i2] - D2[i2][i2];
        swap(tempx[i1][ncol],tempx[i2][ncol]);
    }
    diff = -2.0/nsamp*diff1+1.0/(nsamp*nsamp)*diff2;
    return(diff);
}