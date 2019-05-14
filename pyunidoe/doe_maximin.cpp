//
// Created by Administrator on 2019/4/9.
//

#include "doe_maximin.h"

void Maximin::init_design(vector<vector<double> > init)
{
    int i,j;
    if(!init.empty())
    {
        for(i=0;i<nsamp;i++) for(j=0;j<nv;j++)
            {
                x[i][j]=(init[i][j]-1)/(nlevel-1);
//                x[i][j]=(2*init[i][j]-1)/(2*nlevel);
            }
    }
}

void Maximin::update_design(vector<vector<double> > init)
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

vector<vector<double> > Maximin::get_design()
{
    int i,j;
    for(j=0;j<nv;j++)
    {
        for(i=0;i<nsamp;i++)
        {
            invx[i][j]=x[i][j]*(nlevel-1)+1;
//            invx[i][j]=x[i][j]*nlevel+0.5;
        }
    }
    return(invx);
}

void Maximin::evaluate_criteria()
{
    int i,j,k;
    obj = 0;
    DIST.assign(static_cast<unsigned long long int>(nsamp), vector<double>(nsamp, 0));
    for(i=0;i<nsamp;i++)
    {
        for(j=i+1;j<nsamp;j++)
        {
            for (k=0;k<nv;k++)
            {
                DIST[i][j]+=pow(abs(x[i][k]-x[j][k]), t);
            }
            DIST[i][j] = max(DIST[i][j], pow(EPS, t));
            DIST[j][i] = DIST[i][j];
            obj += pow(DIST[j][i], -p/t);
        }
    }
    surrogate_obj = obj = pow(obj, 1.0/p);
}

double Maximin::columnwise_exchange(int ncol,int ncp, vector<int> idx1,vector<int> idx2)
{
    int i1,i2,i,j,k;
    double temp_obj,diff;
    temp_obj = pow(surrogate_obj, p);
    tempx = x; tempDIST = DIST;
    for(k=0;k<ncp;k++)
    {
        i1=idx1[k]; i2=idx2[k];
        if (i1>i2) swap(i1,i2);
        for (i=0;i<nsamp;i++)
        {
            for (j=i+1;j<nsamp;j++)
            {
                if (tempDIST[i][j]==pow(EPS, t)) tempDIST[i][j] = 0;
                if (j==i1)
                {
                    tempDIST[i][j] += pow(abs(x[i][ncol] - x[i2][ncol]), t) -  pow(abs(x[i][ncol] - x[i1][ncol]),t);
                } else if ((j==i2)&(i!=i1))
                {
                    tempDIST[i][j] += pow(abs(x[i][ncol] - x[i1][ncol]), t) -  pow(abs(x[i][ncol] - x[i2][ncol]),t);
                } else if ((i==i1)&(j!=i2))
                {
                    tempDIST[i][j] += pow(abs(x[j][ncol] - x[i2][ncol]), t) -  pow(abs(x[j][ncol] - x[i1][ncol]),t);
                } else if (i==i2)
                {
                    tempDIST[i][j] += pow(abs(x[j][ncol] - x[i1][ncol]), t) -  pow(abs(x[j][ncol] - x[i2][ncol]),t);
                }
                tempDIST[i][j] = max(tempDIST[i][j], EPS*EPS);
                tempDIST[j][i] = tempDIST[i][j];
                temp_obj += pow(tempDIST[j][i], -p/t) - pow(DIST[j][i], -p/t);
            }
        }
        swap(tempx[i1][ncol],tempx[i2][ncol]);
    }
    diff = pow(temp_obj, 1.0/p) - surrogate_obj;
    return(diff);
}