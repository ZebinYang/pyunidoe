//
// Created by Administrator on 2019/4/11.
//

#include "doe_A2.h"

void A2::init_design(vector<vector<double> > init)
{
    int i,j;
    power = pow((double) 10, ceil(log10((double) nlevel)));
    NPairs.assign(nv, vector<double>(nv, 0));
    if(!init.empty())
    {
        for(i=0;i<nsamp;i++) for(j=0;j<nv;j++)
            {
                x[i][j]=init[i][j];
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
            invx[i][j]=x[i][j];
    }
    return(invx);
}

void A2::evaluate_criteria()
{
    int i,j,k,m;
    vector<double> temp(nsamp,0);
    vector<double> freq_table(nsamp,0);

    obj = 0;
    NPairs.assign(static_cast<unsigned long long int>(nv), vector<double>(nv, 0));
    for(i=0;i<nv;i++)
    {
        for(j=i+1;j<nv;j++)
        {
            for (k=0;k<nsamp;k++)
            {
                temp[k] = power * x[k][i] + x[k][j];
            }
            sort(temp.begin(), temp.end());

            freq_table[0] = 1;
            for (k=1,m=0;k<nsamp;k++)
            {
                freq_table[k] = 0;
                if ((temp[k]-temp[k-1])>EPS)
                {
                    NPairs[i][j]+= pow((double)(freq_table[m]-1.0*nsamp/pow((double)nlevel,2)),2)/(1.0*nsamp/pow((double)nlevel,2));
                    m++;
                }
                freq_table[m] += 1;
            }
            NPairs[i][j]+= pow((double)(freq_table[m]-1.0*nsamp/pow((double)nlevel,2)),2)/(1.0*nsamp/pow((double) nlevel,2));
            obj += NPairs[i][j];
        }
    }
    obj /= (nv*(nv-1)/2);
    surrogate_obj = obj;
}

double A2::columnwise_exchange(int ncol, int ncp, vector<int> idx1,vector<int> idx2)
{
    int i1,i2,i,j,k,m;
    double temp_obj,diff;
    vector<double> temp(nsamp,0);
    vector<double> freq_table(nsamp,0);
    temp_obj = surrogate_obj;
    tempx = x; tempNPairs = NPairs;
    for(j=0;j<ncp;j++)
    {
        i1=idx1[j]; i2=idx2[j];
        swap(tempx[i1][ncol], tempx[i2][ncol]);
        for (i=0;i<nv;i++)
        {
            freq_table[0] = 1;
            for (k=1;k<nsamp;k++) freq_table[k] = 0;
            if (ncol>i) {
                tempNPairs[i][ncol] = 0;
                for (k=0;k<nsamp;k++)
                {
                    temp[k] = power * tempx[k][i] + tempx[k][ncol];
                }
                sort(temp.begin(), temp.end());
                for (k=1,m=0;k<nsamp;k++)
                {
                    if ((temp[k]-temp[k-1])>EPS)
                    {
                        tempNPairs[i][ncol]+= pow((double)(freq_table[m]-1.0*nsamp/pow((double)nlevel,2)),2)/(1.0*nsamp/pow((double)nlevel,2));
                        m++;
                    }
                    freq_table[m] += 1;
                }
                tempNPairs[i][ncol]+= pow((double)(freq_table[m]-1.0*nsamp/pow((double)nlevel,2)),2)/(1.0*nsamp/pow((double)nlevel,2));
                temp_obj += tempNPairs[i][ncol]/(nv*(nv-1)/2) - NPairs[i][ncol]/(nv*(nv-1)/2);
            } else if (ncol<i)
            {
                tempNPairs[ncol][i] = 0;
                for (k=0;k<nsamp;k++)
                {
                    temp[k] = power * tempx[k][i] + tempx[k][ncol];
                }
                sort(temp.begin(), temp.end());
                for (k=1,m=0;k<nsamp;k++)
                {
                    if ((temp[k]-temp[k-1])>EPS)
                    {
                        tempNPairs[ncol][i]+= pow((double)(freq_table[m]-1.0*nsamp/pow((double)nlevel,2)),2)/(1.0*nsamp/pow((double)nlevel,2));
                        m++;
                    }
                    freq_table[m] += 1;
                }
                tempNPairs[ncol][i]+= pow((double)(freq_table[m]-1.0*nsamp/pow((double)nlevel,2)),2)/(1.0*nsamp/pow((double)nlevel,2));
                temp_obj += tempNPairs[ncol][i]/(nv*(nv-1)/2) - NPairs[ncol][i]/(nv*(nv-1)/2);
            }
        }
    }
    diff = temp_obj - surrogate_obj;
    return(diff);
}
