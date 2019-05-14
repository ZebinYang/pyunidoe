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
                count = 1;
            } else if ( abs(abs(CORR[i][j]) - obj ) < EPS)
            {
                count += 1;
            }
            CORR[j][i] = CORR[i][j];
        }
    }
    surrogate_obj = obj * 10000 * (nsamp*(nsamp-1)/2) + count;
}

double MC::columnwise_exchange(int ncol, int ncp, vector<int> idx1,vector<int> idx2)
{
    int i1,i2,i,j, count = 0;
    double temp_obj,diff;
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
                    count = 1;
                } else if (abs(abs(tempCORR[i][j]) - temp_obj)<EPS) {
                    count += 1;
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
                if ( (abs(tempCORR[i][ncol]) - temp_obj ) > EPS)
                {
                    temp_obj = abs(tempCORR[i][ncol]) ;
                    count = 1;
                } else if (abs(abs(tempCORR[i][ncol]) - temp_obj)<EPS)
                {
                    count += 1;
                }

            } else if (ncol<i)
            {
                tempCORR[ncol][i] += x[i1][ncol] * x[i2][i] +  x[i2][ncol] * x[i1][i] -  x[i1][ncol] * x[i1][i] -  x[i2][ncol] * x[i2][i];
                if ( (abs(tempCORR[ncol][i]) - temp_obj ) > EPS)
                {
                    temp_obj = abs(tempCORR[ncol][i]);
                    count = 1;
                } else if (abs(abs(tempCORR[ncol][i]) - temp_obj)<EPS)
                {
                    count += 1;
                }
            }
        }
        swap(tempx[i1][ncol], tempx[i2][ncol]);
    }
    diff =  temp_obj * 10000 * (nsamp*(nsamp-1)/2) + count - surrogate_obj;
    return(diff);
}

double MC::calculate_lower_bound()
{
    double lb1 = 0, lb2 = 0, k1, k2;
    if (nv>nsamp)
    {
        lb1 = sqrt((nv-nsamp+1.0)/((nsamp-1.0)*(nv-1.0)));
        if (((nv%2 == 0)&(nlevel%2 == 0))|(nlevel%2 == 1))
        {
            k1 = floor(-2*nv*A/(nsamp*(nsamp-1)));
            lb2 = sqrt(-k1*nsamp*(nsamp-1)*(k1+1.0)/(4*nv*(nv-1)*A*A) - (2*k1+1)/(2*(nv-1)*A) +1.0*(nv-nsamp)/(nsamp*(nv-1)));
        } else {
            k2 = floor((-4*nv*A-nsamp*(nsamp-1))/(2*nsamp*(nsamp-1)));
            lb2 = sqrt(-nsamp*(nsamp-1)*(4*k2*k2+8*k2+3)/(16*nv*(nv-1)*A*A)-(k2+1)/((nv-1)*A)+1.0*(nv-nsamp)/(nsamp*(nv-1)));
        }
    }
    cout<<"lower bound:"<<lb1<<" "<<lb2<<endl;
    return(lb2);
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