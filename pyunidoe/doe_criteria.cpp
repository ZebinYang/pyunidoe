//
// Created by Administrator on 2019/4/1.
//
#include "doe_criteria.h"

Criteria::Criteria(int nsamp_init,int nv_init,int nlevel_init)
{
    nv = nv_init;
    nsamp = nsamp_init;
    nlevel = nlevel_init;
    x.assign(static_cast<unsigned long long int>(nsamp), vector<double>(nv, 0));
    invx.assign(static_cast<unsigned long long int>(nsamp), vector<double>(nv, 0));
}

double Criteria::get_columnwise_exchange(int ncol,int ncp, vector<int> idx1,vector<int> idx2)
{
    double diff;
    diff = columnwise_exchange(ncol, ncp, idx1, idx2);
    evaluate_criteria();
    return(surrogate_obj+diff);
}

void Criteria::set_columnwise_exchange(int ncol, int ncp, vector<int> idx1,vector<int> idx2)
{
    int i1, i2, j;
    tempx = x;
    for(j=0;j<ncp;j++) {
        i1 = idx1[j];
        i2 = idx2[j];
        swap(tempx[i1][ncol], tempx[i2][ncol]);
    }
    update_design(tempx);
    evaluate_criteria();
}

double Criteria::get_criteria() {
    return (obj);
}

double Criteria::get_surrogate_criteria()
{
    return(surrogate_obj);
}