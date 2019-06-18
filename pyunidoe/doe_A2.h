//
// Created by Administrator on 2019/4/14.
//
#include "doe_criteria.h"

using namespace std;

#ifndef CPP_VERSION_DOE_EXPSQUARE_H
#define CPP_VERSION_DOE_EXPSQUARE_H

class A2: public Criteria{

    vector<vector<double> > CORR, tempCORR;
    double A;
public:
    A2(vector<vector<double> > init,int nsamp_init,int nv_init, int nlevel_init): Criteria(nsamp_init,nv_init, nlevel_init) {
        init_design(init);
        evaluate_criteria();
    }
    vector<vector<double> > get_design();
    void init_design(vector<vector<double> >);
    void update_design(vector<vector<double> >);
    void evaluate_criteria();
    double columnwise_exchange(int ncol,int ncp, vector<int> idx1,vector<int> idx2);
};

#endif //CPP_VERSION_DOE_EXPSQUARE_H
