//
// Created by Administrator on 2019/4/8.
//
#include "doe_criteria.h"

#ifndef CPP_VERSION_WD2_H
#define CPP_VERSION_WD2_H

class WD2: public Criteria {

    vector<vector<double> > xc, D2, tempD2;

public:
    WD2(vector<vector<double> > init,int nsamp_init, int nv_init, int nlevel_init): Criteria(nsamp_init, nv_init, nlevel_init) {
        init_design(init);
        evaluate_criteria();
    }
    vector<vector<double> > get_design();
    void init_design(vector<vector<double> >);
    void update_design(vector<vector<double> >);
    void evaluate_criteria();
    double columnwise_exchange(int ncol,int ncp, vector<int> idx1,vector<int> idx2);
};

#endif //CPP_VERSION_WD2_H
