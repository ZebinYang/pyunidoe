//
// Created by Administrator on 2019/4/2.
//

#include <iostream>
#include <vector>// std::vector
#include <cmath>
#include <numeric>
#include <algorithm>

#include "doe_criteria.h"
#include "doe_CD2.h"
#include "doe_MD2.h"
#include "doe_WD2.h"
#include "doe_maximin.h"
#include "doe_MC.h"
#include "doe_A2.h"

using namespace std;

#ifndef CPP_VERSION_OPTIMIZE_H
#define CPP_VERSION_OPTIMIZE_H

class Optimizer {
    //  The basic problem setting
    int nsamp,nnew,np,nv,nlevel;

    // The optimization setting
    int crit,maxiter,maxpairs,maxcol;
    double hits_ratio;
    double th0,surrogate_obj,global_obj;
    vector<bool> level_permt;
    double alpha1, alpha2;

    // variables used during optimization
    int best_pair; // the best pair obtained
    vector<int> group_num, group_size;
    vector<int> ind1, ind2;
    vector<int> optimize_columns;
    vector<int> nlevel_pairs, nelement_pairs;
    vector<vector<int> > sorted_value_index; //sorted value index for each column
    vector<vector<int> > sorted_level_index; //sorted level index in the new generated design

    // Output Design
    vector<vector<double> > best_design;

public:
    // Criteria class
    Criteria *c;

public:
    Optimizer(vector<vector<double> > x_init, int nnew_init, int np_init, int nv_init, int nlevel_init,vector<int> optimize_columns_init,
              int crit_init, int maxiter_init, double hits_ratio_init, bool level_permt_init);
    vector<vector<double> > get_design();

    int select_column();
    void get_pairs_element(int ncol);
    void get_pairs_level(int ncol);
    double level_permutation(int ncol);
    double elementwise_exchange(int ncol);
    static vector<int> get_sorted_index(int n, vector<double> arr);
    vector<double> SATA_Optimize();
};
#endif //CPP_VERSION_SECH_H
