//
// Created by Administrator on 2019/4/1.
//

#ifndef CPP_VERSION_CRIT_H
#define CPP_VERSION_CRIT_H

#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <numeric>
#include <algorithm>

#define EPS 1.0e-10

using namespace std;

class Criteria {

protected:
    int nv;
    int nsamp;
    int nlevel;

    double obj;
    double surrogate_obj;

    vector<vector<double> > x, invx, tempx;

public:
    Criteria(int,int,int);
    double get_criteria();
    double get_surrogate_criteria();
    double get_columnwise_exchange(int ncol,int ncp, vector<int> idx1,vector<int> idx2);
    void set_columnwise_exchange(int ncol,int ncp, vector<int> idx1,vector<int> idx2);

    virtual void init_design(vector<vector<double> > init) = 0;
    virtual void update_design(vector<vector<double> > init) = 0;
    virtual void evaluate_criteria() = 0;
    virtual vector<vector<double> > get_design() = 0;
    virtual double columnwise_exchange(int ncol,int ncp, vector<int> idx1,vector<int> idx2) = 0;

    virtual double calculate_lower_bound(){return(0.0);} // only available for some criteria.
    virtual vector<double> get_criteria_matrix() {vector<double> col_count(nv, 1); return (col_count);};
};

#endif //CPP_VERSION_CRIT_H