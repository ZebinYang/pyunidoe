#ifndef WRAPPER_H__
#define WRAPPER_H__

#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <string.h>
#include "doe_optimizer.h"

using namespace std;

struct List {
    vector<vector<double> > Init_Design;
    vector<vector<double> > Final_Design;
    double Init_Obj;
    double Final_Obj;
    double Time_Second;
    vector<double> Criterion_history;
};

int criteria_selector(char* crit);
vector<vector<double> > Generate_init_matrix(char* init_method, int nsamp, int nv, int nlevel, vector<vector<double> > initX);
vector<vector<double> > Generate_Aug_matrix(char* init_method, vector< vector < double > > xp, int nnew, int nv,
                                  int nlevel, vector< vector < double > > initX);
double CritEval(vector<vector<double> > x0, int nlevel, char* crit);
List SATA_UD(int nsamp, int nv, int nlevel, char* init_method, vector<vector<double> > initX,
             char* crit, int maxiter, double hits_ratio, bool levelpermt, int rand_seed);
List SATA_AUD(vector<vector<double> > xp,int nnew, int nv, int nlevel, char* init_method, vector<vector<double> > initX,
              char* crit, int maxiter, double hits_ratio, bool levelpermt, int rand_seed);
List SATA_AUD_COL(vector<vector<double> > xp, int nvnew, int nlevel, char* init_method, vector<vector<double> > initX,
                  char* crit, int maxiter, double hits_ratio, bool levelpermt, int rand_seed);

#endif /* WRAPPER_H__ */