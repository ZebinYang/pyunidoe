//
// Created by Administrator on 2019/4/2.
//

#include "doe_optimizer.h"

vector<vector<double> > Optimizer::get_design()
{
    return (best_design);
}

int Optimizer::select_column()
{
    int i, selected_index;
    vector<int> rank(nv, 0);
    vector<double> prob(nv, 0);
    for (i = 0;i<nv; i ++)
    {
        prob[i] = ((double) (rand()-1) / (RAND_MAX))*optimize_columns[i];
    }
    rank = get_sorted_index(nv, prob);
    selected_index = rank.back();
    return (selected_index);
}

typedef pair<double,unsigned int> mypair_double;
bool comparator_double(const mypair_double& l, const mypair_double& r) {return l.first < r.first;}
vector<int> Optimizer::get_sorted_index(int n, vector<double> arr)
{
    int i;
    int size = (int)(arr.size());
    vector<int> indx(size, 0);
    vector<mypair_double> arr_temp(size);
    for (i = 0; i<size; i++)
    {
        arr_temp[i].first = arr[i]; arr_temp[i].second = static_cast<unsigned int>(i);
    }
    sort(arr_temp.begin(), arr_temp.begin() + n, comparator_double);
    for (i = 0; i < n; i++) {
        indx[i] = arr_temp[i].second;
    }
    return(indx);
}

Optimizer::Optimizer(vector<vector<double> > x_init, int nnew_init, int np_init, int nv_init, int nlevel_init,vector<int> optimize_columns_init,
                   int crit_init, int maxiter_init, double hits_ratio_init, bool level_permt_init)
{
    int i, j, k;
    vector<int> idx;
    vector<double> temp;
    vector<int> freq_table(nlevel_init, 0);
    int max_freq, min_freq, valid_nlevel, valid_nvalues, allpairs, allpairs_temp;

    np = np_init;
    nv = nv_init;
    nnew = nnew_init;
    nsamp = np + nnew;
    nlevel = nlevel_init;
    optimize_columns = optimize_columns_init;
    level_permt.assign(nv, level_permt_init);

    //  Initialize some important parameters
    maxcol = 0;
    maxpairs = 0;
    idx.assign(nnew, 0);
    temp.assign(nnew, 0);
    group_num.assign(nv, 0);
    group_size.assign(nv, 0);
    nlevel_pairs.assign(nv, 0);
    nelement_pairs.assign(nv, 0);
    sorted_value_index.assign(static_cast<unsigned long long int>(nsamp), vector<int>(nv, 0));
    sorted_level_index.assign(static_cast<unsigned long long int>(nnew), vector<int>(nv, 0));
    allpairs = nnew * (nnew -1)/2 ;
    for (i=0;i<nv;i++)
    {
        maxcol += 5*optimize_columns[i];
        allpairs_temp = allpairs;
        max_freq = 0; min_freq = nnew;
        valid_nlevel = 0, valid_nvalues = 0;
        for (j=0;j<nnew;j++) temp[j] = x_init[j+np][i];
        idx = get_sorted_index(nnew, temp);
        for (j = 0; j < nnew; j++) {
            sorted_value_index[j][i] = (int) idx[j] + np;
        }
        for (j=0;j<nlevel;j++) freq_table[j] = 0;
        for (j=0;j<nlevel;j++)
        {
            freq_table[j] = count(temp.begin(), temp.end(), j+1);
            allpairs_temp -= freq_table[j]*(freq_table[j]-1)/2;
            if (freq_table[j]>0)
            {
                valid_nlevel +=1;
                valid_nvalues +=freq_table[j];
                for (k=freq_table[j];k>0;k--) sorted_level_index[valid_nvalues-k][i] = valid_nlevel;
            }
            if ((max_freq<freq_table[j]) & (freq_table[j]!=0)) max_freq = freq_table[j];
            if ((min_freq>freq_table[j]) & (freq_table[j]!=0)) min_freq = freq_table[j];
        }
        nelement_pairs[i] = (int) (min(max(0.2 * allpairs_temp, 1.0), 50.0));
        if (min_freq!=max_freq)
        {
            nlevel_pairs[i] = (int) (min(max(0.2 * valid_nlevel * (valid_nlevel - 1) / 2, 1.0), 50.0));
        } else {
            level_permt[i] = false;
            nlevel_pairs[i] = 0;
        }
        if (level_permt[i]) group_num[i] = valid_nlevel; else group_num[i] = nnew;
        group_size[i] = nnew / group_num[i];
        if (maxpairs<nlevel_pairs[i]) maxpairs=nlevel_pairs[i];
        if (maxpairs<nelement_pairs[i]) maxpairs=nelement_pairs[i];
    }
    //  Initialize the exchange index list
    ind1.assign(maxpairs, 0);
    ind2.assign(maxpairs, 0);

    //  Initialize the criteria
    crit = crit_init;
    switch(crit)
    {
        case 1:
            c = new CD2(x_init, nsamp, nv, nlevel);
            break;
        case 2:
            c = new MD2(x_init, nsamp, nv, nlevel);
            break;
        case 3:
            c = new WD2(x_init, nsamp, nv, nlevel);
            break;
        case 4:
            c = new Maximin(x_init, nsamp, nv, nlevel);
            break;
        case 5:
            c = new MC(x_init, nsamp, nv, nlevel);
            break;
        case 6:
            c = new A2(x_init, nsamp, nv, nlevel);
            break;
        default:
            c = new CD2(x_init, nsamp, nv, nlevel);
            break;
    }

    // other optimization settings
    maxiter = maxiter_init;
    hits_ratio = hits_ratio_init;

    best_design = c->get_design();
    surrogate_obj = c->get_surrogate_criteria();
    global_obj = c->get_criteria();
    th0 = 0.005 * surrogate_obj;
}

vector<double> Optimizer::SATA_Optimize()
{
    int i,j, ii1, ii2, ncol;
    vector<int> i1(maxpairs, 0), i2(maxpairs, 0);
    vector<double> critobj_vector;
    double iteration, new_obj, th, cprob, hits;
    vector<double> col_weight;

    th = th0;
    critobj_vector.push_back(global_obj);
    for (iteration = 0; iteration < maxiter; iteration++)
    {
        hits = 0;
        for (i = 0; i < maxcol; i++)
        {
            ncol = select_column();
            new_obj = (level_permt[ncol]) ? level_permutation(ncol) : elementwise_exchange(ncol);
            cprob = 1 - min(1.0, max(0.0, (new_obj - surrogate_obj) * 1.0 / th));
            hits = hits + cprob;
            if (cprob > ((double) (rand()-1) / (RAND_MAX)))
            {
                surrogate_obj = new_obj;
                ii1 = group_size[ncol]*ind1[best_pair];
                ii2 = group_size[ncol]*ind2[best_pair];
                for (j = 0; j < group_size[ncol]; j++)
                {
                    i1[j] = sorted_value_index[ii1 + j][ncol];
                    i2[j] = sorted_value_index[ii2 + j][ncol];
                    swap(sorted_value_index[ii1+j][ncol], sorted_value_index[ii2+j][ncol]);
                }
                c->set_columnwise_exchange(ncol, group_size[ncol], i1, i2);
                col_weight = c->get_criteria_matrix();
                for (j = 0; j < nv; j++) optimize_columns[j] = col_weight[j];
                if (global_obj - c->get_criteria() > EPS) {
                    global_obj = c->get_criteria();
                    best_design = c->get_design();
                }
            }
        }
        if ((hits / maxcol) < hits_ratio) th = th * alpha1; else th = th * alpha2;
        critobj_vector.push_back(global_obj);
    }
    return (critobj_vector);
}

void Optimizer::get_pairs_element(int ncol)
{
    int i = 0, j = 0, tp = 0, allparis;
    allparis = nnew * (nnew - 1) / 2;

    vector<int> ind(allparis, 0);
    for (i=0;i<allparis;i++) ind[i]=i+1;

    random_shuffle(ind.begin(), ind.end());
    for (i = 0, j = 0; (i < allparis) & (j < nelement_pairs[ncol]); i++) {
        tp = (int) (sqrt(ind[i] * 2));
        if ((ind[i] * 2) <= tp * (tp + 1)) {
            ind2[j] = tp;
        } else {
            ind2[j] = tp + 1;
        }
        ind1[j] = ind[i] - (ind2[j] -1) * (ind2[j]) / 2 - 1;
        if (sorted_level_index[ind1[j]][ncol] != sorted_level_index[ind2[j]][ncol]) j++;
    }
}

double Optimizer::elementwise_exchange(int ncol)
{
    int i;
    vector<int> i1(1, 0), i2(1, 0);
    double critobj,critobj_candidate;

    get_pairs_element(ncol);
    critobj_candidate = pow(10,10);
    for (i = 0; i < nelement_pairs[ncol]; i++) {
        i1[0] = sorted_value_index[ind1[i]][ncol];
        i2[0] = sorted_value_index[ind2[i]][ncol];
        critobj = c->get_columnwise_exchange(ncol, group_size[ncol], i1, i2);
        if (critobj_candidate > critobj) {
            critobj_candidate = critobj;
            best_pair = i;
        }
    }
    return (critobj_candidate);
}

void Optimizer::get_pairs_level(int ncol)
{
    int i,j, tp, allparis;
    allparis = (group_num[ncol]) * (group_num[ncol] - 1) / 2;

    vector<int> ind(allparis, 0);
    for (i=0;i<allparis;i++) ind[i]=i+1;
    random_shuffle(ind.begin(), ind.end());

    for (i = 0, j = 0; (i < allparis) & (j < nlevel_pairs[ncol]); i++) {
        tp = (int) (sqrt(ind[i] * 2));
        if ((ind[i] * 2) <= tp * (tp + 1)) {
            ind2[j] = tp;
        } else {
            ind2[j] = tp + 1;
        }
        ind1[j] = ind[i] - (ind2[j] -1) * (ind2[j]) / 2 - 1;
        if (ind1[j] != ind2[j]) j++;
    }
}

double Optimizer::level_permutation(int ncol)
{
    int i, j, ii1, ii2;
    vector<int> i1(nlevel_pairs[ncol], 0), i2(nlevel_pairs[ncol], 0);
    double critobj, critobj_candidate;

    get_pairs_level(ncol);
    critobj_candidate = pow(10,10);
    for (i = 0; i < nlevel_pairs[ncol]; i++)
    {
        ii1 = group_size[ncol] * ind1[i];
        ii2 = group_size[ncol] * ind2[i];
        for (j = 0; j < group_size[ncol]; j++) {
            i1[j] = sorted_value_index[ii1 + j][ncol];
            i2[j] = sorted_value_index[ii2 + j][ncol];
        }
        critobj = c->get_columnwise_exchange(ncol, group_size[ncol], i1, i2);
        if (critobj_candidate>critobj) {
            critobj_candidate = critobj;
            best_pair = i;
        }
    }
    return (critobj_candidate);
}