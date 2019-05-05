#include "wrapper.h"

int criteria_selector(char* crit)
{
  int critopt;
  critopt = 1;
  if(!strcmp(crit,"CD2")) critopt=1;
  else if(!strcmp(crit,"MD2")) critopt=2;
  else if(!strcmp(crit,"WD2")) critopt=3;
  else if(!strcmp(crit,"maximin")) critopt=4;
  else if(!strcmp(crit,"MC")) critopt=5;
  else if(!strcmp(crit,"A2")) critopt=6;
  return critopt;
}

double CritEval(vector<vector<double> > x0, int nlevel, char* crit)
{
  Criteria *c;
  double criteria = 0;
  int i, j, nv = (int) x0[0].size(), nsamp= (int) x0.size();
  int critopt = criteria_selector(crit);
  vector<vector<double> > x(nsamp, vector<double>(nv, 0));

  for(i=0;i<nsamp;i++) for(j=0;j<nv;j++) x[i][j] = x0[i][j];
    
  switch(critopt)
  {
  case 1:
    c = new CD2(x, nsamp, nv, nlevel);
    break;
  case 2:
    c = new MD2(x, nsamp, nv, nlevel);
    break;
  case 3:
    c = new WD2(x, nsamp, nv, nlevel);
    break;
  case 4:
    c = new Maximin(x, nsamp, nv, nlevel);
    break;
  case 5:
    c = new MC(x, nsamp, nv, nlevel);
    break;
  case 6:
    c = new A2(x, nsamp, nv, nlevel);
    break;
  default:
    c = new CD2(x, nsamp, nv, nlevel);
  break;
  }
  c->evaluate_criteria();
  criteria = c->get_criteria();
  return(criteria);
}

vector<vector<double> > Generate_init_matrix(char* init_method, int nsamp, int nv, int nlevel, vector<vector<double> > initX)
{
  int i,j;
  vector<double> col;
  vector<vector<double> > return_matrix(nsamp, vector<double>(nv, 0));
    
  if ((!strcmp(init_method,"input")) && (initX.size()>1))
  {
    for(i=0;i<(int) initX[0].size();i++) {
      for(j=0;j<(int) initX.size();j++) {
        return_matrix[j][i] = initX[j][i];
      }
    }
  }
  else if ((!strcmp(init_method,"rand")))
  {
    for(i=1;i<=nsamp;i++) col.push_back((i%nlevel)+1.0);
    for(i=0;i<nv;i++)
    {
      random_shuffle (col.begin(), col.end());
      for(j=0;j<nsamp;j++)  return_matrix[j][i] = col[j];
    }
  }
  return return_matrix;
}

vector<vector<double> > Generate_Aug_matrix(char* init_method, vector< vector < double > > xp, int nnew, int nv,
                                  int nlevel, vector< vector < double > > initX)
{
    int i,j,k, np, nsamp, fill_size, all_fill_size, max_k = 0;
    np = (int) xp.size();
    nsamp = np + nnew;
    vector<double> temp(np,0);

    vector<vector<int> > freq_table;
    vector<vector<double> > return_matrix(nnew, vector<double>(nv, 0));
    freq_table.assign(nlevel, vector<int>(nv, 0));

    if ((!strcmp(init_method,"input")) && (initX.size()>1))
    {
      for(i=0;i<(int) initX[0].size();i++) {
        for(j=0;j<(int) initX.size();j++) {
          return_matrix[j][i] = initX[j][i];
        }
      }
    }
    else if ((!strcmp(init_method,"rand")))
    {
      for (i=0;i<nv;i++)
      {
        for (j=0;j<np;j++) temp[j] = xp[j][i];
        for (j=0;j<nlevel;j++)
        {
          freq_table[j][i] = (int) count(temp.begin(), temp.end(), j+1);
          if (max_k<freq_table[j][i]) max_k = freq_table[j][i];
        }
      }
      max_k = max(max_k, nsamp/nlevel);
      for (i=0;i<nv;i++)
      {
        all_fill_size = 0;
        for (j=0;j<nlevel;j++)
        {
          fill_size = max_k-freq_table[j][i];
          for (k=all_fill_size;k<all_fill_size+fill_size;k++)
          {
            return_matrix[k][i] = j+1;
          }
          all_fill_size += fill_size;
        }
      }
    }
    return return_matrix;
}


List SATA_UD(int nsamp, int nv, int nlevel, char* init_method, vector<vector<double> > initX,
             char* crit, int maxiter, double hits_ratio, bool levelpermt, int rand_seed)
{
  List lst;
  int i,j;
  clock_t start_time;
  vector<double> critobj_vector;
  vector<int> optimize_columns(nv,1);
  vector<vector<double> > final_design;
  double critobj, critobj0, search_time;

  int critopt = criteria_selector(crit);
  vector<vector<double> > Init_matrix, return_matrix(nsamp, vector<double>(nv, 0));
  vector<vector<double> > x(nsamp, vector<double>(nv, 0));

  srand(rand_seed);
  start_time = clock();
  Init_matrix = Generate_init_matrix(init_method,nsamp,nv,nlevel,initX);
  for(i=0;i<nsamp;i++) for(j=0;j<nv;j++) x[i][j] = Init_matrix[i][j];
  Optimizer opt(x, nsamp, 0, nv, nlevel, optimize_columns, critopt, maxiter, hits_ratio, levelpermt);
  critobj_vector = opt.SATA_Optimize();
  final_design = opt.get_design();
  critobj0 = critobj_vector.front();
  critobj = critobj_vector.back();
  for(i=0;i<nsamp; i++) for(j=0;j<nv;j++) return_matrix[i][j] = final_design[i][j];

  search_time = (double)(clock()-start_time)/CLOCKS_PER_SEC;

  lst.Init_Design = Init_matrix;
  lst.Final_Design = return_matrix;
  lst.Init_Obj = critobj0;
  lst.Final_Obj = critobj;
  lst.Time_Second= search_time;
  lst.Criterion_history = critobj_vector;
  return lst;
}

List SATA_AUD(vector<vector<double> > xp,int nnew, int nv, int nlevel, char* init_method, vector<vector<double> > initX,
              char* crit, int maxiter, double hits_ratio, bool levelpermt, int rand_seed)
{  
  List lst;
  int i,j;
  int np= (int) xp.size();
  int nsamp = np+nnew;
  clock_t start_time;
  vector<double> critobj_vector;
  vector<int> optimize_columns(nv,1);
  vector<vector<double> > final_design;
  double critobj, critobj0, search_time;
  int critopt = criteria_selector(crit);
  vector<vector<double> > x(nsamp, vector<double>(nv, 0));
  vector<vector<double> > InputX(nnew, vector<double>(nv, 0));
  vector<vector<double> > Init_matrix(nsamp, vector<double>(nv, 0));
  vector<vector<double> > return_matrix(nsamp, vector<double>(nv, 0));

  srand(rand_seed);
  start_time = clock();
  InputX = Generate_Aug_matrix(init_method,xp,nnew,nv,nlevel,initX);
  for(j=0;j<nv;j++)
  {
    for(i=0;i<np;i++) {x[i][j] =Init_matrix[i][j]= xp[i][j];}
    for(i=0;i<nnew;i++) {x[i+np][j] = Init_matrix[i+np][j]=InputX[i][j];}
  }
  Optimizer opt(x, nnew, np, nv, nlevel, optimize_columns, critopt, maxiter, hits_ratio, levelpermt);
  critobj_vector = opt.SATA_Optimize();
  final_design = opt.get_design();
  critobj0 = critobj_vector[0];
  critobj = critobj_vector.back();
  for(i=0;i<nsamp; i++) for(j=0;j<nv;j++) return_matrix[i][j] = final_design[i][j];
  search_time = (double)(clock()-start_time)/CLOCKS_PER_SEC;

  lst.Init_Design = Init_matrix;
  lst.Final_Design = return_matrix;
  lst.Init_Obj = critobj0;
  lst.Final_Obj = critobj;
  lst.Time_Second= search_time;
  lst.Criterion_history = critobj_vector;
  return lst;
}

List SATA_AUD_COL(vector<vector<double> > xp, int nvnew, int nlevel, char* init_method, vector<vector<double> > initX,
                  char* crit, int maxiter, double hits_ratio, bool levelpermt, int rand_seed)
{
  List lst;
  int i,j;
  int nsamp = (int)xp.size();
  int nvp = (int)xp[0].size();
  int nv = nvnew+nvp;
  clock_t start_time;
  vector<double> critobj_vector;
  vector<int> optimize_columns(nv,1);
  vector<vector<double> > final_design;
  double critobj, critobj0, search_time;
  int critopt = criteria_selector(crit);
  vector<vector<double> > x(nsamp, vector<double>(nv, 0));
  vector<vector<double> > InputX(nsamp, vector<double>(nvnew, 0));
  vector<vector<double> > Init_matrix(nsamp, vector<double>(nv, 0));
  vector<vector<double> > return_matrix(nsamp, vector<double>(nv, 0));

  srand(rand_seed);
  start_time = clock();
  InputX = Generate_init_matrix(init_method,nsamp,nvnew,nlevel,initX);
  for(j=0;j<nvp;j++) optimize_columns[j] = 0;
  for(i=0;i<nsamp;i++)
  {
    for(j=0;j<nvp;j++) x[i][j] =Init_matrix[i][j]= xp[i][j];
    for(j=0;j<nvnew;j++) x[i][j+nvp] =Init_matrix[i][j+nvp]= InputX[i][j];
  }

  Optimizer opt(x, nsamp, 0, nv, nlevel, optimize_columns, critopt, maxiter, hits_ratio, levelpermt);
  critobj_vector = opt.SATA_Optimize();
  final_design = opt.get_design();
  critobj0 = critobj_vector[0];
  critobj = critobj_vector.back();
  for(i=0;i<nsamp; i++) for(j=0;j<nv;j++) return_matrix[i][j] = final_design[i][j];
  search_time = (double)(clock()-start_time)/CLOCKS_PER_SEC;

  lst.Init_Design = Init_matrix;
  lst.Final_Design = return_matrix;
  lst.Init_Obj = critobj0;
  lst.Final_Obj = critobj;
  lst.Time_Second= search_time;
  lst.Criterion_history = critobj_vector;
  return lst;
}