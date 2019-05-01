import json
import numpy as np
import pandas as pd
import pkg_resources
from seaborn import pairplot
import matplotlib.pyplot as plt
from .pyUniDOE_swig import CritEval, SATA_UD, SATA_AUD, SATA_AUD_COL

__all__ = ["DesignPairsPlot", "DesignEval", "DesignQuery", "GenUD", "GenAUD", "GenAUD_COL", "GenUD_MS", "GenAUD_MS", "GenAUD_COL_MS"]

DATA_PATH = pkg_resources.resource_filename('pyUniDOE', 'data/')

def DesignPairsPlot(x):
    """
    This function draws a pairs plot for checking the design.
    
    Arguments:
    
    -x: an integer matrix object, representing the design matrix.
    """
    pairplot(pd.DataFrame(x))


def DesignQuery(n,s,q,crit="CD2", ShowCrit = True):
    """
    This function takes size of desired design,criterion crit. If the required design exists in database, then return the design, else return NULL.

    Arguments:
    
    -n: an integer object. Run of Experiment.

    -s: an integer object. Factor of Experiment.

    -q: an integer object. Level of Experiment.

    - crit: a character object. Currently, we only support the following two types of criterion:

             "CD2"  --Centered L2 Discrepancy (default) ;

             "MD2"  --Mixture L2  Discrepancy ;
    """
    if (isinstance(n,int)&isinstance(s,int)&isinstance(q,int) == False): 
        return "Wrong types of n,s,q."
    elif ((n%q) != 0): 
        return "n should be multiple of q."
    elif ((s<=1) | (n<=2) | (q <=1)): 
        return "The size of design should be larger than 2*2."

    if (crit=="CD2"): 
        DataX = json.load(open(DATA_PATH+'UD_CD2.json'))
    if (crit=="MD2"): 
        DataX = json.load(open(DATA_PATH+'UD_MD2.json'))

    idx = np.where((np.array(DataX['n']) == n) & (np.array(DataX['s']) == s) & (np.array(DataX['q']) == q))[0]
    if (idx.shape[0] == 0):
        print("No design found.")
        return None
    else:
        D = np.array(DataX['Design'][idx[0]], dtype = int)  
        if(ShowCrit):
            print("CD2 =", DesignEval(D, "CD2"),"MD2 =", DesignEval(D, "MD2"),"Maximin =", DesignEval(D, "maximin"))
    return D


def DesignEval(x, crit="CD2"):
    """
    This function takes matrix X0,q and crit to output the criterion value.

    Arguments:
    - x: an integer matrix object. Representing the design to be evaluated.

    - crit: a character object. Type of criterion to use:

             "CD2"  --Centered L2 Discrepancy (default) ;

             "WD2" -- Wrap-around L2 Discrepancy ;

             "MD2"  --Mixture L2  Discrepancy ;

             "maximin" -- Maximin Discrepancy ;

             "MC" -- Minimize Coherence ;

             "A2" -- Minimize Average Chi-Square.

    Example:
      x = np.array([[1, 2],
               [3, 3],
               [2, 1]])
      crit = "MD2"
      obj = DesignEval(x,crit)
    """
    if (not isinstance(x, np.ndarray)):
        return "The design matrix must be a numpy array."
    elif ((np.min(x)<=0) | (x.dtype != np.int)):
        return "The values in design matrix should be integers: 1,2,3,..."
    
    nlevel = int(np.max(x) - np.min(x) + 1) 
    return CritEval(x.tolist(),nlevel,crit)


def GenUD(n,s,q,init="rand",initX=np.array([[]]),crit="CD2", maxiter=10000,hits_ratio = 0.1,levelpermt=False,vis=False):
    """
    This function takes n,s,q and other arguments to output a list(described below).

    Arguments:
    
    -n: an integer object. Run of Experiment.

    -s: an integer object. Factor of Experiment.

    -q: an integer object. Level of Experiment.

    -init: a string vector object. Initialization method for the design:

              "rand" --Randomly generate initial design (default);

              "input" --User specified.

    -initX: a user-defined numpy integer matrix object. This is the user-defined initial design matrix, and will be used when init="input".

    - crit: a character object. Type of criterion to use:

             "CD2"  --Centered L2 Discrepancy (default) ;

             "WD2" -- Wrap-around L2 Discrepancy ;

             "MD2"  --Mixture L2  Discrepancy ;

             "maximin" -- Maximin Discrepancy ;

             "MC" -- Minimize Coherence ;

             "A2" -- Minimize Average Chi-Square.

    -maxiter: a positive integer  object. Maximum iteration number in outer while loop of SATA algorithm.

    -levelpermt: a boolean object. It controls whether to use level permutation.

    -hits_ratio: a float object. Default value is 0.1, which is the ratio to accept changes of design in inner for loop. 

    -vis: a boolean object. If true, plot the criterion value sequence.

    Examples:
      ## 1
      n=12 #(must be multiples of q)
      s=3
      q=4
      crit = "CD2"#(Centered L2 criteria)
      stat = GenUD(n,s,q,crit=crit,maxiter=100)

      ## 2
      n=10 
      s=3
      q=5
      init = "rand"
      crit = "MD2" #(Mixture L2 criteria)
      vis=TRUE
      stat=GenUD(n,s,q,init=init,crit=crit,maxiter=100,vis=vis)

      ## 3
      #If init="input", algorithm will search for better a better design with same size as initX (balanced design).      
      n=3 
      s=2
      q=3
      initX = np.array([[1, 1],
                  [2, 2],
                  [3, 3]])
      stat = GenUD(n,s,q, init="input", initX = initX, maxiter=100)
    """
   #check the arguments
    if (isinstance(n,int)&isinstance(s,int)&isinstance(q,int) == False): 
        return "Wrong types of n,s,q."
    elif ((n%q) != 0): 
        return "n should be multiple of q."
    elif ((s<=1) | (n<=2) | (q <=1)): 
        return "The size of design should be larger than 2*2."
    elif (init=="input"):
        if ((not isinstance(initX,np.ndarray))):
            return "initX must be a numpy array."
        elif ((n != initX.shape[0]) | (s != initX.shape[1])):
            return "The size of the input design matrix does not match the given n,s."
        elif ((1 != (np.min(initX))) | (q != (np.max(initX))) | (initX.dtype != np.int)):
            return "The values of the input design matrix should be integers within: 1,2,3...,q."
        
    results = SATA_UD(n,s,q,init,initX.tolist(),crit,maxiter,hits_ratio,levelpermt)
    stat = {"initial_design": np.array(results.Init_Design, dtype = int), 
           "final_design": np.array(results.Final_Design, dtype = int), 
          "initial_criterion": results.Init_Obj, 
          "final_criterion": results.Final_Obj, 
          "time_consumed": results.Time_Second, 
          "criterion_history": np.array(results.Criterion_history)}
    
    if vis:
        plt.figure(figsize = (10,6))
        plt.plot(stat['criterion_history'], color = "black", linewidth =1)
        bst_score = round(stat['final_criterion'],5)
        min_index = np.argmin(stat['criterion_history']) 
        plt.axvline(min_index,color="red", linewidth =1)
        plt.axhline(stat['criterion_history'][min_index], linewidth =1 )
        plt.title("Best value = "+ str(bst_score) + " in "+ str(round(stat['time_consumed'],3)) +" sec", fontsize = 20)
        plt.xlabel("Iterations", fontsize = 18)
        plt.ylabel("Criteria", fontsize = 18)
        yticks = np.round(np.linspace(np.min(stat['criterion_history']), np.max(stat['criterion_history']), 4), 4)
        plt.yticks(yticks, fontsize=14)
        plt.xticks(fontsize=14)
    return stat


def GenAUD(xp,n,s,q,init="rand",initX=np.array([[]]),crit="CD2", maxiter=10000,hits_ratio = 0.1,levelpermt=False,vis=False):
    """
    This function takes n,s,q; a unchanged initial design and other arguments to output a list (described below).
    
    Arguments:
    
    -xp: a numpy integer matrix object. Representing the previous existing design matrix.

    -n: an integer object. Run of Experiment (including the previous design xp).

    -s: an integer object. Factor of Experiment.

    -q: an integer object. Level of Experiment.

    -init: a string vector object. Initialization method for the run-augmented design:

              "rand" --Randomly generate initial design (default);

              "input" --User specified.

    -initX: a user-defined numpy integer matrix object. This is the user-defined initial augmentation matrix, and will be used when init="input".

    - crit: a character object. Type of criterion to use:

             "CD2"  --Centered L2 Discrepancy (default) ;

             "WD2" -- Wrap-around L2 Discrepancy ;

             "MD2"  --Mixture L2  Discrepancy ;

             "maximin" -- Maximin Discrepancy ;

             "MC" -- Minimize Coherence ;

             "A2" -- Minimize Average Chi-Square.

    -maxiter: a positive integer  object. Maximum iteration number in outer while loop of SATA algorithm.

    -levelpermt: a boolean object. It controls whether to use level permutation.

    -hits_ratio: a float object. Default value is 0.1, which is the ratio to accept changes of design in inner for loop. 

    -vis: a boolean object. If true, plot the criterion value sequence.

    Examples:
      n=6 
      s=2
      q=3
      xp = np.array([[1, 1],
                [2, 2],
                [3, 3]])
      crit = "CD2"
      res = GenAUD(xp,n,s,q,crit=crit,maxiter=100,vis = True)
    """
    
   #check the arguments
    if (isinstance(n,int)&isinstance(s,int)&isinstance(q,int) == False): 
        return "Wrong types of n,s,q."
    elif ((n%q) != 0): 
        return "n should be multiple of q."
    elif ((s<=1) | (n<=2) | (q <=1)): 
        return "The size of design should be larger than 2*2."
    elif (not isinstance(xp,np.ndarray)):
        return "xp must be a numpy array."
    elif ((n <= xp.shape[0]) | (s != xp.shape[1])):
        return "The size of the existing design matrix xp does not match the given n,s."
    elif ((xp.dtype != np.int) | (1 > np.min(xp)) | (q < np.max(xp))):
        return "The values of the existing design matrix x0 should be integers within: 1,2,3...,q."
    elif (any(np.array([v for i in range(xp.shape[1]) for v in np.unique(xp[:,i], return_counts=True)[1]])>n/q)):
        return "xp does not follow a balanced design, please increase the number of n or remove duplicated elements (per column) in xp."
    elif (init=="input"):
        if ((not isinstance(initX,np.ndarray))):
            return "initX must be a numpy array."
        elif ((n <= initX.shape[0]) | (s != initX.shape[1])):
            return "The size of the input design matrix does not match the given n,s."
        elif ((np.min(initX)<1) | (q < np.max(initX)) | (initX.dtype != np.int)):
            return "The values of the input design matrix should be integers within: 1,2,3...,q."

    nnp = xp.shape[0]
    results = SATA_AUD(xp.tolist(),n-nnp,s,q,init,initX.tolist(),crit,maxiter,hits_ratio,levelpermt)
    stat = {"initial_design": np.array(results.Init_Design, dtype = int), 
           "final_design": np.array(results.Final_Design, dtype = int), 
          "initial_criterion": results.Init_Obj, 
          "final_criterion": results.Final_Obj, 
          "time_consumed": results.Time_Second, 
          "criterion_history": np.array(results.Criterion_history)}
    if vis:
        plt.figure(figsize = (10,6))
        plt.plot(stat['criterion_history'], color = "black", linewidth =1)
        bst_score = round(stat['final_criterion'],5)
        min_index = np.argmin(stat['criterion_history']) 
        plt.axvline(min_index,color="red", linewidth =1)
        plt.axhline(stat['criterion_history'][min_index], linewidth =1 )
        plt.title("Best value = "+ str(bst_score) + " in "+ str(round(stat['time_consumed'],3)) +" sec", fontsize = 20)
        plt.xlabel("Iterations", fontsize = 18)
        plt.ylabel("Criteria", fontsize = 18)
        yticks = np.round(np.linspace(np.min(stat['criterion_history']), np.max(stat['criterion_history']), 4), 4)
        plt.yticks(yticks, fontsize=14)
        plt.xticks(fontsize=14)
    return stat


def GenAUD_COL(xp,n,s,q,init="rand",initX=np.array([[]]),crit="CD2", maxiter=10000,hits_ratio = 0.1,levelpermt=False,vis=False):
        
    """
    This function takes n,s,q; a unchanged initial design and other arguments to output a list (described below). 
    
    Arguments:
    
    -xp: a numpy integer matrix object, representing the previous existing design matrix.

    -n: an integer object. Run of Experiment.

    -s: an integer object. Factor of Experiment (including the number of factors in previous design xp).

    -q: an integer object. Level of Experiment.

    -init: a string vector object. Initialization method for the run-augmented design:

              "rand" --Randomly generate initial design (default);

              "input" --User specified.

    -initX: a user-defined numpy integer matrix object. This is the user-defined initial augmentation matrix, and will be used when init="input".

    - crit: a character object. Type of criterion to use:

             "CD2"  --Centered L2 Discrepancy (default) ;

             "WD2" -- Wrap-around L2 Discrepancy ;

             "MD2"  --Mixture L2  Discrepancy ;

             "maximin" -- Maximin Discrepancy ;

             "MC" -- Minimize Coherence ;

             "A2" -- Minimize Average Chi-Square.

    -maxiter: a positive integer  object. Maximum iteration number in outer while loop of SATA algorithm.

    -levelpermt: a boolean object. It controls whether to use level permutation.

    -hits_ratio: a float object. Default value is 0.1, which is the ratio to accept changes of design in inner for loop. 

    -vis: a boolean object. If true, plot the criterion value sequence.

    Examples:
      n=3 
      s=4
      q=3
      xp = np.array([[1, 1],
                [2, 2],
                [3, 3]])
      crit = "CD2"
      res = GenAUD_COL(xp,n,s,q,crit=crit,maxiter=100,vis = True)
    """
    #check the arguments
    if (isinstance(n,int)&isinstance(s,int)&isinstance(q,int) == False): 
        return "Wrong types of n,s,q."
    elif ((n%q) != 0): 
        return "n should be multiple of q."
    elif ((s<=1) | (n<=2) | (q <=1)): 
        return "The size of design should be larger than 2*2."
    elif (not isinstance(xp,np.ndarray)):
        return "xp must be a numpy array."
    elif ((n != xp.shape[0]) | (s <= xp.shape[1])):
        return "The size of the existing design matrix xp does not match the given n,s."
    elif ((xp.dtype != np.int)| (1 > np.min(xp)) | (q < np.max(xp)) ):
        return "The values of the existing design matrix x0 should be integers within: 1,2,3...,q."
    elif (any(np.array([v for i in range(xp.shape[1]) for v in np.unique(xp[:,i], return_counts=True)[1]])>n/q)):
        return "xp does not follow a balanced design."
    elif (init=="input"):
        if ((not isinstance(initX,np.ndarray))):
            return "initX must be a numpy array."
        elif ((n != initX.shape[0]) | (s <= initX.shape[1])):
            return "The size of the input design matrix does not match the given n,s."
        elif ((1 != np.min(initX)) | (q != np.max(initX)) | (initX.dtype != np.int)):
            return "The values of the input design matrix initX should be integers within: 1,2,3...,q."

    nvp = xp.shape[1]
    results = SATA_AUD_COL(xp.tolist(), s - nvp, q, init, initX.tolist(), crit, maxiter, hits_ratio, levelpermt)
    stat = {"initial_design": np.array(results.Init_Design, dtype = int), 
           "final_design": np.array(results.Final_Design, dtype = int), 
          "initial_criterion": results.Init_Obj, 
          "final_criterion": results.Final_Obj, 
          "time_consumed": results.Time_Second, 
          "criterion_history": np.array(results.Criterion_history)}
    if vis:
        plt.figure(figsize = (10,6))
        plt.plot(stat['criterion_history'], color = "black", linewidth =1)
        bst_score = round(stat['final_criterion'],5)
        min_index = np.argmin(stat['criterion_history']) 
        plt.axvline(min_index,color="red", linewidth =1)
        plt.axhline(stat['criterion_history'][min_index], linewidth =1 )
        plt.title("Best value = "+ str(bst_score) + " in "+ str(round(stat['time_consumed'],3)) +" sec", fontsize = 20)
        plt.xlabel("Iterations", fontsize = 18)
        plt.ylabel("Criteria", fontsize = 18)
        yticks = np.round(np.linspace(np.min(stat['criterion_history']), np.max(stat['criterion_history']), 4), 4)
        plt.yticks(yticks, fontsize=14)
        plt.xticks(fontsize=14)
    return stat


def GenUD_MS(n, s, q, crit="CD2", maxiter=30, nshoot = 5, vis=False):
    """
    This function generates Uniform Design of Experiments using diffrent initializations.
    
    Arguments:
    
    -n: an integer object. Run of Experiment.

    -s: an integer object. Factor of Experiment (including the number of factors in previous design xp).

    -q: an integer object. Level of Experiment.

    - crit: a character object. Type of criterion to use:

             "CD2"  --Centered L2 Discrepancy (default) ;

             "WD2" -- Wrap-around L2 Discrepancy ;

             "MD2"  --Mixture L2  Discrepancy ;

             "maximin" -- Maximin Discrepancy ;

             "MC" -- Minimize Coherence ;

             "A2" -- Minimize Average Chi-Square.

    -maxiter: a positive integer  object. Maximum iteration number in outer while loop of SATA algorithm.

    -nshoot: Total counts to try different initial designs.
    
    -vis: a boolean object. If true, plot the criterion value sequence.
    """
    if (isinstance(n,int)&isinstance(s,int)&isinstance(q,int) == False): 
        return "Wrong types of n,s,q."
    elif ((n%q) != 0): 
        return "n should be multiple of q."
    elif ((s<=1) | (n<=2) | (q <=1)): 
        return "The size of design should be larger than 2*2."
    
    crit_list = []
    shoot_idx = []
    time_list = []
    bestcrit = 1e10
    for i in range(nshoot):
        stat = GenUD(n=n,s=s,q=q, crit=crit, maxiter = maxiter)
        crit_list.append(list(stat["criterion_history"]))
        shoot_idx.append(len(crit_list))
        time_list.append(stat["time_consumed"])
        tmp = DesignEval(stat["final_design"], crit = crit)
        if (tmp < bestcrit):
            bestcrit = tmp
            bestdesign = np.array(stat["final_design"], dtype = int)
    if vis: 
        fig, axes = plt.subplots(int(np.ceil(1.0*nshoot/5)),5, figsize = (15,4*np.ceil(1.0*nshoot/5)), sharex=True, sharey=True)
        plt.setp(axes, xticks=[], yticks=[])
        axbig = fig.add_subplot(111) 
        axbig.set_xticks([])
        axbig.set_yticks([])
        axbig.set_xlabel("Iteration",fontsize = 14)
        axbig.set_ylabel("Criteria Value",fontsize = 14)
        axbig.xaxis.set_label_coords(0.5, -0.1/np.ceil(1.0*nshoot/5))
        axbig.yaxis.set_label_coords(-0.05, 0.5)
        plt.grid(False)
        for i in range(nshoot):
            ax = fig.add_subplot(np.ceil(1.0*nshoot/5),5, i+1)
            ax.plot(crit_list[i]) 
            ax.set_title("Shoot: "+str(i+1))
            ax.title.set_position((0.5,0.9))
            ax.set_xticks(np.linspace(0, maxiter, 3, dtype=int))
            if i % 5!= 0:
                ax.set_yticks([])
            if np.min(crit_list[i])==bestcrit:
                ax.axvline(np.argmin(crit_list[i]), color="red", linewidth =1)
        fig.suptitle("Best value = "+ str(bestcrit) + " in " + str(round(np.sum(time_list),3))+" sec", fontsize = 20,
                 x=0.5, y = 1 - 0.03*np.ceil(1.0*nshoot/5))
        fig.subplots_adjust(wspace = 0)
    return bestdesign


def GenAUD_MS(xp, n, s, q, crit="CD2", maxiter=30, nshoot = 5, vis=False):
    """
    This function generates sequential Uniform Design of Experiments (Augmenting Runs) using diffrent initializations.
    
    Arguments:
    
    -xp: a numpy integer matrix object, representing the previous existing design matrix.

    -n: an integer object. Run of Experiment.

    -s: an integer object. Factor of Experiment (including the number of factors in previous design xp).

    -q: an integer object. Level of Experiment.

    - crit: a character object. Type of criterion to use:

             "CD2"  --Centered L2 Discrepancy (default) ;

             "WD2" -- Wrap-around L2 Discrepancy ;

             "MD2"  --Mixture L2  Discrepancy ;

             "maximin" -- Maximin Discrepancy ;

             "MC" -- Minimize Coherence ;

             "A2" -- Minimize Average Chi-Square.

    -maxiter: a positive integer  object. Maximum iteration number in outer while loop of SATA algorithm.

    -nshoot: Total counts to try different initial designs.
    
    -vis: a boolean object. If true, plot the criterion value sequence.
    """
    if (isinstance(n,int)&isinstance(s,int)&isinstance(q,int) == False): 
        return "Wrong types of n,s,q."
    elif ((n%q) != 0): 
        return "n should be multiple of q."
    elif ((s<=1) | (n<=2) | (q <=1)): 
        return "The size of design should be larger than 2*2."
    elif (not isinstance(xp,np.ndarray)):
        return "xp must be a numpy array."
    elif ((n <= xp.shape[0]) | (s != xp.shape[1])):
        return "The size of the existing design matrix xp does not match the given n,s."
    elif ((xp.dtype != np.int) | (1 > np.min(xp)) | (q < np.max(xp))):
        return "The values of the existing design matrix x0 should be integers within: 1,2,3...,q."
    elif (any(np.array([v for i in range(xp.shape[1]) for v in np.unique(xp[:,i], return_counts=True)[1]])>n/q)):
        return "xp does not follow a balanced design, please increase the number of n or remove duplicated elements (per column) in xp."
    
    crit_list = []
    shoot_idx = []
    time_list = []
    bestcrit = 1e10
    for i in range(nshoot):
        stat = GenAUD(xp=xp, n=n,s=s,q=q, crit=crit, maxiter = maxiter)
        crit_list.append(list(stat["criterion_history"]))
        shoot_idx.append(len(crit_list))
        time_list.append(stat["time_consumed"])
        tmp = DesignEval(stat["final_design"], crit = crit)
        if (tmp < bestcrit):
            bestcrit = tmp
            bestdesign = np.array(stat["final_design"], dtype = int)
    if vis: 
        fig, axes = plt.subplots(int(np.ceil(1.0*nshoot/5)),5, figsize = (15,4*np.ceil(1.0*nshoot/5)), sharex=True, sharey=True)
        plt.setp(axes, xticks=[], yticks=[])
        axbig = fig.add_subplot(111) 
        axbig.set_xticks([])
        axbig.set_yticks([])
        axbig.set_xlabel("Iteration",fontsize = 14)
        axbig.set_ylabel("Criteria Value",fontsize = 14)
        axbig.xaxis.set_label_coords(0.5, -0.1/np.ceil(1.0*nshoot/5))
        axbig.yaxis.set_label_coords(-0.05, 0.5)
        plt.grid(False)
        for i in range(nshoot):
            ax = fig.add_subplot(np.ceil(1.0*nshoot/5),5, i+1)
            ax.plot(crit_list[i]) 
            ax.set_title("Shoot: "+str(i+1))
            ax.title.set_position((0.5,0.9))
            ax.set_xticks(np.linspace(0, maxiter, 3, dtype=int))
            if i % 5!= 0:
                ax.set_yticks([])
            if np.min(crit_list[i])==bestcrit:
                ax.axvline(np.argmin(crit_list[i]), color="red", linewidth =1)
        fig.suptitle("Best value = "+ str(bestcrit) + " in " + str(round(np.sum(time_list),3))+" sec", fontsize = 20,
                 x=0.5, y = 1 - 0.03*np.ceil(1.0*nshoot/5))
        fig.subplots_adjust(wspace = 0)
    return bestdesign


def GenAUD_COL_MS(xp, n, s, q, crit="CD2", maxiter=30, nshoot = 5, vis=False):
    """
    This function generates sequential Uniform Design of Experiments (Augmenting Factors) using diffrent initializations.
    
    Arguments:
    
    -xp: a numpy integer matrix object, representing the previous existing design matrix.

    -n: an integer object. Run of Experiment.

    -s: an integer object. Factor of Experiment (including the number of factors in previous design xp).

    -q: an integer object. Level of Experiment.

    - crit: a character object. Type of criterion to use:

             "CD2"  --Centered L2 Discrepancy (default) ;

             "WD2" -- Wrap-around L2 Discrepancy ;

             "MD2"  --Mixture L2  Discrepancy ;

             "maximin" -- Maximin Discrepancy ;

             "MC" -- Minimize Coherence ;

             "A2" -- Minimize Average Chi-Square.

    -maxiter: a positive integer  object. Maximum iteration number in outer while loop of SATA algorithm.

    -nshoot: Total counts to try different initial designs.
    
    -vis: a boolean object. If true, plot the criterion value sequence.
    """
    
    if (isinstance(n,int)&isinstance(s,int)&isinstance(q,int) == False): 
        return "Wrong types of n,s,q."
    elif ((n%q) != 0): 
        return "n should be multiple of q."
    elif ((s<=1) | (n<=2) | (q <=1)): 
        return "The size of design should be larger than 2*2."
    elif (not isinstance(xp,np.ndarray)):
        return "xp must be a numpy array."
    elif ((n != xp.shape[0]) | (s <= xp.shape[1])):
        return "The size of the existing design matrix xp does not match the given n,s."
    elif ((xp.dtype != np.int)| (1 > np.min(xp)) | (q < np.max(xp)) ):
        return "The values of the existing design matrix x0 should be integers within: 1,2,3...,q."
    elif (any(np.array([v for i in range(xp.shape[1]) for v in np.unique(xp[:,i], return_counts=True)[1]])>n/q)):
        return "xp does not follow a balanced design."

    crit_list = []
    shoot_idx = []
    time_list = []
    bestcrit = 1e10
    for i in range(nshoot):
        stat = GenAUD_COL(xp=xp, n=n,s=s,q=q, crit=crit, maxiter = maxiter)
        crit_list.append(list(stat["criterion_history"]))
        shoot_idx.append(len(crit_list))
        time_list.append(stat["time_consumed"])
        tmp = DesignEval(stat["final_design"], crit = crit)
        if (tmp < bestcrit):
            bestcrit = tmp
            bestdesign = np.array(stat["final_design"], dtype = int)
    if vis: 
        fig, axes = plt.subplots(int(np.ceil(1.0*nshoot/5)),5, figsize = (15,4*np.ceil(1.0*nshoot/5)), sharex=True, sharey=True)
        plt.setp(axes, xticks=[], yticks=[])
        axbig = fig.add_subplot(111) 
        axbig.set_xticks([])
        axbig.set_yticks([])
        axbig.set_xlabel("Iteration",fontsize = 14)
        axbig.set_ylabel("Criteria Value",fontsize = 14)
        axbig.xaxis.set_label_coords(0.5, -0.1/np.ceil(1.0*nshoot/5))
        axbig.yaxis.set_label_coords(-0.05, 0.5)
        plt.grid(False)
        for i in range(nshoot):
            ax = fig.add_subplot(np.ceil(1.0*nshoot/5),5, i+1)
            ax.plot(crit_list[i]) 
            ax.set_title("Shoot: "+str(i+1))
            ax.title.set_position((0.5,0.9))
            ax.set_xticks(np.linspace(0, maxiter, 3, dtype=int))
            if i % 5!= 0:
                ax.set_yticks([])
            if np.min(crit_list[i])==bestcrit:
                ax.axvline(np.argmin(crit_list[i]), color="red", linewidth =1)
        fig.suptitle("Best value = "+ str(bestcrit) + " in " + str(round(np.sum(time_list),3))+" sec", fontsize = 20,
                 x=0.5, y = 1 - 0.03*np.ceil(1.0*nshoot/5))
        fig.subplots_adjust(wspace = 0)
    return bestdesign