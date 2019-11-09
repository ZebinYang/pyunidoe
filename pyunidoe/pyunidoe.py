import json
import numpy as np
import pandas as pd
import pkg_resources
from seaborn import pairplot
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from .pyunidoe_swig import CritEval, SATA_UD, SATA_AUD, SATA_AUD_COL

__all__ = ["design_pairs_plot",
           "design_update",
           "design_eval",
           "design_query",
           "gen_ud",
           "gen_aud",
           "gen_aud_col",
           "gen_ud_ms",
           "gen_aud_ms",
           "gen_aud_col_ms"]

DATA_PATH = pkg_resources.resource_filename('pyunidoe', 'data/')


def design_pairs_plot(x):
    """
    This function draws a pairs plot for checking the design.

    Parameters
    ----------
    :type  x: an integer numpy matrix
    :param x: representing the design matrix

    """

    pairplot(pd.DataFrame(x))


def design_eval(x, crit="CD2"):
    """
    This function takes matrix X0,q and crit to output the criterion value.

    Parameters
    ----------
    :type  x: an integer numpy matrix
    :param x: representing the design matrix:

    :type  crit: a character object, default="CD2"
    :param crit: criterion to be evaluated:

             "CD2" -- Centered L2 Discrepancy;

             "WD2" -- Wrap-around L2 Discrepancy;

             "MD2" -- Mixture L2 Discrepancy;

             "maximin" -- Maximin Discrepancy;

             "MC" -- Minimum Coherence;

             "A2" -- Mean Squared Correlation.

    Example
    -------
    >>> x = np.array([[1, 2],
    >>>          [3, 3],
    >>>          [2, 1]])
    >>> crit = "MD2"
    >>> obj = design_eval(x,crit)

    """

    if (not isinstance(x, np.ndarray)):
        raise ValueError("The design matrix must be a numpy array.")
    elif ((np.min(x) <= 0) | (x.dtype != np.int)):
        raise ValueError("The values in design matrix should be integers: 1,2,3,...")

    nlevel = int(round(np.max(x) - np.min(x) + 1))
    return CritEval(x.tolist(), nlevel, crit)


def design_query(n, s, q, crit="CD2", show_crit=True):
    """
    This function takes size of desired design,criterion crit. If the required design exists in database, then return the design, else return NULL.

    Parameters
    ----------
    :type  n: an integer object
    :param n: run of experiments

    :type  s: an integer object
    :param s: number of experimental factors

    :type  q: an integer object
    :param q: number of experimental levels for each factor

    :type  crit: a character object, default="CD2"
    :param crit: criterion of the query:

             "CD2": Centered L2 Discrepancy;

             "MD2": Mixture L2 Discrepancy.

    :type  show_crit: boolean
    :param show_crit: choose to print the criteria value

    """

    if ((isinstance(n, int) & isinstance(s, int) & isinstance(q, int)) is False):
        raise ValueError("Wrong types of n,s,q.")
    elif ((n % q) != 0):
        raise ValueError("n should be multiple of q.")
    elif (n > q**s):
        raise ValueError("n should not be greater than q^s.")
    elif ((s < 1) | (n < 2) | (q < 2)):
        raise ValueError("Invalid design table.")

    if (crit == "CD2"):
        db = json.load(open(DATA_PATH + 'ud_cd2.json'))
    if (crit == "MD2"):
        db = json.load(open(DATA_PATH + 'ud_md2.json'))
    else:
        raise ValueError("Invalid criterion.")

    if (str(n) + "_" + str(s) + "_" + str(q) in db.keys()):
        design_table = np.array(np.round(db[str(n) + "_" + str(s) + "_" + str(q)]), dtype=int)
        if(show_crit):
            print("CD2 = ", design_eval(design_table, "CD2"), "MD2 = ",
                design_eval(design_table, "MD2"), "Maximin = ", design_eval(design_table, "maximin"))
    else:
        design_table = None

    return design_table


def design_update(n, s, q, x, crit="CD2"):
    """
    This function takes size of desired design,criterion crit. If the required design exists in database, then return the design, else return NULL.

    Parameters
    ----------
    :type  n: an integer object
    :param n: run of experiments

    :type  s: an integer object
    :param s: number of experimental factors

    :type  q: an integer object
    :param q: number of experimental levels for each factor

    :type  x: an np array
    :param x: the design table that will be saved in Database.

    :type  crit: a character object, default="CD2"
    :param crit: criterion of the query:

             "CD2": Centered L2 Discrepancy;

             "MD2": Mixture L2 Discrepancy.

    :type  show_crit: boolean
    :param show_crit: choose to print the criteria value

    """

    if ((isinstance(n, int) & isinstance(s, int) & isinstance(q, int)) is False):
        raise ValueError("Wrong types of n,s,q.")
    elif ((n % q) != 0):
        raise ValueError("n should be multiple of q.")
    elif (n > q**s):
        raise ValueError("n should not be greater than q^s.")
    elif ((s < 1) | (n < 2) | (q < 2)):
        raise ValueError("Invalid design table.")
    if (not isinstance(x, np.ndarray)):
        raise ValueError("The design matrix must be a numpy array.")
    elif ((np.min(x) <= 0) | (x.dtype != np.int)):
        raise ValueError("The values in design matrix should be integers: 1,2,3,...")

    success_flag = False
    if (crit == "CD2"):
        db = json.load(open(DATA_PATH + 'ud_cd2.json'))
        if (str(n) + "_" + str(s) + "_" + str(q) in db.keys()):
            current_obj = design_eval(np.array(db[str(n) + "_" + str(s) + "_" + str(q)]), crit="CD2")
            new_obj = design_eval(x, crit="CD2")
            if new_obj <= current_obj:
                db.update({str(n) + "_" + str(s) + "_" + str(q): x.tolist()})
                with open(DATA_PATH + 'ud_cd2.json', 'w') as fp:
                    json.dump(db, fp)
                success_flag = True
        else:
            db.update({str(n) + "_" + str(s) + "_" + str(q): x.tolist()})
            with open(DATA_PATH + 'ud_cd2.json', 'w') as fp:
                json.dump(db, fp)
            success_flag = True

    if (crit == "MD2"):
        db = json.load(open(DATA_PATH + 'ud_md2.json'))
        if (str(n) + "_" + str(s) + "_" + str(q) in db.keys()):
            current_obj = design_eval(np.array(db[str(n) + "_" + str(s) + "_" + str(q)]), crit="MD2")
            new_obj = design_eval(x, crit="CD2")
            if new_obj <= current_obj:
                db.update({str(n) + "_" + str(s) + "_" + str(q): x.tolist()})
                with open(DATA_PATH + 'ud_md2.json', 'w') as fp:
                    json.dump(db, fp)
                success_flag = True
        else:
            db.update({str(n) + "_" + str(s) + "_" + str(q): x.tolist()})
            with open(DATA_PATH + 'ud_md2.json', 'w') as fp:
                json.dump(db, fp)
            success_flag = True
    return success_flag


def gen_ud(n, s, q, init="rand", initX=np.array([[]]), crit="CD2", maxiter=100, hits_ratio=0.1, levelpermt=False, rand_seed=0, vis=False):
    """
    This function takes n,s,q and other arguments to output a list(described below).

    Parameters
    ----------
    :type  n: an integer object
    :param n: run of experiments

    :type  s: an integer object
    :param s: number of experimental factors

    :type  q: an integer object
    :param q: number of experimental levels for each factor

    :type  crit: a character object, default="CD2"
    :param crit: criterion to be optimized:

             "CD2" -- Centered L2 Discrepancy;

             "WD2" -- Wrap-around L2 Discrepancy;

             "MD2" -- Mixture L2 Discrepancy;

             "maximin" -- Maximin Discrepancy;

             "MC" -- Minimum Coherence;

             "A2" -- Mean Squared Correlation.

    :type  init: a string vector object, default="rand"
    :param init: initialization method for the design:

              "rand": randomly generate initial design;

              "input": user specified.

    :type  initX: a user-defined numpy integer matrix object, default=np.array([[]])
    :param initX: This is the user-defined initial design matrix, and will be used when init="input"

    :type  maxiter: a positive integer object, default=100
    :param maxiter: maximum iteration number in outer while loop of SATA algorithm.

    :type  levelpermt: a boolean object, default=False
    :param levelpermt: it controls whether to use level permutation

    :type  hits_ratio: a float object, default=0.1
    :param hits_ratio: Default value is 0.1, which is the ratio to accept changes of design in inner for loop

    :type  rand_seed: an integer object, default=0
    :param rand_seed: random seed

    :type vis: a boolean object, default=False
    :param vis: if true, plot the criterion value sequence

    Examples
    ---------
    >>> ## 1
    >>> n=12 #(must be multiples of q)
    >>> s=3
    >>> q=4
    >>> crit = "CD2"#(Centered L2 criteria)
    >>> stat = gen_ud(n,s,q,crit=crit,maxiter=100)

    ## 2
    >>>  n=10
    >>>  s=3
    >>>  q=5
    >>>  init = "rand"
    >>>  crit = "MD2" #(Mixture L2 criteria)
    >>>  vis=TRUE
    >>>  stat = gen_ud(n,s,q,init=init,crit=crit,maxiter=100,vis=vis)

    ## 3
    >>>  #If init="input", algorithm will search for better a better design with same size as initX (balanced design).
    >>>  n=3
    >>>  s=2
    >>>  q=3
    >>>  initX = np.array([[1, 1],
    >>>              [2, 2],
    >>>              [3, 3]])
    >>>  stat = gen_ud(n,s,q, init="input", initX = initX, maxiter=100)

    """

    # check the arguments
    if ((isinstance(n, int) & isinstance(s, int) & isinstance(q, int)) is False):
        raise ValueError("Wrong types of n,s,q.")
    elif ((n % q) != 0):
        raise ValueError("n should be multiple of q.")
    elif (n > q**s):
        raise ValueError("n should not be greater than q^s.")
    elif ((s < 1) | (n < 2) | (q < 2)):
        raise ValueError("Invalid design table.")
    elif (init == "input"):
        if ((not isinstance(initX, np.ndarray))):
            raise ValueError("initX must be a numpy array.")
        elif ((n != initX.shape[0]) | (s != initX.shape[1])):
            raise ValueError("The size of the input design matrix does not match the given n,s.")
        elif ((1 != (np.min(initX))) | (q != (np.max(initX))) | (initX.dtype != np.int)):
            raise ValueError("The values of the input design matrix should be integers within: 1,2,3...,q.")
        elif (any(np.array([v for i in range(initX.shape[1]) for v in np.unique(initX[:, i], return_counts=True)[1]]) > n / q)):
            raise ValueError("initX does not follow a balanced design.")

    results = SATA_UD(n, s, q, init, initX.tolist(), crit, maxiter, hits_ratio, levelpermt, rand_seed)
    stat = {"initial_design": np.array(np.round(results.Init_Design), dtype=int),
            "final_design": np.array(np.round(results.Final_Design), dtype=int),
            "initial_criterion": results.Init_Obj,
            "final_criterion": results.Final_Obj,
            "time_consumed": results.Time_Second,
            "criterion_history": np.array(results.Criterion_history)}

    if vis:
        plt.figure(figsize=(10, 6))
        plt.plot(stat['criterion_history'], color="black", linewidth=1)
        bst_score = round(stat['final_criterion'], 5)
        min_index = np.argmin(stat['criterion_history'])
        plt.axvline(min_index, color="red", linewidth=1)
        plt.axhline(stat['criterion_history'][min_index], linewidth=1)
        plt.title("Best value = " + str(bst_score) + " in " + str(round(stat['time_consumed'], 3)) + " sec", fontsize=20)
        plt.xlabel("Iterations", fontsize=18)
        plt.ylabel("Criteria", fontsize=18)
        yticks = np.round(np.linspace(np.min(stat['criterion_history']), np.max(stat['criterion_history']), 4), 4)
        plt.yticks(yticks, fontsize=14)
        plt.xticks(fontsize=14)
    return stat


def gen_aud(xp, n, s, q, init="rand", initX=np.array([[]]), crit="CD2", maxiter=100, hits_ratio=0.1, levelpermt=False, rand_seed=0, vis=False):
    """
    This function takes n,s,q; a unchanged initial design and other arguments to output a list (described below).

    Parameters
    ----------
    :type  xp: a numpy integer matrix object
    :param xp: representing the previous existing design matrix

    :type  n: an integer object
    :param n: run of experiments, including the previous design in xp

    :type  s: an integer object
    :param s: number of experimental factors

    :type  q: an integer object
    :param q: number of experimental levels for each factor

    :type  crit: a character object, default="CD2"
    :param crit: criterion to be optimized:

             "CD2" -- Centered L2 Discrepancy;

             "WD2" -- Wrap-around L2 Discrepancy;

             "MD2" -- Mixture L2 Discrepancy;

             "maximin" -- Maximin Discrepancy;

             "MC" -- Minimum Coherence;

             "A2" -- Mean Squared Correlation.

    :type  init: a string vector object, default="rand"
    :param init: initialization method for the run-augmented design:

              "rand": randomly generate initial design;

              "input": user specified.

    :type  initX: a user-defined numpy integer matrix object, default=np.array([[]])
    :param initX: This is the user-defined initial design matrix, and will be used when init="input"

    :type  maxiter: a positive integer object, default=100
    :param maxiter: maximum iteration number in outer while loop of SATA algorithm.

    :type  levelpermt: a boolean object, default=False
    :param levelpermt: it controls whether to use level permutation

    :type  hits_ratio: a float object, default=0.1
    :param hits_ratio: Default value is 0.1, which is the ratio to accept changes of design in inner for loop

    :type  rand_seed: an integer object, default=0
    :param rand_seed: random seed

    :type vis: a boolean object, default=False
    :param vis: if true, plot the criterion value sequence

    Example
    -------
    >>> n=6
    >>> s=2
    >>> q=3
    >>> xp = np.array([[1, 1],
    >>>           [2, 2],
    >>>           [3, 3]])
    >>> crit = "CD2"
    >>> res = gen_aud(xp,n,s,q,crit=crit,maxiter=100,vis = True)

    """

    # check the arguments
    if ((isinstance(n, int) & isinstance(s, int) & isinstance(q, int)) is False):
        raise ValueError("Wrong types of n,s,q.")
    elif ((n % q) != 0):
        raise ValueError("n should be multiple of q.")
    elif (n > q**s):
        raise ValueError("n should not be greater than q^s.")
    elif ((s < 1) | (n < 2) | (q < 2)):
        raise ValueError("Invalid design table.")
    elif (not isinstance(xp, np.ndarray)):
        raise ValueError("xp must be a numpy array.")
    elif ((n <= xp.shape[0]) | (s != xp.shape[1])):
        raise ValueError("The size of the existing design matrix xp does not match the given n,s.")
    elif ((xp.dtype != np.int) | (1 > np.min(xp)) | (q < np.max(xp))):
        raise ValueError("The values of the existing design matrix x0 should be integers within: 1,2,3...,q.")
    elif (any(np.array([v for i in range(xp.shape[1]) for v in np.unique(xp[:, i], return_counts=True)[1]]) > n / q)):
        raise ValueError("xp does not follow a balanced design, please increase the number of n or remove duplicated elements (per column) in xp.")
    elif (init == "input"):
        if ((not isinstance(initX, np.ndarray))):
            raise ValueError("initX must be a numpy array.")
        elif ((n <= initX.shape[0]) | (s != initX.shape[1])):
            raise ValueError("The size of the input design matrix does not match the given n,s.")
        elif ((np.min(initX) < 1) | (q < np.max(initX)) | (initX.dtype != np.int)):
            raise ValueError("The values of the input design matrix should be integers within: 1,2,3...,q.")
        elif (any(np.array([v for i in range(initX.shape[1]) for v in np.unique(initX[:, i], return_counts=True)[1]]) > n / q)):
            raise ValueError("initX does not follow a balanced design.")

    nnp = xp.shape[0]
    results = SATA_AUD(xp.tolist(), n - nnp, s, q, init, initX.tolist(), crit, maxiter, hits_ratio, levelpermt, rand_seed)
    stat = {"initial_design": np.array(np.round(results.Init_Design), dtype=int),
            "final_design": np.array(np.round(results.Final_Design), dtype=int),
            "initial_criterion": results.Init_Obj,
            "final_criterion": results.Final_Obj,
            "time_consumed": results.Time_Second,
            "criterion_history": np.array(results.Criterion_history)}
    if vis:
        plt.figure(figsize=(10, 6))
        plt.plot(stat['criterion_history'], color="black", linewidth=1)
        bst_score = round(stat['final_criterion'], 5)
        min_index = np.argmin(stat['criterion_history'])
        plt.axvline(min_index, color="red", linewidth=1)
        plt.axhline(stat['criterion_history'][min_index], linewidth=1)
        plt.title("Best value = " + str(bst_score) + " in " + str(round(stat['time_consumed'], 3)) + " sec", fontsize=20)
        plt.xlabel("Iterations", fontsize=18)
        plt.ylabel("Criteria", fontsize=18)
        yticks = np.round(np.linspace(np.min(stat['criterion_history']), np.max(stat['criterion_history']), 4), 4)
        plt.yticks(yticks, fontsize=14)
        plt.xticks(fontsize=14)
    return stat


def gen_aud_col(xp, n, s, q, init="rand", initX=np.array([[]]), crit="CD2", maxiter=100, hits_ratio=0.1, levelpermt=False, rand_seed=0, vis=False):
    """
    This function takes n,s,q; a unchanged initial design and other arguments to output a list (described below).

    Parameters
    ----------
    :type  xp: a numpy integer matrix object
    :param xp: representing the previous existing design matrix

    :type  n: an integer object
    :param n: run of experiments

    :type  s: an integer object
    :param s: number of experimental factors, including the previous design in xp

    :type  q: an integer object
    :param q: number of experimental levels for each factor

    :type  crit: a character object, default="CD2"
    :param crit: criterion to be optimized:

             "CD2" -- Centered L2 Discrepancy;

             "WD2" -- Wrap-around L2 Discrepancy;

             "MD2" -- Mixture L2 Discrepancy;

             "maximin" -- Maximin Discrepancy;

             "MC" -- Minimum Coherence;

             "A2" -- Mean Squared Correlation.

    :type  init: a string vector object, default="rand"
    :param init: initialization method for the factor-augmented design:

              "rand": randomly generate initial design;

              "input": user specified.

    :type  initX: a user-defined numpy integer matrix object, default=np.array([[]])
    :param initX: This is the user-defined initial design matrix, and will be used when init="input"

    :type  maxiter: a positive integer object, default=100
    :param maxiter: maximum iteration number in outer while loop of SATA algorithm.

    :type  levelpermt: a boolean object, default=False
    :param levelpermt: it controls whether to use level permutation

    :type  hits_ratio: a float object, default=0.1
    :param hits_ratio: Default value is 0.1, which is the ratio to accept changes of design in inner for loop

    :type  rand_seed: an integer object, default=0
    :param rand_seed: random seed

    :type vis: a boolean object, default=False
    :param vis: if true, plot the criterion value sequence

    Example
    --------
    >>> n=3
    >>> s=4
    >>> q=3
    >>> xp = np.array([[1, 1],
    >>>           [2, 2],
    >>>           [3, 3]])
    >>> crit = "CD2"
    >>> res = gen_aud_col(xp,n,s,q,crit=crit,maxiter=100,vis = True)

    """

    # check the arguments
    if ((isinstance(n, int) & isinstance(s, int) & isinstance(q, int)) is False):
        raise ValueError("Wrong types of n,s,q.")
    elif ((n % q) != 0):
        raise ValueError("n should be multiple of q.")
    elif (n > q**s):
        raise ValueError("n should not be greater than q^s.")
    elif ((s < 1) | (n < 2) | (q < 2)):
        raise ValueError("Invalid design table.")
    elif (not isinstance(xp, np.ndarray)):
        raise ValueError("xp must be a numpy array.")
    elif ((n != xp.shape[0]) | (s <= xp.shape[1])):
        raise ValueError("The size of the existing design matrix xp does not match the given n,s.")
    elif ((xp.dtype != np.int) | (1 != np.min(xp)) | (q != np.max(xp))):
        raise ValueError("The values of the existing design matrix x0 should be integers within: 1,2,3...,q.")
    elif (any(np.array([v for i in range(xp.shape[1]) for v in np.unique(xp[:, i], return_counts=True)[1]]) > n / q)):
        raise ValueError("xp does not follow a balanced design.")
    elif (init == "input"):
        if ((not isinstance(initX, np.ndarray))):
            raise ValueError("initX must be a numpy array.")
        elif ((n != initX.shape[0]) | (s <= initX.shape[1])):
            raise ValueError("The size of the input design matrix does not match the given n,s.")
        elif ((1 != np.min(initX)) | (q != np.max(initX)) | (initX.dtype != np.int)):
            raise ValueError("The values of the input design matrix initX should be integers within: 1,2,3...,q.")
        elif (any(np.array([v for i in range(initX.shape[1]) for v in np.unique(initX[:, i], return_counts=True)[1]]) > n / q)):
            raise ValueError("initX does not follow a balanced design.")

    nvp = xp.shape[1]
    results = SATA_AUD_COL(xp.tolist(), s - nvp, q, init, initX.tolist(), crit, maxiter, hits_ratio, levelpermt, rand_seed)
    stat = {"initial_design": np.array(np.round(results.Init_Design), dtype=int),
            "final_design": np.array(np.round(results.Final_Design), dtype=int),
            "initial_criterion": results.Init_Obj,
            "final_criterion": results.Final_Obj,
            "time_consumed": results.Time_Second,
            "criterion_history": np.array(results.Criterion_history)}
    if vis:
        plt.figure(figsize=(10, 6))
        plt.plot(stat['criterion_history'], color="black", linewidth=1)
        bst_score = round(stat['final_criterion'], 5)
        min_index = np.argmin(stat['criterion_history'])
        plt.axvline(min_index, color="red", linewidth=1)
        plt.axhline(stat['criterion_history'][min_index], linewidth=1)
        plt.title("Best value = " + str(bst_score) + " in " + str(round(stat['time_consumed'], 3)) + " sec", fontsize=20)
        plt.xlabel("Iterations", fontsize=18)
        plt.ylabel("Criteria", fontsize=18)
        yticks = np.round(np.linspace(np.min(stat['criterion_history']), np.max(stat['criterion_history']), 4), 4)
        plt.yticks(yticks, fontsize=14)
        plt.xticks(fontsize=14)
    return stat


def gen_ud_ms(n, s, q, crit="CD2", maxiter=100, nshoot=5, rand_seed=0, n_jobs=1, vis=False):
    """
    This function generates Uniform Design of Experiments using diffrent initializations.

    Parameters
    ----------
    :type  xp: a numpy integer matrix object
    :param xp: representing the previous existing design matrix

    :type  n: an integer object
    :param n: run of experiments

    :type  s: an integer object
    :param s: number of experimental factors

    :type  q: an integer object
    :param q: number of experimental levels for each factor
    
    :type  crit: a character object, default="CD2"
    :param crit: criterion to be optimized:

             "CD2" -- Centered L2 Discrepancy;

             "WD2" -- Wrap-around L2 Discrepancy;

             "MD2" -- Mixture L2 Discrepancy;

             "maximin" -- Maximin Discrepancy;

             "MC" -- Minimum Coherence;

             "A2" -- Mean Squared Correlation.

    :type  maxiter: a positive integer object, default=100
    :param maxiter: maximum iteration number in outer while loop of SATA algorithm.

    :type  nshoot: a positive integer object, default=5
    :param nshoot: total counts to try different initial designs

    :type  rand_seed: an integer object, default=0
    :param rand_seed: random seed

    :type  hits_ratio: a float object, default=0.1
    :param hits_ratio: Default value is 0.1, which is the ratio to accept changes of design in inner for loop

    :type  rand_seed: an integer object, default=0
    :param rand_seed: random seed

    :type  n_jobs: an integer object, default=1
    :param n_jobs: the number of cores to be used for parallelization

    :type vis: a boolean object, default=False
    :param vis: if true, plot the criterion value sequence
    """


    if ((isinstance(n, int) & isinstance(s, int) & isinstance(q, int)) is False):
        raise ValueError("Wrong types of n,s,q.")
    elif ((n % q) != 0):
        raise ValueError("n should be multiple of q.")
    elif (n > q**s):
        raise ValueError("n should not be greater than q^s.")
    elif ((s < 1) | (n < 2) | (q < 2)):
        raise ValueError("Invalid design table.")

    crit_list = []
    shoot_idx = []
    time_list = []
    best_crit = 1e10
    
    stats = Parallel(n_jobs=n_jobs)(delayed(gen_ud)(n=n, s=s, q=q, crit=crit, 
                                    maxiter=maxiter, rand_seed=rand_seed + i) for i in range(nshoot))
    for i in range(nshoot):
        stat = stats[i]
        crit_list.append(list(stat["criterion_history"]))
        shoot_idx.append(len(crit_list))
        time_list.append(stat["time_consumed"])
        tmp = design_eval(stat["final_design"], crit=crit)
        if (tmp < best_crit):
            best_crit = tmp
            best_design = np.array(np.round(stat["final_design"]), dtype=int)
    if vis:
        fig, axes = plt.subplots(int(np.ceil(1.0 * nshoot / 5)), 5, figsize=(15, 4 * np.ceil(1.0 * nshoot / 5)), sharex=True, sharey=True)
        plt.setp(axes, xticks=[], yticks=[])
        axbig = fig.add_subplot(111)
        axbig.set_xticks([])
        axbig.set_yticks([])
        axbig.set_xlabel("Iteration", fontsize=14)
        axbig.set_ylabel("Criteria Value", fontsize=14)
        axbig.xaxis.set_label_coords(0.5, -0.1 / np.ceil(1.0 * nshoot / 5))
        axbig.yaxis.set_label_coords(-0.05, 0.5)
        plt.grid(False)
        for i in range(nshoot):
            ax = fig.add_subplot(np.ceil(1.0 * nshoot / 5), 5, i + 1)
            ax.plot(crit_list[i])
            ax.set_title("Shoot: " + str(i + 1))
            ax.title.set_position((0.5, 0.9))
            ax.set_xticks(np.linspace(0, maxiter, 3, dtype=int))
            if i % 5 != 0:
                ax.set_yticks([])
            if np.min(crit_list[i]) == best_crit:
                ax.axvline(np.argmin(crit_list[i]), color="red", linewidth=1)
        fig.suptitle("Best value = " + str(best_crit) + " in " + str(round(np.sum(time_list), 3)) + " sec",
                     fontsize=20,
                     x=0.5, y=1 - 0.03 * np.ceil(1.0 * nshoot / 5))
        fig.subplots_adjust(wspace=0)
    return best_design


def gen_aud_ms(xp, n, s, q, crit="CD2", maxiter=100, nshoot=5, rand_seed=0, n_jobs=1, vis=False):
    """
    This function generates sequential Uniform Design of Experiments (Augmenting Runs) using diffrent initializations.

    Parameters
    ----------
    :type  xp: a numpy integer matrix object
    :param xp: representing the previous existing design matrix

    :type  n: an integer object
    :param n: run of experiments, including the previous design in xp

    :type  s: an integer object
    :param s: number of experimental factors

    :type  q: an integer object
    :param q: number of experimental levels for each factor

    :type  crit: a character object, default="CD2"
    :param crit: criterion to be optimized:

             "CD2" -- Centered L2 Discrepancy;

             "WD2" -- Wrap-around L2 Discrepancy;

             "MD2" -- Mixture L2 Discrepancy;

             "maximin" -- Maximin Discrepancy;

             "MC" -- Minimum Coherence;

             "A2" -- Mean Squared Correlation.

    :type  maxiter: a positive integer object, default=100
    :param maxiter: maximum iteration number in outer while loop of SATA algorithm.

    :type  nshoot: a positive integer object, default=5
    :param nshoot: total counts to try different initial designs

    :type  rand_seed: an integer object, default=0
    :param rand_seed: random seed

    :type  hits_ratio: a float object, default=0.1
    :param hits_ratio: Default value is 0.1, which is the ratio to accept changes of design in inner for loop

    :type  rand_seed: an integer object, default=0
    :param rand_seed: random seed
    
    :type  n_jobs: an integer object, default=1
    :param n_jobs: the number of cores to be used for parallelization

    :type vis: a boolean object, default=False
    :param vis: if true, plot the criterion value sequence

    """

    if ((isinstance(n, int) & isinstance(s, int) & isinstance(q, int)) is False):
        raise ValueError("Wrong types of n,s,q.")
    elif ((n % q) != 0):
        raise ValueError("n should be multiple of q.")
    elif (n > q**s):
        raise ValueError("n should not be greater than q^s.")
    elif ((s < 1) | (n < 2) | (q < 2)):
        raise ValueError("Invalid design table.")
    elif (not isinstance(xp, np.ndarray)):
        raise ValueError("xp must be a numpy array.")
    elif ((n <= xp.shape[0]) | (s != xp.shape[1])):
        raise ValueError("The size of the existing design matrix xp does not match the given n,s.")
    elif ((xp.dtype != np.int) | (1 > np.min(xp)) | (q < np.max(xp))):
        raise ValueError("The values of the existing design matrix x0 should be integers within: 1,2,3...,q.")
    elif (any(np.array([v for i in range(xp.shape[1]) for v in np.unique(xp[:, i], return_counts=True)[1]]) > n / q)):
        raise ValueError("xp does not follow a balanced design, please increase the number of n or remove duplicated elements (per column) in xp.")

    crit_list = []
    shoot_idx = []
    time_list = []
    best_crit = 1e10
    
    stats = Parallel(n_jobs=n_jobs)(delayed(gen_aud)(xp=xp, n=n, s=s, q=q, crit=crit, 
                                     maxiter=maxiter, rand_seed=rand_seed + i) for i in range(nshoot))
    for i in range(nshoot):
        stat = stats[i]
        crit_list.append(list(stat["criterion_history"]))
        shoot_idx.append(len(crit_list))
        time_list.append(stat["time_consumed"])
        tmp = design_eval(stat["final_design"], crit=crit)
        if (tmp < best_crit):
            best_crit = tmp
            best_design = np.array(np.round(stat["final_design"]), dtype=int)
    if vis:
        fig, axes = plt.subplots(int(np.ceil(1.0 * nshoot / 5)), 5, figsize=(15, 4 * np.ceil(1.0 * nshoot / 5)), sharex=True, sharey=True)
        plt.setp(axes, xticks=[], yticks=[])
        axbig = fig.add_subplot(111)
        axbig.set_xticks([])
        axbig.set_yticks([])
        axbig.set_xlabel("Iteration", fontsize=14)
        axbig.set_ylabel("Criteria Value", fontsize=14)
        axbig.xaxis.set_label_coords(0.5, -0.1 / np.ceil(1.0 * nshoot / 5))
        axbig.yaxis.set_label_coords(-0.05, 0.5)
        plt.grid(False)
        for i in range(nshoot):
            ax = fig.add_subplot(np.ceil(1.0 * nshoot / 5), 5, i + 1)
            ax.plot(crit_list[i])
            ax.set_title("Shoot: " + str(i + 1))
            ax.title.set_position((0.5, 0.9))
            ax.set_xticks(np.linspace(0, maxiter, 3, dtype=int))
            if i % 5 != 0:
                ax.set_yticks([])
            if np.min(crit_list[i]) == best_crit:
                ax.axvline(np.argmin(crit_list[i]), color="red", linewidth=1)
        fig.suptitle("Best value = " + str(best_crit) + " in " + str(round(np.sum(time_list), 3)) + " sec",
                     fontsize=20,
                     x=0.5, y=1 - 0.03 * np.ceil(1.0 * nshoot / 5))
        fig.subplots_adjust(wspace=0)
    return best_design


def gen_aud_col_ms(xp, n, s, q, crit="CD2", maxiter=100, nshoot=5, rand_seed=0, n_jobs=1, vis=False):
    """
    This function generates sequential Uniform Design of Experiments (Augmenting Factors) using diffrent initializations.

    Parameters
    ----------
    :type  xp: a numpy integer matrix object
    :param xp: representing the previous existing design matrix

    :type  n: an integer object
    :param n: run of experiments

    :type  s: an integer object
    :param s: number of experimental factors, including the number of factors in previous design xp

    :type  q: an integer object
    :param q: number of experimental levels for each factor

    :type  crit: a character object, default="CD2"
    :param crit: criterion to be optimized:

             "CD2" -- Centered L2 Discrepancy;

             "WD2" -- Wrap-around L2 Discrepancy;

             "MD2" -- Mixture L2 Discrepancy;

             "maximin" -- Maximin Discrepancy;

             "MC" -- Minimum Coherence;

             "A2" -- Mean Squared Correlation.

    :type  maxiter: a positive integer object, default=100
    :param maxiter: maximum iteration number in outer while loop of SATA algorithm.

    :type  nshoot: a positive integer object, default=5
    :param nshoot: total counts to try different initial designs

    :type  rand_seed: an integer object, default=0
    :param rand_seed: random seed

    :type  hits_ratio: a float object, default=0.1
    :param hits_ratio: Default value is 0.1, which is the ratio to accept changes of design in inner for loop

    :type  rand_seed: an integer object, default=0
    :param rand_seed: random seed

    :type  n_jobs: an integer object, default=1
    :param n_jobs: the number of cores to be used for parallelization
    
    :type vis: a boolean object, default=False
    :param vis: if true, plot the criterion value sequence

    """

    if ((isinstance(n, int) & isinstance(s, int) & isinstance(q, int)) is False):
        raise ValueError("Wrong types of n,s,q.")
    elif ((n % q) != 0):
        raise ValueError("n should be multiple of q.")
    elif (n > q**s):
        raise ValueError("n should not be greater than q^s.")
    elif ((s < 1) | (n < 2) | (q < 2)):
        raise ValueError("Invalid design table.")
    elif (not isinstance(xp, np.ndarray)):
        raise ValueError("xp must be a numpy array.")
    elif ((n != xp.shape[0]) | (s <= xp.shape[1])):
        raise ValueError("The size of the existing design matrix xp does not match the given n,s.")
    elif ((xp.dtype != np.int) | (1 != np.min(xp)) | (q != np.max(xp))):
        raise ValueError("The values of the existing design matrix x0 should be integers within: 1,2,3...,q.")
    elif (any(np.array([v for i in range(xp.shape[1]) for v in np.unique(xp[:, i], return_counts=True)[1]]) > n / q)):
        raise ValueError("xp does not follow a balanced design.")

    crit_list = []
    shoot_idx = []
    time_list = []
    best_crit = 1e10
    
    
    stats = Parallel(n_jobs=n_jobs)(delayed(gen_aud_col)(xp=xp, n=n, s=s, q=q, crit=crit, 
                                     maxiter=maxiter, rand_seed=rand_seed + i) for i in range(nshoot))
    for i in range(nshoot):
        stat = stats[i]
        crit_list.append(list(stat["criterion_history"]))
        shoot_idx.append(len(crit_list))
        time_list.append(stat["time_consumed"])
        tmp = design_eval(stat["final_design"], crit=crit)
        if (tmp < best_crit):
            best_crit = tmp
            best_design = np.array(np.round(stat["final_design"]), dtype=int)
    if vis:
        fig, axes = plt.subplots(int(np.ceil(1.0 * nshoot / 5)), 5, figsize=(15, 4 * np.ceil(1.0 * nshoot / 5)), sharex=True, sharey=True)
        plt.setp(axes, xticks=[], yticks=[])
        axbig = fig.add_subplot(111)
        axbig.set_xticks([])
        axbig.set_yticks([])
        axbig.set_xlabel("Iteration", fontsize=14)
        axbig.set_ylabel("Criteria Value", fontsize=14)
        axbig.xaxis.set_label_coords(0.5, -0.1 / np.ceil(1.0 * nshoot / 5))
        axbig.yaxis.set_label_coords(-0.05, 0.5)
        plt.grid(False)
        for i in range(nshoot):
            ax = fig.add_subplot(np.ceil(1.0 * nshoot / 5), 5, i + 1)
            ax.plot(crit_list[i])
            ax.set_title("Shoot: " + str(i + 1))
            ax.title.set_position((0.5, 0.9))
            ax.set_xticks(np.linspace(0, maxiter, 3, dtype=int))
            if i % 5 != 0:
                ax.set_yticks([])
            if np.min(crit_list[i]) == best_crit:
                ax.axvline(np.argmin(crit_list[i]), color="red", linewidth=1)
        fig.suptitle("Best value = " + str(best_crit) + " in " + str(round(np.sum(time_list), 3)) + " sec",
                     fontsize=20,
                     x=0.5, y=1 - 0.03 * np.ceil(1.0 * nshoot / 5))
        fig.subplots_adjust(wspace=0)
    return best_design
