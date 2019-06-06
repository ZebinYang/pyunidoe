Examples
===============
Here we give more example usage of this package.


Get an existing design from database
---------------------------------------------------
.. code-block::
        pydoe.design_query(n=12,s=4,q=6,crit="CD2", show_crit = True)
        
        
Evaluate existing designs
---------------------------------------------------
.. code-block::
        x = np.array([[1, 2],
              [3, 3],
              [2, 1]])
        pydoe.design_eval(x,crit="CD2")


Generate uniform design from random initialization
---------------------------------------------------

.. code-block::

        import numpy as np 
        import pyunidoe as pydoe
        stat=pydoe.gen_ud(n=12,s=4,q=6,init="rand",crit="CD2",maxiter=100,vis=True)
        print("The initial design: ")
        print(stat["initial_design"])

        print("The final design: ")
        print(stat["final_design"])
        pydoe.design_pairs_plot(stat["final_design"])


Augment uniform design (Runs)
-----------------------------------

.. code-block::

        import numpy as np 
        import pyunidoe as pydoe
        stat=pydoe.gen_ud(n=12,s=4,q=6,init="rand",crit="CD2",maxiter=100,vis=True)
        xp = stat["final_design"]

        stat = pydoe.gen_aud(xp = xp, n = 24, s = 4, q = 6, crit="CD2", maxiter=100, vis = True)
        print("The initial design: ")
        print(stat["initial_design"])

        print("The final design: ")
        print(stat["final_design"])
        pydoe.design_pairs_plot(stat["final_design"])


Augment uniform design (Factors)
-----------------------------------

.. code-block::

        import numpy as np 
        import pyunidoe as pydoe
        stat=pydoe.gen_ud(n=12,s=4,q=6,init="rand",crit="CD2",maxiter=100,vis=True)
        xp = stat["final_design"]

        stat = pydoe.gen_aud_col(xp=xp, n = 12, s = 5 ,q = 6, crit="CD2", maxiter=100, vis = True)
        print("The initial design: ")
        print(stat["initial_design"])

        print("The final design: ")
        print(stat["final_design"])
        pydoe.design_pairs_plot(stat["final_design"])


Multi-shoot Strategy
-----------------------------------

.. code-block::

        import numpy as np 
        import pyunidoe as pydoe
        x1_multi = pydoe.gen_ud_ms(n=12, s=4, q=6, crit="CD2", maxiter=100, nshoot = 5, vis=False)
        pydoe.design_eval(x1_multi,crit="CD2")
        
        x2_multi = pydoe.gen_aud_ms(x1_multi, n=24, s=4, q=6, crit="CD2", maxiter=100, nshoot = 5, vis=False)
        pydoe.design_eval(x2_multi,crit="CD2")
        
        x3_multi = pydoe.gen_aud_col_ms(x1_multi, n=12, s=5, q=6, crit="CD2", maxiter=100, nshoot = 5, vis=False)
        pydoe.design_eval(x3_multi,crit="CD2")