# pyUniDOE

[![Build Status](https://travis-ci.com/ZebinYang/pyUniDOE.svg?branch=master)](https://travis-ci.org/joerick/cibuildwheel)

## Installation

- Enviroment: 
    - Python 3
    - C++ compiler
    - swig 3
    
Assume you have figured out the above environment, the most convenient way for installation is via the pip command. 
```sheel
pip install git+https://github.com/ZebinYang/pyunidoe.git.
```

More details can be found in [documentation](https://zebinyang.github.io/pyunidoe/build/html/index.html).

## Examples

### Evaluate existing designs
```python
x = np.array([[1, 2],
              [3, 3],
              [2, 1]])
pydoe.design_eval(x,crit="CD2")
```

### Get an existing design from database
```python 
import numpy as np 
import pyunidoe as pydoe

pydoe.design_query(n=12, s=4, q=6, crit="CD2", show_crit=True)
```

### Generate uniform design from random initialization
```python 
stat=pydoe.gen_ud(n=12, s=4, q=6, init="rand", crit="CD2", maxiter=100, vis=True)
print("The initial design: ")
print(stat["initial_design"])

print("The final design: ")
print(stat["final_design"])
```

### Augment uniform design (Runs)
```python 
stat = pydoe.gen_aud(xp=x1, n=24, s=4, q=6, crit="CD2", maxiter=100, vis=True)
print("The initial design: ")
print(stat["initial_design"])

print("The final design: ")
print(stat["final_design"])
```

### Augment uniform design (Factors)
```python 
stat = pydoe.gen_aud_col(xp=x1, n=12, s=5, q=6, crit="CD2", maxiter=100, vis=True)
print("The initial design: ")
print(stat["initial_design"])

print("The final design: ")
print(stat["final_design"])
```

### Multi-shoot Strategy
```python 
x1_multi = pydoe.gen_ud_ms(n=12, s=4, q=6, crit="CD2", maxiter=100, nshoot=1000, n_jobs=10, vis=False)
print(pydoe.design_eval(x1_multi,crit="CD2"))

x2_multi = pydoe.gen_aud_ms(x1_multi, n=24, s=4, q=6, crit="CD2", maxiter=100, nshoot=1000, n_jobs=10, vis=False)
print(pydoe.design_eval(x2_multi,crit="CD2"))

x3_multi = pydoe.gen_aud_col_ms(x1_multi, n=12, s=5, q=6, crit="CD2", maxiter=100, nshoot=1000, n_jobs=10, vis=False)
print(pydoe.design_eval(x3_multi,crit="CD2"))
```

More examples can be referred to the [documentation](https://zebinyang.github.io/pyunidoe/build/html/examples.html)

## Contact:
If you find any bugs or have any suggestions, please contact us via email: yangzb2010@hku.hk or ajzhang@hku.hk.

## Reference:
Zebin Yang, Aijun Zhang and Ji Zhu. (2019) Hyperparameter Optimization via Sequential Uniform Designs. Submitted. 
