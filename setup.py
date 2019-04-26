import os 
from distutils.sysconfig import *
from setuptools import setup, Extension

cfg_vars = get_config_vars()
cfg_vars['OPT'] = cfg_vars['OPT'].replace("-Wstrict-prototypes", "")
cfg_vars['PY_CFLAGS'] = cfg_vars['PY_CFLAGS'].replace("-Wstrict-prototypes", "")
cfg_vars['CFLAGS'] = cfg_vars['CFLAGS'].replace("-Wstrict-prototypes", "")

pyUniDOE_module=Extension('_pyUniDOE_swig', sources=['pyUniDOE/pyUniDOE_swig_wrap.cxx','pyUniDOE/wrapper.cpp','pyUniDOE/doe_optimizer.cpp','pyUniDOE/doe_criteria.cpp',
                                'pyUniDOE/doe_CD2.cpp','pyUniDOE/doe_MD2.cpp','pyUniDOE/doe_WD2.cpp', 
                                'pyUniDOE/doe_maximin.cpp','pyUniDOE/doe_MC.cpp','pyUniDOE/doe_A2.cpp'])

setup(name='pyUniDOE',
      version='0.1',
      description='Efficient procedures for constructing uniform design of experiments under various space-filling criteria. It is based on a stochastic and adaptive threshold accepting algorithm with flexible initialization, adaptive threshold, and stochastic evolution. The package may also construct the augmented uniform designs in a sequential manner.',
      author='Zebin Yang',
      author_email='yangzb2010@hku.hk',
      license='GPL',
      packages=['pyUniDOE'],
      ext_modules=[pyUniDOE_module],
      package_dir = {'pyUniDOE': '/'},
      package_data = {'pyUniDOE': ['data/*.json']},
      install_requires=['matplotlib', 'numpy', 'pandas', 'seaborn'], zip_safe=False)