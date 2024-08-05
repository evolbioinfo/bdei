from __future__ import absolute_import, with_statement, print_function, division

import os

from setuptools import setup, Extension, find_packages
import numpy as np

os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"

# the C++ extension module
bdei_module = Extension('_pybdei',
                        sources=['pybdei/_pybdei.cpp', 'bdei/BDEI.cpp'],
                        include_dirs=['bdei', 'pybdei', np.get_include()],
                        depends=['pybdei/_pybdei.hpp', 'pybdei/_python.hpp', 'bdei/BDEI.hpp', 'bdei/pool.hpp', 'bdei/queue.hpp', 'bdei/semaphore.hpp'],
                        libraries=['nlopt'],
                        extra_compile_args=['-pthread', '-O3', '-std=c++17', '-lnlopt_cxx', '-c'],
                        )

setup(
    name='pybdei',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3 :: Only'
    ],
    packages=find_packages(),
    install_requires=['numpy>=1.24.0', 'ete3>=3.1.3', 'six>=1.16.0', 'scipy>=1.11.1', 'treesimulator>=0.1.22',
                      'nlopt>=2.6.1'],
    setup_requires=['numpy>=1.24.0', 'nlopt>=2.6.1'],
    include_package_data=True,
    package_data={'pybdei': ['*.hpp'],
                  'bdei': ['*.hpp']},
    version='0.9',
    description='Fast and accurate epidemiological parameter estimation from phylogenetic trees '
                'with the Birth-Death Exposed-Infectious (BDEI) model.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Frédéric Heicht',
    author_email='frederic.hecht@sorbonne-universite.fr',
    maintainer='Anna Zhukova',
    maintainer_email='anna.zhukova@pasteur.fr',
    url='https://github.com/evolbioinfo/BDEI',
    download_url='https://github.com/evolbioinfo/BDEI',
    keywords=['BDEI', 'phylodynamics', 'epidemiological parameters'],
    ext_modules=[bdei_module],
    entry_points={
                'console_scripts': [
                    'bdei_infer = pybdei.inference:main',
                    'bdei_loglikelihood = pybdei.loglikelihood:main',
                    'bdei_u = pybdei.u_calculator:main',
                ]
        },
)