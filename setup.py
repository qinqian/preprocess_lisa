"""
MARGE2
-------

TF enrichment ranking through cis-element knock-out

"""

from setuptools import setup, find_packages
from glob import glob, os

os.system("cp lisa.conf lib/m2")

setup(
    name='LISA',
    version='0.0',
    packages=find_packages('lib'),
    package_dir={'': 'lib'},
    package_data={'': ['*.conf', 'template/*.css', "template/*.html", "template/*.js", 'template/*genome', 'template/Snakemake', 'template/*py', 'template/conda.yaml', 'template/*gz']},
    include_package_data=True,
    py_modules=['marge2_conf', 'marge2'],
    install_requires=['numpy', 'matplotlib', 'sklearn', 'theano', 'h5py', 'pandas', 'seaborn', 'scipy', 'snakemake', 'PyYAML', 'statsmodels'],
    author='Qian Qin',
    scripts=glob('scripts/*'),
    classifiers=[
        'Environment::Console',
        'Operating System:: POSIX',
        "Programming Language :: Python :: 2",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules"
    ],
    zip_safe = False
)
