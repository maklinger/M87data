# M87data
The code to plot various data of M87: 

- jet collimation 
- jet speed
- jet magnetic field strength

based on figs 8 and 16 of [Hada et al. 2024](https://link.springer.com/10.1007/s00159-024-00155-y)

- SED
    - quiet from [Prieto et al. 2016](https://academic.oup.com/mnras/article-lookup/doi/10.1093/mnras/stw166)
    - [EHT MWL 2017](https://iopscience.iop.org/article/10.3847/2041-8213/abef71)
    - [EHT MWL 2018](https://www.aanda.org/10.1051/0004-6361/202450497)

## Installation with pip
```shell
pip install git+https://github.com/maklinger/M87data.git
```


## Changes and updating version number

1. Update version in M87data/__init__.py
2. Update version in pyproject.toml
3. push news to git
4. new git tag: `git tag v1.1.0 && git push origin v1.1.0`