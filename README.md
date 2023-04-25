- [1. A repo that contains miscellaneous Python/C++ modules/programs, a standalone Python module 'palos' by the yfish group (can be installed by pip).](#1-a-repo-that-contains-miscellaneous-pythonc-modulesprograms-a-standalone-python-module-palos-by-the-yfish-group-can-be-installed-by-pip)
- [2. Prerequisites to run Python programs in Pymodule](#2-prerequisites-to-run-python-programs-in-pymodule)
  - [2.1. Palos](#21-palos)
  - [2.2. Optional prerequisites](#22-optional-prerequisites)
  - [2.3. Optional C++ libraries](#23-optional-c-libraries)
- [3. Example on how to run some pymodule programs](#3-example-on-how-to-run-some-pymodule-programs)

Yu S. Huang, polyactis@gmail.com

# 1. A repo that contains miscellaneous Python/C++ modules/programs, a standalone Python module 'palos' by the yfish group (can be installed by pip).

This repository is mixture of a python module 'palos' and other standalone programs developed and used by the yfish group, http://www.yfish.org/.

It contains code related to bioinformatics projects focusing on next-generation sequencing data, population genetics, genome-wide association studies, pedigree genetics, etc.

[palos/](palos/) is the source of the [https://pypi.org/project/palos](https://pypi.org/project/palos) module. 

[palos/algorithm/](palos/algorithm/) contains pure algorithms, not specific to Bioinformatics.


[GADA/](GADA/) contains a faster algorithm than the original GADA (2008/2009) by using a Red-Black tree. Now in an independent repo https://github.com/polyactis/eGADA.

[ngs/](ngs/) contains programs analyzing next-generation sequencing data.

# 2. Prerequisites to run Python programs in Pymodule
Most programs in pymodule is dependent on the `palos` module, which is housed in [palos/](palos/). Installation of `palos` will trigger installation of other dependencies.

## 2.1. Palos
Palos supports Python3 primarily, but is ported to Python2 via https://github.com/asottile/future-fstrings, because some pymodule programs are Python2-only.

```sh
pip3 install --upgrade palos
```

```sh
# to run some Python2 pymodule programs
pip install --upgrade palos
```


Package future-fstrings allows the use of f-string in Python2.
```python
# -*- coding: future_fstrings -*-
thing = 'world'
print(f'hello {thing}')
```

## 2.2. Optional prerequisites

The following pakcages are optional, only needed for some functions.

1. mysqldb
2. biopython
3. pegaflow https://pypi.org/project/Pegaflow/
4. psycopg2 http://initd.org/psycopg/
5. matplotlib basemap toolkit http://matplotlib.sourceforge.net/basemap/doc/html/
6. python imaging library http://www.pythonware.com/products/pil/
7. python-scientific http://www.scipy.org/
8. biopython
9. python-rpy2
10. networkx https://networkx.lanl.gov/wiki
11. hcluster
12. python-h5py
13. python-tables

## 2.3. Optional C++ libraries

Required if you plan to compile all binaries in pymodule by typeing 'make all'.

apt-get install libhdf5-dev libhdf5-serial-dev libhdf5-cpp-100 hdf5-tools \
       libarmadillo-dev libboost-program-options-dev libboost-iostreams-dev \
       libboost-python-dev python-dev



# 3. Example on how to run some pymodule programs

```sh
./ngs/DownsampleWorkflow.py  -h
```

