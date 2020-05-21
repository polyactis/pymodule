# A repo that contains miscellaneous Python programs and a 'palos' module by the yfish group.

This repository is mixture of a python module 'palos' and other standalone programs developed and used by the yfish group, http://www.yfish.org/.

It contains code related to bioinformatics projects focusing on next-generation sequencing data, population genetics, genome-wide association studies, pedigree genetics, etc.

[palos/](palos/) is the source of the [https://pypi.org/project/palos](https://pypi.org/project/palos) module. 

[palos/algorithm/](palos/algorithm/) contains pure algorithms, not specific to Bioinformatics.


[GADA/](GADA/) contains a fast segmentation algorithm/program. It contains some bugfixes and improvement over the original GADA (2008/2009) by using a Red-Black tree to accelerate its speed.

[ngs/](ngs/) contains programs analyzing next-generation sequencing data.

# Prerequisites
The `palos` module is universally used by many programs and must be installed before running any program in pymodule.

## Optional prerequisites

Following pakcages are optional, only needed for some functions.

1. python-mysqldb
1. python-psycopg2 http://initd.org/psycopg/
1. matplotlib basemap toolkit http://matplotlib.sourceforge.net/basemap/doc/html/
1. python imaging library http://www.pythonware.com/products/pil/
1. python-scientific http://www.scipy.org/
1. biopython
1. python-rpy
1. networkx https://networkx.lanl.gov/wiki
1. hcluster
1. python-h5py
1. python-tables

## Optional C++ libraries

Required if you plan to compile all binaries in pymodule by typeing 'make all'.

apt-get install libhdf5-dev libhdf5-serial-dev libhdf5-cpp-100 hdf5-tools \
       libarmadillo-dev libboost-program-options-dev libboost-iostreams-dev \
       libboost-python-dev python-dev

## Install palos

```sh
pip3 install --upgrade palos
```

# Usage

```sh
./ngs/DownsampleWorkflow.py  -h
```

