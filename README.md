# A repo that contains miscellaneous Python programs and a 'palos' module by the yfish group.

This repository is mixture of a python module 'palos' and other standalone programs developed and used by the yfish group, http://www.yfish.org/.

It contains code related to bioinformatics projects focusing on next-generation sequencing data, population genetics, genome-wide association studies, pedigree genetics, etc.

[palos/](palos/) is the source of the [https://pypi.org/project/palos](https://pypi.org/project/palos) module. 

[palos/algorithm/](palos/algorithm/) contains pure algorithms, not specific to Bioinformatics.


[GADA/](GADA/) contains a fast segmentation algorithm/program. It contains some bugfixes and improvement over the original GADA (2008/2009) by using a Red-Black tree to accelerate its speed.

[ngs/](ngs/) contains programs analyzing next-generation sequencing data.

# Prerequisites
The `palos` module is universally used by many programs and must be installed before running the standalone programs.

## Install palos

```sh
pip3 install --upgrade palos
```

# Usage


```sh
./ngs/DownsampleWorkflow.py  -h
```
