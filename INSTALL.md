# Installation Directions for BU-COVID pipeline


## R Requirements

The code has been tested with R 3.6 and 4.0 versions.
The following packages are used in R scripts: igraph, readr, dplyr, ggplot2, tidyr, scales, Matrix

## Install bu_covid Python code

This code has been tested with Python 3.6 and higher.  The `--user` flag here installs bu_covid to your home directory. This flag can be omitted if you wish to install to the default Python package location. This will install the latest version of the Covasim library, the Python igraph library, and all required dependencies.  

```
cd BU_COVID/python/bu_covid
python setup.py install --user
```
You can check that the bu_covid and Covasim libraries are loading correctly with an import test:
```
python -c 'import bu_covid'
```
 
