# Installation Directions for BU-COVID pipeline


## R Requirements

List of required R packages and minimum versions here.


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
 
