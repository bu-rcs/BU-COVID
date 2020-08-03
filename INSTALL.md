# Installation Directions for BU-COVID pipeline


## R Requirements

List of required R packages and minimum versions here.


## Install bu_covid Python code

This code has only been tested with Python 3.6 and higher.  The `--user` flag here installs bu_covid to your home directory. This flag can be 

```
cd BU_COVID/python/bu_covid
python setup.py --user
```
You can check that the bu_covid and Covasim libraries are loading correctly with an import test:
```
python -c 'import bu_covid'
```

This will install the latest version of the Covasim library, the Python igraph library, and all required dependencies.  The Covasim code requires a minor patch to improve how it handles plots.  To apply the patch to Covasim's run.py file use the following command:
```
patch `python -c 'import covasim ; import sys ; sys.stderr.write(covasim.run.__file__)' > /dev/null` < covasim_patch.txt
```



