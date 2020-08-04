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
The Covasim code requires a minor patch to improve how it handles plots.  By default it will print a vertical dashed line in its own plots wherever there is an intervention applied.  However, due to the way we create our simulations we have interventions occuring on every day to handle the classroom schedules so the vertical lines need to be able to be turned off.  On a Linux system, you can apply a patch to Covasim's run.py to allow the disabling of the vertical lines file using the following command.  This may also work on Mac OSX systems:
```
patch `python -c 'import covasim ; import sys ; sys.stderr.write(covasim.run.__file__)' > /dev/null` < covasim_patch.txt
```
On Windows, it may be easier to edit the run.py file by hand.  To do so print out the location of the file:
```
python -c 'import covasim ; print(covasim.run.__file__)'
```
Edit that file and change line 350 from this:
```
fig = self.base_sim.plot(to_plot=to_plot, colors=colors, **kwargs)
```
to this:
```
fig = self.base_sim.plot(to_plot=to_plot, colors=colors, show_args=show_args, **kwargs)
```


