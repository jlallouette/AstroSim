# AstroSim

AstroSim is an astrocyte network simulation software. It implements astrocyte calcium dynamics models described in:
- [Sparse short-distance connections enhance calcium wave propagation in a 3D model of astrocyte networks](http://journal.frontiersin.org/article/10.3389/fncom.2014.00045/full)
- [Glutamate Mediated Astrocytic Filtering of Neuronal Activity](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003964)

## Dependencies

The following packages are required:
- xutils-dev
- libgsl0-dev
- libboost-filesystem-dev
- libalglib-dev

## Compilation

```
make depend
make
```

## Usage

### Parameters

AstroSim is a command-line software, simulation parameters are passed through arguments or by providing a file containing a list of arguments. All the available arguments can be displayed using:
```
./AstroSim -help
```
This list of arguments is not meant to be humanly readable, it's better used in conjunction with a small GUI program that parses it and allows you to fill in parameters values.
This program is located in ./utility/. It depends on [Qt](https://wiki.qt.io/Main) and can be built by:
```
qmake
make
```
You then need to create the parameter description file with e.g.:
```
./AstroSim -help > ParamDescr.ini
```
After having launched the GUI program, go to File -> Open Parameters Description File and select ParamDescr.ini. This initializes the GUI with all possible arguments to AstroSim. Since no description is provided, it's better to first look in the code which parameters to set. Parameter sets can then be saved and loaded with File -> Save / Load Parameters. Example of parameter files are provided in utility/parameters.

To launch a simulation using a parameter file, just type:
```
./AstroSim -useParams Path/To/ParamFile
```
### Saved Data
Using the -aM (add metric) and -sD (save metric) options, one can choose which values will be computed and which values will be saved. The GUI program shows a full list of available values for these options. To know in details which metric does what, the best is to check the code in source files named \*Metrics.cpp (NetworkMetrics.cpp, PropagationMetrics.cpp, etc.).

The created files are saved by default in `./data`, this can be modified using the -Path option.

