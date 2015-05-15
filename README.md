# dna_ode

Simulating DNA using coarse-grained Langevin molecular dynamics with ODE (Open Dynamics Engine).

This program allows to simulate a DNA molecule under optical or magnetic tweezers. The DNA molecule is modelled as chain of cylinders connected by rigid articulations (joints), and using the Open Dynamics Engine.

## Dependencies

* ODE: Open Dynamics Engine, version >= 0.13. Currently maintained at [https://bitbucket.org/odedevs/ode/](https://bitbucket.org/odedevs/ode/). The documentation is available at the [ODE homepage](http://ode-wiki.org/wiki/index.php?title=Main_Page).
* GSL: Gnu Scientific Library, version >= 0.16, available at [http://www.gnu.org/software/gsl/](http://www.gnu.org/software/gsl/).
* libconfig: Simple C library for reading configuration files, version >= 1.4.9. Available at [www.hyperrealm.com/libconfig/](www.hyperrealm.com/libconfig/).

## Installation
To compile and install this package, other than the external libraries mentioned above you normally need only a C99-compliant C compiler:

1. Download the full distribution from `dists/dna_ode-(version).tar.gz`

2. Unpack it in a new directory with `tar -xf dna_ode-(version).tar.gz`

3. Run `./configure` (with your favorite options... documented elsewhere)

4. Run `make`. The package should compile without complaining if everything is okay.

5. Run `make install` if you want to install the program system-wide.

If you clone the code directly from GitHub, you also need libtool and autotools. In such case, before step 3, you should run `./bootstrap.sh` to generate the necessary initialization files.

## Usage
The basic command line invocation for running a simulation is:
```
dna_ode run <configuration_file> <seed>
```
where the `<configuration_file>` contains all the parameters of the simulation. Please see the "examples" directory after compiling the package for example configuration files. The `<seed>` is the seed for the random number generator used in the program.

After a simulation is finished, you can run a few analyses by typing
```
dna_ode analyze <configuration_file> ...
```
Run this command without further arguments to see a full list of possible analysis types. This command looks at the simulation trajectory file specified in the configuration file, and performs the requested analysis type. Output varies depending on the analysis type.

It is useful when scripting to be able to extract the value of a parameter from a configuration file, so we provided a
```
dna_ode getp <configuration_file> <parameter_name>
```
which returns the value of `<parameter_name>` from `<configuration_file>`.


## References
Currently two papers were published using this software:

1. P. Carrivain, M. Barbi, and J.-M. Victor, "In silico single-molecule manipulation of DNA with rigid body dynamics", *PLoS Computational Biology* **10**(2), e1003456 (2014). DOI: [10.1371/journal.pcbi.1003456](https://dx.doi.org/10.1371/journal.pcbi.1003456)

2. R. Cortini, B. Caré, J.-M. Victor and M. Barbi, "Theory and simulations of toroidal and rod-like structures in single-molecule DNA condensation", *The Journal of Chemical Physics*, **142**, 105102-1–9, (2015). DOI: [10.1063/1.4914513](https://dx.doi.org/10.1063/1.4914513)

