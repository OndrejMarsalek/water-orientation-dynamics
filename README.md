Orientational structure and dynamics of water around solutes
============================================================

Orientation-resolved distribution functions
-------------------------------------------

Program: `bin/or-df`

Plotting tool: `bin/plot-or-df`


Conditional time correlation functions
--------------------------------------

Program: `bin/ctcf`

Plotting tool: `bin/plot-ctcf`


Build
-----

Building the analysis tools requires an installation of GROMACS 5.1 and a C++
compiler. Note that earlier versions of GROMACS (including 5.0.x) will not
work, as there have been changes in header files and the signatures of some
functions. The same goes for newer versions of GROMACS, unfortunately. If/when
the constant changes are over, I would consider an update.

Edit `makefile` to set the compiler and compiler flags under the "user build
settings" section or pick one of the suggested settings. Make sure that you
have sourced `GMXRC` to have the paths available.

Compile with `make`. There is no installation, the programs are used from the
directory where they are built. If you want to call them conveniently, source
the `env.sh` file, which updates your executable path.
