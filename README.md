# randomice
Generate orientationally disordered configurations of water ice

## What is randomice?
Randomice is a C++ code that generates proton-disordered (or ordered) configurations of water ice that obey the [Bernal-Fowler ice rules](http://https://en.wikipedia.org/wiki/Ice_rules) using algorithms by [Rick](https://aip.scitation.org/doi/full/10.1063/1.1853351) and [Buch *et al*](https://pubs.acs.org/doi/abs/10.1021/jp980866f). It can also generate a slab model with surfaces with an adjustable degree of proton ordering, as described by [Pan *et al*](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.101.155703).

## Dependencies and compiling
Paths to the compiler and GSL and Boost C++ libraries should be specified in Makefile.inc. Compiling should then be as simple as typing "make".

## Usage
    randomice -h
will output a list of command line options. The required in input file should be formatted as a [CASTEP .cell file](https://www.tcm.phy.cam.ac.uk/castep/documentation/WebHelp/content/modules/castep/keywords/k_main_structure.htm) (only the lattice and position sections are required). To generate a supercell with a random ice Ih configuration, the input file should contain the oxygen positions only, then the command,

    randomice -i ice.cell -s 9 5 4

will generate a proton disordered supercell of dimension 9 x 5 x 4.
