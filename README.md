# Stellar Spin Rates Are Awesome

This package includes routines and notes from my research on the spin evolution of low mass stars. 

The jupyter-notebook `SpinEvolutionModel.ipynb` contains my research notes on the implementation  Matt et al. 2015 spin evolution model along with modifications I added to it in order to include the rotational evolution during the early pre-Main Sequence. 

The package `SpinRatesAreHere.py` allows anyone to load and manipulate datasets of spin rates of low-mass stars.

This initial version contains datasets for the following clusters:

- NGC6530

- NGC2264

- NGC2362

- UpperSco

- hPer

- Pleiades

- Praesepe

- NGC6811

I will soon include the following other clusters (tables and literature review are ready, just need formating):


- ONC

- CepOB3b

- Orion OB1

- Orion OB1a

- Orion OB1b

- 25 Ori

- Rho Oph


In the future I will also include:

- Cygnus OB2

- Taurus

- older clusters and the Kepler Field.

You can load the data for a cluster by doing:

 `import SpinRatesAreHere as spin`

`ngc2264=spin.NGC2264()`

The keyword `mass_type` can be set as: `0` for the masses in the source paper, `1` for mass estimations using MESA and `2` for mass estimations using Baraffe et al. 98.

The table with the full dataset can be found in `ngc2264.data`.

`Disk` will give you a flag say if the star has disk or not.

the directory tables/ will include all the tables in the fomat required for my package to work.

the directory ipynb contains jupyter notebooks with notes about how I built this package. tables/info_init/ includes the notes about how I built the database and some of the references for the source papers. tables/StellarParameters/ includes notes about how I did the mass transformations and estimated reddenings.

Note that! this package is still under construction. If you found out that I am missing some reference or if you find any sort of mistake, please let me know at jt574@exeter.ac.uk.
