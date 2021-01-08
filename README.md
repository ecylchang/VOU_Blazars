# VOU_Blazars

VOU-BLazars is a tool designed to identify blazars in a given sky region based on multi-frequency information retrieved using Virtual Observatory protocols. This is the first version of the tool; a technical paper describing the tool could be found on Astronomy and Computing, 2020, 30, 100350. Any comments and suggestions are welcomed.

## Update

Go to the folder where you download the VOU-Blazars

```bash
$ git pull
```

to update the VOU-Blazars for latest version.

## Install

### Dependencies

See below for dependencies setup.

* [PGPlot](http://www.astro.caltech.edu/~tjp/pgplot/)
* GFortran
* Python 2.7
* [EADA and conesearch pipline](https://github.com/chbrandt/eada)


### Setup EADA/Conesearch pipeline

VOU-Blazars required conesearch pipeline, External Archive Data Access (EADA) Tools, to access VO data services.
If you use [Anaconda Python Distribution](https://www.anaconda.com/download/), you can install EADA in its own (virtual) environment as following.
These steps allow you to create an python conda environment called eada, and all the conesearch/eada jobs are installed and running insides this environment.
See [EADA and conesearch pipline](https://github.com/chbrandt/eada) for more details.

```bash
$ conda create -n eada python=2 pip ipython
$ conda activate eada
$ pip install https://github.com/chbrandt/eada/archive/0.9.7.4.zip
```
Note that the Conessearch pipeline till now only support runing with Python 2.7. And, here is just an example to install eada tool version 0.9.7.4, in future, there might be new version of eada updated.


### Setup PGPlot - MACOS

PGPlot are the main plotting tool in VOU-Blazars. MacOS users are advised to use [brew](https://brew.sh/)
 
```bash
$ brew install pgplot
```

### Compiling Fortran programs

Most of the subroutines in VOU-Blazars are written in Fortran.
Note that you need to install PGPlot before compiling these subroutines.
**Once** the dependencies have been satisfied, go inside `VOU_Blazars/bin/fort` and run the `compile.sh` script:

```bash
$ git clone https://github.com/ecylchang/VOU_Blazars.git
$ cd VOU_Blazars
$ cd bin/fort
$ ./compile.sh
$ cd ../..
```
**Note that if the compile.sh script doesn't work, maybe you need to change the path of the PGPLOT library in compile.sh by seting the path parameter, PGPLOT_DIR**
```bash
$ export PGPLOT_DIR='YOUR_PGPLOT_LIBRARY_DIRECTORY'
```

**Otherwise, you could compile all Fortran programs by**
First compile mylib.f with command.
```bash
$ gfortran -c mylib.f -ffixed-line-length-500
```
and then all the other programs (with .f file name) with command
```bash
$ gfortran -o $XXX $XXX.f -ffixed-line-length-500 mylib.o -L${HOME}/pgplot -lpgplot
```
**Change the PGPLOT path to where you put the PGPLOT library.**

**Always compile all the Fortran routines/subroutines inside the folder "fort"**


## Running VOU-Blazars

> Guarantee EADA's `conesearch` tool is available in your environment -- check this note's section [Setup EADA/Conesearch pipeline](#setup-eadaconesearch-pipeline).

**Before running the tool, make sure that you are under the main folder where the VOU-Blazars are installed and have activated the eada environment** Bin folder should inside the VOU_Blazars folder.

```
$ cd VOU_Blazars
$ conda activate eada
$ ./bin/vou-blazars --ra RA --dec DEC --FOV SEARCH_AREA
```

Usage: vou-blazars { --ra <degrees> --dec <degrees> --area <arcmin> }

ARGUMENTS: <br />
 --ra     : Right Ascension (in DEGREES) <br />
 --dec    : Declination (in DEGREES) <br />
 --FOV   : Searchin Radius (in ARC-MINUTES) around RA,DEC to search for observations <br />

OPTIONS: <br />
--mode    : Running mode <br />
      Options are 'f' find candidate mode (default): finding interesing candidates within a specified region <br />
                          's' SED mode: obtaining SED for a specified source with given R.A. Dec. <br />
                          'l' Light curve mode: obtaining light curve for a specified source with given R.A. Dec. <br />

--nh      : nH column density (in cm^2). Default is 5.e20 cm^2 <br />
            (If the user has installed Heasoft and did not specify the nh, it will use the value calculated by Heasoft) <br />

First and larger error region, circle (--radius) or elliptical (--major --minor --angle) <br />
--radius  : Error circle radius (in ARC-MINUTES). Default is 0 <br />
--major   : Error elliptical major axis (in ARC-MINUTES). Default is 0 <br />
--minor   : Error elliptical minor axis (in ARC-MINUTES). Default is 0 <br />
--angle   : Position angle of the error elliptical (in DEGREES). Default is 0 <br />

Second and smaller error region, circle (--radius2) or elliptical (--major2 --minor2 --angle2) <br />
--radius2 : Second error circle radius (in ARC-MINUTES). Default is 0 <br />
--major2  : Second error elliptical major axis (in ARC-MINUTES). Default is 0 <br />
--minor2  : Second error elliptical minor axis (in ARC-MINUTES). Default is 0 <br />
--angle2  : Position angle of the second error elliptical (in DEGREES). Default is 0



## Examples

```bash
$ ./bin/vou-blazars --ra 153.76 --dec 49.43 --FOV 30
(Simplest input)

$ ./bin/vou-blazars --ra 153.76 --dec 49.43 --FOV 30 --nh 3.e21 --mode s
(Specify nh and running mode)

$ ./bin/vou-blazars --ra 153.76 --dec 49.43 --FOV 30 --radius 15 --radius2 7
(Specify 2 error circles, one radius 15, the other radius is 7)

$ ./bin/vou-blazars --ra 153.76 --dec 49.43 --FOV 30. --mode s --radius 15
(Specify one error circle radius and running mode)

$ ./bin/vou-blazars --ra 153.76 --dec 49.43 --FOV 30. --nh 3.e21 -- major 15 --minor 10 --angle 120
(Specify the nh, and the error elliptical, major axis, minor axis, position angle)

$ ./bin/vou-blazars --ra 153.76 --dec 49.43 --FOV 30. -- radius 20 -- major2 15 --minor2 10 --angle2 120
$ ./bin/vou-blazars --ra 153.76 --dec 49.43 --FOV 30. -- major 35 --minor 10 --angle 120 -- radius2 20 
(Specify one error cirle and one error elliptical)

$ ./bin/vou-blazars --ra 153.76 --dec 49.43 --FOV 30. --mode l --major 15 --minor 10 --angle 120 --major2 5 --minor2 2 --angle2 90
(Specify two error ellipticals and running mode)
```

