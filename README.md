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
$ gfortran -c mylib.f -ffixed-line-length-132
```
and then all the other programs (with .f file name) with command
```bash
$ gfortran -o $XXX $XXX.f -ffixed-line-length-132 mylib.o -L${HOME}/pgplot -lpgplot
```
**Change the PGPLOT path to where you put the PGPLOT library.**

**Always compile all the Fortran routines/subroutines inside the folder "fort"**


## Running VOU-Blazars

> Guarantee EADA's `conesearch` tool is available in your environment -- check this note's section [Setup EADA/Conesearch pipeline](#setup-eadaconesearch-pipeline).

To run the tool, `RA`, `DEC`, `SEARCH_RADIUS` are mandatory to input., and the other parameters, nH and error regions are optionals.
The input parameters are sorted by R.A., Dec., Searching Radius, nH ( if specified ), and error circles or/and ellipticals ( if specified ).

**Before running the tool, make sure that you are under the main folder where the VOU-Blazars are installed and have activated the eada environment** Bin folder should inside the VOU_Blazars folder.


```
$ cd VOU_Blazars
$ conda activate eada
$ ./bin/vou-blazars RA DEC SEARCH_RADIUS
```

The optional error radius could be :
* one circle or two circle. (error radius)
* one elliptical or two elliptical (major axis, minor axis, angle)
* one circle and one elliptical (circle first)
* nothing

The R.A. Dec. are in degrees, and searching radius is in arcmin.
nH is in cm^2.The radius of the error circles and the axes of the error ellipticals are in arcmin. The position angle is in degree, north-east on sky.
The nh are set to 5.e20 if no specify and no Heasoft installed. (If you have installed the Heasoft, and did not specify the nh, it will use the value calculated by Heasoft). The error radius are set to 0. If not specify.
The last cardinal thing to mention is that if two error circles or ellipticals are specified, the one with a larger size should input first.


## Examples

```bash
$ ./bin/vou-blazars 153.76 49.43 30.

$ ./bin/vou-blazars 153.76 49.43 30. 3.e21

$ ./bin/vou-blazars 153.76 49.43 30. 3.e21 15 7
(With 2 error circles, one radius 15, the other radius is 7)

$ ./bin/vou-blazars 153.76 49.43 30. 15 7
(Same as above, just no nh specified)

$ ./bin/vou-blazars 153.76 49.43 30. 15
(Only specify one error circle radius, and the other is set to 0.)

$ ./bin/vou-blazars 153.76 49.43 30. 3.e21 15 10 120
(Specify the nh, and the error elliptical, major axis, minor axis, position angle)

$ ./bin/vou-blazars 153.76 49.43 30. 20 15 10 120
(Specify two error regions, one is circle with radius 20, the other is elliptical 15 10 120 degrees)
Note that if you specify one circle and one elliptical, circle go first.

$ ./bin/vou-blazars 153.76 49.43 30. 3.e21 20 15 10 120
(Same as above, but also specify the nh value)

$ ./bin/vou-blazars 153.76 49.43 30. 3.e21 15 10 120 5 2 90
(Specify two error ellipticals)
```

