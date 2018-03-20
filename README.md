# VOU_Blazars

VOU-Blazars is an extremely alpha version of blazars search tool.

## Install

### Dependencies

See below for dependencies setup.

* [PGPlot](http://www.astro.caltech.edu/~tjp/pgplot/)
* GFortran
* Python 2.7
* [EADA and conesearch pipline](https://github.com/chbrandt/eada)

**Once** the dependencies have been satisfied, go inside `VOU_Blazars/bin/fort` and run the `compile.sh` script:

```bash
$ git clone https://github.com/ecylchang/VOU_Blazars.git
$ cd VOU_Blazars
$ cd bin/fort
$ ./compile.sh
$ cd ../..
```

### Setup EADA/Conesearch pipeline

If you use [Anaconda Python Distribution](https://www.anaconda.com/download/), you can install EADA in its own (virtual) environment:

```bash
$ conda create -n eada2 python=2 pip
$ source activate eada2
$ pip install astropy
$ pip install pyvo
$ pip install pyyaml
$
$ git clone https://github.com/chbrandt/eada
$ cd eada
$ pip install .
```

### Setup PGPlot - MACOS

MacOS users are advised to use [brew](https://brew.sh/)
 
```bash
$ brew install pgplot
```


## Running VOU-Blazars

> Guarantee EADA's `conesearch` tool is available in your environment -- check this note's section [Setup EADA/Conesearch pipeline](#Setup-EADA/Conesearch-pipeline).

To run the tool,  RA, DEC, SEARCHING RADIUS are mandatory to input., and the other parameters, nH and error regions are optionals. The input could be RA, Dec, Radius, nH (if specified), error radius (if specified)

```
$ cd VOU_Blazars
$ ./bin/vou-blazars RA DEC RADIUS
```

The optional error radius could be :
*one circle or two circle. (error radius)
*one elliptical or two elliptical (major axis, minor axis, angle)
*One circle and one elliptical (circle first)
*nothing

The R.A. Dec. are in degrees, and searching radius is in arcmin.
nH is in cm^2. The radius of the error circles and the axes of the error ellipticals are in arcmin. The position angle is in degree, north-east on sky.

The nh are set to 5.e20 if no specify and no Heasoft installed. (If you have installed the Heasoft, and did not specify the nh, it will use the value calculated by Heasoft). The error radius are set to 0. If not specify.

## output examples
```bash
$ ./bin/vou-blazars 153.76 49.43 30.

$ ./bin/vou-blazars 153.76 49.43 30. 3.e21

$ ./bin/vou-blazars 153.76 49.43 30. 3.e21 15 7
(With 2 error circles, one radius 15, the other radius is 7)

$ ./bin/vou-blazars 153.76 49.43 30. 15 7
(Same as above, just no specify nh)

$ ./bin/vou-blazars 153.76 49.43 30. 15
(Only specify one error circle radius, and the other is set to 0.)

$ ./bin/vou-blazars 153.76 49.43 30. 3.e21 15 10 120
(Specify the nh, and the error elliptical, major axis, minor axis, position angle)

$ ./bin/vou-blazars 153.76 49.43 30. 20 15 10 120
(Specify two error regions, one is circle with radius 20, the other is elliptical 15 10 120degree)
Note that if you specify one circle and one elliptical, circle go first.

$ ./bin/vou-blazars 153.76 49.43 30. 3.e21 20 15 10 120
(Same as above, but also specify the nh value)

$ ./bin/vou-blazars 153.76 49.43 30. 3.e21 15 10 120 5 2 90
(Specify two error ellipticals)
```

