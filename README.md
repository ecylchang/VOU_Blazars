# VOU_Blazars

To run the tool, 
Gfortran comopiler, PGPlot, Python, are need to be installed.
Anaconda and conesearch pipline also are need.

All the files should be put to the same folder. The requried files are ;
vou-blazars, vueue.sh
find_candidates1, find_candidates2, find_candidates_int, readcat, plot_sed, gnomo_plot_types, connvert_sed.f
mylib.f nhdeabsorb2.f
cats1.ini, cats2.ini, count_rate folder, status.code

############CONESEARCH
Eada environment and conesearch pipline (need anaconda and python)
See setup eada to install the conesearch and set up the eada environment.
Also, eada environment should be activated before each runing.

############COPPILE FORTRAN PROGRAMS
After set up the conesearch and eada, next is to compile all the fortran programs.
First is to compile the mylib.f and nhdeabsorb2.f into object file, using commend : gfortran -c XXX.f -ffixed-line-length-132
then compile all the others program (find_candidates1, find_candidates2, find_candidates_int, readcat, plot_sed, gnomo_plot_types, connvert_sed.f),
using command : source compile XXX

##########To RUN THE Program
#To run the tool,  RA, DEC, SEARCHING RADIUS are mandatory to input., and the other parameters, nH and error regions are optionals.
#The input could be RA, Dec, Radius, nH (if specified), error radius (if specified)

#The error radius could be :
#-one circle or two circle. (error radius)
#-one elliptical or two elliptical (major axis, minor axis, angle)
#-One circle and one elliptical (circle first)
#-nothing

#The R.A. Dec. are in degrees, and searching radius is in arcmin.
#nH is in cm^2
#The radius of the error circles and the axes of the error ellipticals are in arcmin.
#The position angle is in degree, north-east on sky.

#The nh are set to 5.e20 if no specify and no Heasoft installed. (If you have installed the Heasoft, and did not specify the nh, it will use the value calculated by Heasoft)
#The error radius are set to 0. If not specify.

#output EXAMPLES
#e.g.: source vou-blazars 153.76 49.43 30.
#e.g.: source vou-blazars 153.76 49.43 30. 3.e21
#e.g.: source vou-blazars 153.76 49.43 30. 3.e21 15 7 (With 2 error circles, one radius 15, the other radius is 7)
#e.g.: source vou-blazars 153.76 49.43 30. 15 7 (Same as above, just no specify nh)
#e.g.: source vou-blazars 153.76 49.43 30. 15 (Only specify one error circle radius, and the other is set to 0.)
#e.g.: source vou-blazars 153.76 49.43 30. 3.e21 15 10 120 (Specify the nh, and the error elliptical, major axis, minor axis, position angle)
#e.g.: source vou-blazars 153.76 49.43 30. 20 15 10 120(Specify two error regions, one is circle with radius 20, the other is elliptical 15 10 120degree)
#Note that if you specify one circle and one elliptical, circle go first.
#e.g. source vou-blazars 153.76 49.43 30. 3.e21 20 15 10 120(Same as above, but also specify the nh value)
#e.g. source vou-blazars 153.76 49.43 30. 3.e21 15 10 120 5 2 90(Specify two error ellipticals)


