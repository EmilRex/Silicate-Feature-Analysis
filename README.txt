README.txt
8/21/2014

AUTHORS:
Tushar Mittal - tmittal2@berkley.edu
Hannah Jang-Condell -  hjangcon@uwyo.edu
Christine Chen - cchen@stsci.edu
Emil Christensen - chris2er@dukes.jmu.edu

DESCRIPTION:
Silicate_Feature_Analysis contains code for fitting silicates and water-ice features of 
circumstellar dust. The code is based on work done for the following papers:

1) “The IRS Debris Disk Catalog: Spectral Feature Analysis” Tushar Mittal, Christine H. 
	Chen, Hannah Jang-Condell, P. Manoj, Benjamin Sargent, Dan Watson, William Forrest, & 
	Carey Lisse 2013 ApJS, in preparation
2) “The IRS Debris Disk Catalog: Continuum Analysis” Christine H. Chen, Tushar Mittal, 
	Marc Kuchner, William Forrest, Carey Lisse, P. Manoj, Benjamin Sargent, & Dan Watson 
	2013, ApJS, in preparation

The code is still in a pretty raw form and may need a good deal of trouble shooting. We 
welcome help with this work. Please contact Emil Christensen or Christine Chen if you wish
to build on this project.

USAGE:
1) Update data (See BetaPic for example)
	a) In input_files/input_param_file.txt, add the object's name, temperature (in Kelvin),
	   minimum grain size (in microns), and distance to the star in parsecs.
	b) Create an IDL savefile as the name of the object (name.sav) with the wavelengths
	   in microns as FINAL_WAVE, the flux in Janskies as FINAL_SPEC, the error of the flux 
	   in Janskies and the name of the source of the data point.
2) Set parameters in main.pro
	a) Set the "home_dir" to the directory containing the repository
	b) Set the names of the objects you wish to simulate along with the types of fits
3) Run! 
	   
DISCLAIMER & Acknowledgements;
This work is being conducted by Emil Christensen at the Space Telescope Science Institute 
(STScI) under the supervision of Dr. Christine Chen. Much of the code was developed by 
Tushar Mittal at STScI and Dr. Hannah Jang-Condell at the University of Wyoming. 

The use of this work for commercial reasons is strictly prohibited. Please contact Emil
Christensen at the email above before using this work.

Copyright 2014 Space Telescope Science Institute 

VERSION LOG:
1.0.0 - First generalized version