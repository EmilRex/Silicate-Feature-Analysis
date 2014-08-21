README.txt
8/21/2014

This is the active research branch of the project.

Instructions for Dr. Chen:
1) To update water-ice q-values:
	a) put data file in the Silicate_Feature_analysis directory
		data should have the format: lwater_ice, m_re, m_im, format='(F,F,F)'
			i.e. [lambda,n,k]
	b) open qwater_ice.pro 
	c) on line 12 change "IOP_2008_ASCIItable.dat" to the new file
	d) run qwater_ice.pro
		This should update qwaterice.sav and output a bunch of lambda
		values to the command line. It takes a few minutes so this is
		a good time to go get a cup of coffee.
	e) once that's done, open up tempdata.pro and run that.
		It will output a lot of "grain sizes's... it will also
		take a while to run. 
		Not that for it to work, you need to have the kurucz/ directory 
		with all the kurucz models. You can take it off of vega. I think it 
		is in both Hannah's and Tushar's directories.
	f) Rejoice... that should have updated everything.
	
2) To run the models
	a) open up main.pro
	b) change the "home_dir" to the directory you cloned this into.
	c) change the other parameters as you see fit: object names, fits to
		be run, output locations etc.
	d) run main.pro
		everything should work fine. Make sure the output directory
		actually exists.


AUTHORS:
Tushar Mittal - tmittal2@berkley.edu
Hannah Jang-Condell -  hjangcon@uwyo.edu
Christine Chen - cchen@stsci.edu
Emil Christensen - chris2er@dukes.jmu.edu

DESCRIPTION:
Silicate_Feature_Analysis contains the library of IDL codes used for the analysis 
presented in "The Spitzer IRS Debris Disk Catalog: II. Silicate Feature Analysis" by 
Mittal et al. To begin please see main.pro


DISCLAIMER & Acknowledgements;
This work is still in development. Hopefully the repository will be ready by late August 
2014. This work is being conducted by Emil Christensen at the Space Telescope
Science Institute (STScI) under the supervision of Dr. Christine Chen. Much of the code 
was developed by Tushar Mittal at STScI and Dr. Hannah Jang-Condell at the University of
Wyoming. 

The use of this work for commercial reasons is strictly prohibited. Please contact Emil
Christensen at the email above before using this work.

DATA:
MIPS SED-mode data sourced from: 
http://dirty.as.arizona.edu/~kgordon/mips/cchen_mips_sed_point/cchen_mips_sed_point.html
