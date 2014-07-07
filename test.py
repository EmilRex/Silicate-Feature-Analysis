import astropy.io.fits as fits
import matplotlib.pyplot as plt 

hdu = fits.open('output_v2/HD117214_chn_mcmc_multi_mips_part.fits')

data = hdu[0].data


#plt.plot(data(1,),data(2,))