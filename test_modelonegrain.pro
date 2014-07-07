pro test_modelonegrain

; Data from most recent run of HD117214
; collected through python astropy.io

; Home directory
home_dir = '~/Summer2014/Silicate_Feature_Analysis/'
CD, home_dir

; Path
!PATH=!PATH+':'+Expand_Path('+'+home_dir)




mcmc_result = readfits('output_v2/HD117214_chn_mcmc_multi_mips_part.fits',EXTEN_NO=51)
params1=mcmc_result[0:5]
params2=mcmc_result[6:11]

mcmc_data = readfits('output_v2/HD117214_chn_mcmc_multi_mips_part.fits',EXTEN_NO=0)
lambda = mcmc_data[0,*]
flux = mcmc_data[1,*]

;plot, lambda, flux

result1 = modelonegrain(lambda, params1, mie=mie)
result2 = modelonegrain(lambda, params2, mie=mie)

plot, lambda, (result1+result2)

end