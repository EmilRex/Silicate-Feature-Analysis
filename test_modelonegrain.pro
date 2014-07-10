pro test_modelonegrain
; Data from most recent run of HD117214


; Set working directory
home_dir = '~/Summer2014/Silicate_Feature_Analysis/'
CD, home_dir

; Set path for auxillary files
!PATH=!PATH+':'+Expand_Path('+'+home_dir)

; Define file_path common structure
COMMON file_path, in_dir, out_dir, fit_name
in_dir = 'savfiles_MIPS_SED_corrected'
out_dir = 'output_v2'
name = 'HD117214'
fit_name = 'multi_mips' ;

; *************************************************** ;

plot_multi,name
stop
mcmc_result = readfits(out_dir+'/'+name+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=51)
params1=mcmc_result[0:5]
params2=mcmc_result[6:11]

mcmc_data = readfits(out_dir+'/'+name+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=0)
lambda = mcmc_data[0,*]
flux = mcmc_data[1,*]

;plot, lambda, flux

result1 = modelonegrain(lambda, params1, mie=mie)
result2 = modelonegrain(lambda, params2, mie=mie)

plot, lambda, (result1+result2)

end