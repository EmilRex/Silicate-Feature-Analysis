pro main

; *************************************************** ;
; Perform setup operations

; Home directory
home_dir = '~/Summer2014/Silicate_Feature_Analysis/'
CD, home_dir

; Path
!PATH=!PATH+':'+Expand_Path('+'+home_dir)

; Functional directories
COMMON file_path, in_dir, out_dir, fit_name
in_dir = 'savfiles_MIPS_SED_corrected'
out_dir = '../Silicate_Feature_Analysis_output'

; *************************************************** ;
; Select object and fit model
; Remember to comment out bad new data points for new runs
;name = 'HD117214' 
name = 'HD106906'

fit_name = 'multi_mips' ;

; *************************************************** ;
; Run program
fits_v1, name=name, fittype=fit_name

; Print plot of result
plot_multi, name

; See how well mcmc did
;mcmc_analytics, name

; *************************************************** ;
; Print results

mcmc_result = readfits(out_dir+'/'+name+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=51)

print, 'Parameter results: '
print, ' '
print, mcmc_result[0:(n_elements(mcmc_result)-2)]
print, ' '
print, 'Chisq: ', -2.0*mcmc_result[12]

end