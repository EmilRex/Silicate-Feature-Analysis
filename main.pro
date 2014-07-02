pro main

; *************************************************** ;
; Perform setup operations

; Home directory
home_dir = '~/Summer2014/Silicate_Feature_Analysis/'
CD, home_dir

; Path
!PATH=!PATH+':'+Expand_Path('+'+home_dir)

; Functional directories
COMMON file_path, in_dir, out_dir
in_dir = 'savfiles_MIPS_SED_corrected'
out_dir = 'output_v2'

; *************************************************** ;
; Select object and fit model
name = 'HD117214'
fittype = 'multi_mips' ;'multi_mips'

; *************************************************** ;
; Run program
fits_v1, name=name, fittype=fittype

; Print plot of result
plot_multi, name

; *************************************************** ;
; Print results

mcmc_result = readfits(out_dir+'/'+name+'_chn_mcmc_multi_part.fits',EXTEN_NO=51)

print, 'Parameter results: '
print, ' '
print, mcmc_result[0:(n_elements(mcmc_result)-2)]
print, ' '
print, 'Chisq: ', -2.0*mcmc_result[12]

end