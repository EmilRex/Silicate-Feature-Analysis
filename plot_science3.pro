; *************************************************** ;
function display_results,name

  COMMON file_path, in_dir, out_dir, fit_name
  mcmc_result = readfits(out_dir+'/'+name+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=51,/silent)
  
  print, ' '
  print, systime()
  print, 'Parameter results for: '
  print, [name,fit_name]
  print, mcmc_result[0:(n_elements(mcmc_result)-2)]
  print, ' '
  print, 'Chisq: ', -2.0*mcmc_result[(n_elements(mcmc_result)-1)]
  print, ' '
  
end
; *************************************************** ;


pro plot_science3

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
  out_dir = '../Science3_output'
  
  ; *************************************************** ;
  ;RUN MULTIPLE
  ; *************************************************** ;
  
  names = ['HD146897','HD117214']
  fit_names = ['single','multi_mips', 'disk_mips']
  
  FOREACH name, names DO BEGIN
    FOREACH fit_name, fit_names DO BEGIN
      ;fits_v1, name=name, fittype=fit_name
      plot_result, name
      null = display_results(name)
    ENDFOREACH
  ENDFOREACH

end