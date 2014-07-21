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

  COMMON disk_benchmarking, run, times, lines
  run = 0
  times = dblarr(8000,11)
  lines = dblarr(8000,8)
  
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
      null = display_historic(name)
    ENDFOREACH
  ENDFOREACH

end