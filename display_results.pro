; Separated from main.pro by EC on 7/21/14
; ****************************************************************************************************** ;
function display_results,name
  ; ****************************************************************************************************** ;
  ; Print simulation results to CL
  
  COMMON file_path, in_dir, out_dir, fit_name, object_name
  mcmc_result = readfits(out_dir+'/'+name+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=51,/silent)
  
  print, ' '
  print, systime()
  print, 'New parameter results for: '
  print, [name,fit_name]
  CASE fit_name OF
    'single': BEGIN
      print,'Temp ','a_grain ','dustmass ','folive ','fcrys ','ffors '
    END
    'multi_mips': BEGIN
      print,'Temp1 ','a_grain1 ','dustmass1 ','folive1 ','fcrys1 ','ffors1 ','Temp2 ','a_grain2 ','dustmass2 ','folive2 ','fcrys2 ','ffors2 '
    END
    'disk_mips': BEGIN
      print, 'rin ','rout ','rlaw ','amin ','amax ','alaw ','diskmass ','folive ','fcrys ','ffors '
    END
  ENDCASE
  
  print, 'Temp','
  print, mcmc_result[0:(n_elements(mcmc_result)-2)]
  print, ' '
  print, 'Chisq: ', -2.0*mcmc_result[(n_elements(mcmc_result)-1)]
  print, ' '
  
end