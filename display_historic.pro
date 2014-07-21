; Separated from main.pro by EC on 7/21/14
; ****************************************************************************************************** ;
function display_historic,name
  ; ****************************************************************************************************** ;
  ; Print Tushar's old simulation results to CL
  
  COMMON file_path, in_dir, out_dir, fit_name
  
  ; To Do: add if statement to deal with non-existent historic data
  ; else print, 'No prior data available'
  
  CASE fit_name OF
    'single': BEGIN
      data_dir = '/Users/echristensen/Summer2014/dust_fit/hannah_model/pro/output_fin_new2/output_fin_new'
      old_result = readfits(data_dir+'/'+name+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=61,/silent)
      head = ['Temp ','a_grain ','dustmass ','folive ','fcrys ','ffors ']
    END
    
    'multi_mips': BEGIN
      data_dir = '/Users/echristensen/Summer2014/dust_fit/hannah_model/pro/output_fin_new2/output_fin_new'
      old_result = readfits(data_dir+'/'+name+'_chn_mcmc_'+'multi'+'_part.fits',EXTEN_NO=51,/silent)
      head = ['Temp1 ','a_grain1 ','dustmass1 ','folive1 ','fcrys1 ','ffors1 ','Temp2 ','a_grain2 ','dustmass2 ','folive2 ','fcrys2 ','ffors2 ']
    END
    
    'disk_mips': BEGIN
      data_dir = '/Users/echristensen/Summer2014/dust_fit/hannah_model/pro/output_fin_new2/output_fin_new'
      old_result = readfits(data_dir+'/'+name+'_chn_mcmc_'+'disk'+'_part.fits',EXTEN_NO=51,/silent)
      head = ['rin ','rout ','rlaw ','amin ','amax ','alaw ','diskmass ','folive ','fcrys ','ffors ']
    END
  ENDCASE
  
  print, ' '
  print, systime()
  print, 'Old parameter results for: '
  print, [name,fit_name]
  print, head
  print, old_result[0:(n_elements(old_result)-2)]
  print, ' '
  print, 'Chisq: ', -2.0*old_result[(n_elements(old_result)-1)]
  print, ' '
  
end