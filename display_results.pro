; Separated from main.pro by EC on 7/21/14
; ****************************************************************************************************** ;
function display_results,name
  ; ****************************************************************************************************** ;
  ; Print simulation results to CL
  
  COMMON file_path, in_dir, out_dir, fit_name, object_name
  mcmc_result = readfits(out_dir+'/'+name+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=51,/silent)
  
  ; Load R-star values
  fmt = 'a,f'
  readcol, 'r_star.csv',F=fmt,r_star_name,c_r_star,/silent
  
  ; Load correct r_star
  for k=0, n_elements(r_star_name)-1 do begin
    if (r_star_name[k] eq name) then begin
      r_star = c_r_star[k]
      break
    endif
  endfor
  
  r_sun = 0.00464913034 ;AU
  
  print, ' -------------------------------------------------------------------------------------------- '
  print, systime()
  print, name,' --- ',fit_name,' --- NEW'
  print, ' '
  CASE fit_name OF
    'single': BEGIN
      print,['Dist[AU] ','a_grain ','dustmass ','folive ','fcrys ','ffors ']
      
      dist = (r_sun*r_star)*(10^mcmc_result[0])
      print, [dist, mcmc_result[1:5]]
      
    END
    'multi_mips': BEGIN
      print,['Dist1[AU] ','a_grain1 ','dustmass1 ','folive1 ','fcrys1 ','ffors1 ','Dist2[AU] ','a_grain2 ','dustmass2 ','folive2 ','fcrys2 ','ffors2 ']
      
      dist1 = (r_sun*r_star)*(10^mcmc_result[0])
      dist2 = (r_sun*r_star)*(10^mcmc_result[6])
      print, [dist1, mcmc_result[1:5], dist2, mcmc_result[7:11]]
      
    END
    'disk_mips': BEGIN
      print, ['rin ','rout ','rlaw ','amin ','amax ','alaw ','diskmass ','folive ','fcrys ','ffors ']
      
      print, mcmc_result[0:(n_elements(mcmc_result)-2)]
      
    END
  ENDCASE
  
  print, ' '
  print, 'Chisq: ', -2.0*mcmc_result[(n_elements(mcmc_result)-1)]
  print, ' '
  
end