pro param_hist
; Plot histograms of rlaw and/or alaw

out_dir = '../Science3_output';/old'
names = ['HD95086','HD106906','HD108257','HD110058','HD111520','HD113556','HD113766','HD114082','HD115600','HD117214','HD145560','HD146181','HD146897']
;fit_name = "disk_mips"
fit_name = "multi_mips"

k = 0

for i=0,(n_elements(names)-1) do begin
  
  ; Skip object if no data exists
  if (file_test(out_dir+'/'+names[i]+'_chn_mcmc_'+fit_name+'_part.fits') eq 0) then begin
    print, "No data for "+names[i]+" ("+fit_name+")"
    continue
  endif
  
  mcmc_result = readfits(out_dir+'/'+names[i]+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=51,/silent)
  param_current = mcmc_result[2]
  
  if (k eq 0) then begin
    param = param_current
  endif
  if (k gt 0) then begin
    param = [param,param_current]
  endif
  
  k = k + 1
endfor

; Plot histogram of data
plothist, param,/autobin

print, param

end