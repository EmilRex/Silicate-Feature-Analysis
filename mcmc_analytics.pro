pro mcmc_analytics,object_name

COMMON file_path, in_dir, out_dir, fit_name

for i=1,50 do begin
  ; Look at final chisq values for all links
  data = readfits(out_dir+'/'+object_name+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=i,/silent)
  last_link = -2*data[12,99,*]
  print,min(last_link)
  plot,findgen(20),last_link
  ;plothist,last_link
  stop
endfor


; Look at total length of each link

run = make_array(5000,1,/integer)
chisq = make_array(5000,1,/float)


for i=0, 4999 do begin
  run[i]= (i+1)
endfor

for k=0,19 do begin
  for i=1,50 do begin
    data = readfits(out_dir+'/'+object_name+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=i,/silent)
    for j=1,100 do begin
      chisq[(j-1)+((i-1)*100)] = -2*data[12,(j-1),k]
    endfor
  endfor
  plot,run,chisq,title='Chain: '+string((k+1))
  print,min(chisq)
  stop
endfor




end