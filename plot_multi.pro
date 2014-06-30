; +
; NAME:
;  plot_multi.pro
;
; PURPOSE:
;  Plot a spectrum of the given object with the best fitting model overlayed and save as post script to /plots
;
; INPUTS:
;   NAME: Name of object
;
; KEYWORDS
;   NONE
;
; OUTPUTS:
;   post script plot
;
; AUTHORS:
;  Tushar Mittal
;  Emil Christensen - chris2er@dukes.jmu.edu
;  Christine Chen
;
; DISCLAIMER
;  This software is provided as is without any warranty whatsoever.
;  Permission to use, copy, modify, and distribute modified or
;  unmodified copies is granted, provided this disclaimer
;  is included unchanged.
;
; MODIFICATION HISTORY:
;  Written by TM (June 2013) 
;  Modified by EC (6/27/2014)
; -
; *************************************************** ;
pro plot_multi, object_name

; -------------------------------------------------- ;
; Replace with new structure for reading in data. 
; Take from _chn_mcmc_multi_part.fits, extension = final output
; generated through mcmc_m_mips_v1

; Get result of simulation to plot model
mcmc_result = readfits('output_v1/'+object_name+'_chn_mcmc_multi_part.fits',EXTEN_NO=51)
link=mcmc_result[0:11]
chisq_best = mcmc_result[12]

print, mcmc_result

; Get Teff, amin and dist_val for the object
fmt='a,f,f,f'
readcol,'input_files/input_param_file.txt',F=fmt,catalog_nameA,c_teff,c_amin,c_dist_val
  
 
for i = 0, size(catalog_nameA,/n_elements)-1 do begin
  if (catalog_nameA[i] eq object_name) then begin
    Teff=c_teff[i]
    amin=c_amin[i]
    dist_val=c_dist_val[i]
  endif
endfor


; Deal with global values
COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs
COMMON stellarprops, temptable, folivine, effectiveTemp, lambdastar, fluxstar
COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, effectiveTempArray, stellar_emit
restore, 'graintempdata.sav'
restore, 'qtables_withcrys2.sav' ; qastrosil, qolivine, qpyroxene     

; *************************************************** ;                           
; Retrieve grainprops

effectiveTemp = Teff

; Calculate based on masses based on Isochrones etc
; find the right grain model for calculating temperatures             
cmd = 'ls modelgrids/Teff*grains.sav'
spawn, cmd, grainfiles
  
; read in temperatures                                                
strbeg = strpos(grainfiles, 'Teff')+3
strend = strpos(grainfiles, 'grains')

tarray = fltarr(n_elements(grainfiles))
for i=0,n_elements(tarray)-1 do begin
  tarray[i] = float(strmid(grainfiles[i], strbeg[i]+1, strend[i]-strbeg[i]-1))
endfor
  
; reorder arrays                                                      
ii = sort(tarray)
grainfiles = grainfiles[ii]
tarray = tarray[ii]
        
kuruczindex = interpol(findgen(n_elements(tarray)),tarray,effectiveTemp)
ki = round(kuruczindex) < (n_elements(tarray)-1) > 0

; next command will retrieve temptable, folivine, Teff                
; as generated in generateallgrid.pro                                 
; folivine = [0.0, 1.0]                                               

restore, grainfiles[ki]
restore,'old_savfiles_mcmc/'+object_name+'.sav'
      
wave_irs = [final_wave,71.42]
fl_diff = [final_spec,mips70_val]
uncer_irs = [final_specerr,mips70_error]

; *************************************************** ;
; Begin plotting

; Set up device
set_plot,'PS'
device, filename ='plots/'+object_name+'_multi.ps',/COLOR,/HELVETICA,XSIZE=15,YSIZE=12.5 & !p.font =0
loadct,39
!p.background=16777215

;Make IRS Spectrum Plot w/ Errorbars                                                                                                    
TVLCT,[0,255,0,0],[0,0,255,0],[0,0,0,255]

; Plot data points
plot,wave_irs,fl_diff,title=object_name+' (Two-Grain Model)', $
     ystyle=1,psym=4,xstyle=1,xtitle='Wavelength ('+cggreek('mu')+'m)', $
     ytitle='F'+cggreek('nu')+' (Jy)',charthick=1, thick=1, $
     xthick=2, ythick=2, charsize=1,color=0

; Connect the points
oplot,wave_irs,fl_diff,color=0

; Add error bars
oploterr,wave_irs,fl_diff,uncer_irs;,psym=1;,color=0

; Calculate model spectrum
link[2]=10^(link[2])
link[8]=10^(link[8])
lines = [1,2,3,4]
out_model = modeltwograin(wave_irs, link)
tmp1 = round(chisq_best*100.)/100.

; Plot the model
oplot,wave_irs,out_model,color=1, thick=5,linestyle=lines[1]

; Create legend
;LEGEND,TEXTOIDL('\chi^{2}/d.o.f. : ')+  strtrim(string(tmp1,format='(f18.2)'),1), /top, /left,color=1,linestyle=lines[1]

device,/close
set_plot,'x'

end