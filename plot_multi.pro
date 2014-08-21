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

COMMON file_path, in_dir, out_dir, fit_name

; -------------------------------------------------- ;
; Replace with new structure for reading in data. 
; Take from _chn_mcmc_multi_part.fits, extension = final output
; generated through mcmc_m_mips_v1

; Get result of simulation to plot model
mcmc_result = readfits(out_dir+'/'+object_name+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=51)
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
restore,in_dir+'/'+object_name+'.sav'
      
wave_irs = final_wave ;[final_wave,71.42]
fl_diff = final_spec ;[final_spec,mips70_val]
uncer_irs = final_specerr ;[final_specerr,mips70_error]

; *************************************************** ;
; Separate out data types and set colors/markers

IRS_wave = wave_irs[where(strmatch(final_source, 'SpitzerIRS') EQ 1)]
IRS_spec = fl_diff[where(strmatch(final_source, 'SpitzerIRS') EQ 1)]
IRS_color = 1 ; red
IRS_psym = 1 ; plus

;MIPS70_wave = wave_irs[where(strmatch(final_source, 'Spitzer_MIPS_Phot') EQ 1)]
;MIPS70_spec = fl_diff[where(strmatch(final_source, 'Spitzer_MIPS_Phot') EQ 1)]
;MIPS70_color = 2 ; red
;MIPS70_psym = 6 ; square

MIPS_SED_wave = wave_irs[where(strmatch(final_source, 'SpitzerMIPS_SED') EQ 1)]
MIPS_SED_spec = fl_diff[where(strmatch(final_source, 'SpitzerMIPS_SED') EQ 1)]
MIPS_SED_color = 3 ; blue
MIPS_SED_psym = 5 ; triangle

; Calculate model spectrum
link[2]=10^(link[2])
link[8]=10^(link[8])
lines = [1,2,3,4]

x_start = min(wave_irs)
x_range = max(wave_irs)-min(wave_irs)
model_x = (findgen(round(x_range)*100)/100 + x_start)
out_model = modeltwograin(model_x, link)
tmp1 = round(chisq_best*100.)/100.

; *************************************************** ;
; Begin plotting

; Set up device
set_plot,'PS'
device, filename ='plots/MIPS_SED_'+object_name+'_'+fit_name+'.ps',/COLOR,/HELVETICA,XSIZE=15,YSIZE=12.5 & !p.font =0
loadct,39
!p.background=16777215

;Make IRS Spectrum Plot w/ Errorbars                                                                                                    
TVLCT,[0,255,0,0],[0,0,255,0],[0,0,0,255]

; Plot data points
plot,wave_irs,fl_diff,title=object_name+' (Two-Grain Model)', $
     ystyle=1,psym=4,xstyle=1,xtitle='Wavelength ('+cggreek('mu')+'m)', $
     ytitle='F'+cggreek('nu')+' (Jy)',charthick=1, thick=1, $
     xthick=2, ythick=2, charsize=1,color=0, $
     yrange=[min([out_model,fl_diff]),max([out_model,fl_diff])]

; Connect the points
oplot,wave_irs,fl_diff,color=0

; Overlay data with different colors and markers
oplot,IRS_wave,IRS_spec,color=IRS_color,psym=IRS_psym
;oplot,MIPS70_wave,MIPS70_spec,color=MIPS70_color,psym=MIPS70_psym
oplot,MIPS_SED_wave,MIPS_SED_spec,color=MIPS_SED_color,psym=MIPS_SED_psym

; Add error bars
oploterr,wave_irs,fl_diff,uncer_irs;,psym=1;,color=0

; Plot the model
oplot,model_x,out_model,color=1, thick=5,linestyle=lines[1]

; Create legend
;LEGEND,TEXTOIDL('\chi^{2}/d.o.f. : ')+  strtrim(string(tmp1,format='(f18.2)'),1), /top, /left,color=1,linestyle=lines[1]

device,/close
set_plot,'x'

end