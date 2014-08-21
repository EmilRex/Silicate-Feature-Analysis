; +
; NAME:
;  plot_result.pro
;
; PURPOSE:
;  Plot a spectrum of the given object with the fitted model overlayed and save as post script to /plots
;
; INPUTS:
;   NAME: Name of object
;
; KEYWORDS
;   SEPARATE: Plot flux from each grain species separately
;
; OUTPUTS:
;   .ps file with plot
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
;  Generalized and renamed by EC (7/14/14)
; -
; *************************************************** ;
pro plot_result, separate=separate

COMMON file_path, in_dir, out_dir, fit_name, object_name

; Check if data file exists. If not exit program.
if (file_test(out_dir+'/'+object_name+'_chn_mcmc_'+fit_name+'_part.fits') eq 0) then begin
  print, "Data not found for object: "+object_name
  print, "Continuing to next object..."
  return
endif

; Get result of simulation to plot model
mcmc_result = readfits(out_dir+'/'+object_name+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=51,/silent)
link=mcmc_result[0:(n_elements(mcmc_result)-2)]
chisq_best = mcmc_result[(n_elements(mcmc_result)-1)]

; Get Teff, amin and dist_val for the object
fmt='a,f,f,f'
readcol,'input_files/input_param_file.txt',F=fmt,catalog_nameA,c_teff,c_amin,c_dist_val,/silent
  
for i = 0, size(catalog_nameA,/n_elements)-1 do begin
  if (catalog_nameA[i] eq object_name) then begin
    Teff=c_teff[i]
    amin=c_amin[i]
    dist_val=c_dist_val[i]
  endif
endfor


; Deal with global values
COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, Qwaterice, crystallineabs
COMMON stellarprops, temptable, folivine, effectiveTemp, lambdastar, fluxstar
COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, waterice_emit, effectiveTempArray, stellar_emit
restore, 'graintempdata.sav'
restore, 'qtables_withcrys2.sav' ; qastrosil, qolivine, qpyroxene     
restore, 'qwaterice.sav'

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

                                              

restore, grainfiles[ki]
restore,in_dir+'/'+object_name+'.sav'
      
wave_irs = final_wave
fl_diff = final_spec
uncer_irs = final_specerr

x_start = min(wave_irs)
x_range = max(wave_irs)-min(wave_irs)
model_x = (findgen(round(x_range)*100)/100 + x_start)
lines = [1,2,3,4]

; *************************************************** ;
; Fittype specifics
; Calculate model spectrum

; Define constants
m_moon = 7.34767309e22 ; in g, from google
r_sun = 0.00464913034 ;AU

plot_old = 0

IF (fit_name eq 'single') THEN BEGIN
  
  out_model = modelsinglespectrum(transpose(model_x), link, /single)
  
  if keyword_set(separate) then begin
    out_model_separate = modelsinglespectrum(transpose(model_x), link, /single, /separate)
    n_models = 0
  endif
  
ENDIF



IF (fit_name eq 'multi') THEN BEGIN

  out_model = modelsinglespectrum(transpose(model_x), link, /multi)

  if keyword_set(separate) then begin
    out_model_separate = modelsinglespectrum(transpose(model_x), link, /multi, /separate)
    n_models = 1
  endif

ENDIF 



IF (fit_name eq 'disk') THEN BEGIN
  
  out_model = modelsinglespectrum(transpose(model_x), link, /disk)

  if keyword_set(separate) then begin
    out_model_separate = modelsinglespectrum(transpose(model_x), link, /disk, /separate)
    n_models = 0
  endif

ENDIF

tmp1 = round(chisq_best*100.)/100.

; *************************************************** ;
; Begin plotting

sep_label = ""
if keyword_set(separate) then begin
  sep_label = "_separate"
endif

; Set up device
set_plot,'PS'
device, filename ='plots/'+object_name+'_'+fit_name+sep_label+'.ps',/COLOR,/HELVETICA,XSIZE=15,YSIZE=12.5 & !p.font =0
loadct,39,/silent
!p.background=16777215

;Make IRS Spectrum Plot w/ Errorbars                                                                                                    
TVLCT,[0,255,0,0],[0,0,255,0],[0,0,0,255]

; Plot data points
plot,wave_irs,fl_diff,title=object_name+' ('+fit_name+' Model)', $
     ystyle=1,psym=0,xstyle=1,xtitle='Wavelength ('+cggreek('mu')+'m)', $
     ytitle='F'+cggreek('nu')+' (Jy)',charthick=1, thick=1, $
     xthick=2, ythick=2, charsize=1,color=0;, $
     ;yrange=[1.0e-6,1.0],/ylog

; Add error bars
oploterr,wave_irs,fl_diff,uncer_irs,0;,psym=1;,color=0

; Plot the models
oplot,model_x,out_model,color=3, thick=5,linestyle=lines[1]

legend_names = ['Model']
legend_psyms = [0]
legend_colors=[3]
legend_linestyle=[2]
legend_textcolors=[0]

if keyword_set(separate) then begin
  for j=0,n_models do begin
    for i=0,4 do begin
      oplot,model_x,out_model_separate[*,(j*5+i)],color=(1+j),linestyle=(i)
      legend_psyms = [legend_psyms,0]
      legend_colors = [legend_colors,1+j]
      legend_linestyle = [legend_linestyle,i]
      legend_textcolors = [legend_textcolors,0]
    endfor
    legend_names = [legend_names,'Olivine','Pyroxene','Forsterite','Enstatite','Water-ice']
  endfor
endif

; Create legend
legend,[legend_names],psym=legend_psyms,colors=legend_colors,linestyle=legend_linestyle,textcolors=legend_textcolors;,corners=[150.0,0.8,200.0,1.6]

xyouts,35,0.1,cggreek('chi')+'!E2!N!X / d.o.f. : '+strtrim(string(-2.0*chisq_best,format='(f18.2)'),1),/data


device,/close
set_plot,'x'

end