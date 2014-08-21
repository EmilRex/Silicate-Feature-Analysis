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
;  Generalized and renamed by EC (7/14/14)
; -
; *************************************************** ;
pro plot_result_PACS, separate=separate

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
IRS_psym = 2 ; asterisk

;MIPS70_wave = wave_irs[where(strmatch(final_source, 'Spitzer_MIPS_Phot') EQ 1)]
;MIPS70_spec = fl_diff[where(strmatch(final_source, 'Spitzer_MIPS_Phot') EQ 1)]
;MIPS70_color = 2 ; red
;MIPS70_psym = 6 ; square

PACS_wave = wave_irs[where(strmatch(final_source, 'Herschel_PACS') EQ 1)]
PACS_spec = fl_diff[where(strmatch(final_source, 'Herschel_PACS') EQ 1)]
PACS_color = 3 ; blue
PACS_psym = 6 ; square

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

  fmt = 'a,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f'
  readcol, 'multi_new.csv',F=fmt,db_name,chisq,temp1,temp2,Loc,Loc2,amin1,amin2,mass1,mass2,fcryst1,fcryst2,oliv1,oliv2,ffost1,ffost2,/silent

  FOR i=0,(n_elements(db_name)-1) DO BEGIN
    IF (object_name eq db_name[i]) THEN BEGIN

      plot_old = 1
      data2 = [temp1[i],amin1[i],(mass1[i]*m_moon),oliv1[i],fcryst1[i],ffost1[i],temp2[i],amin2[i],(mass2[i]*m_moon),oliv2[i],fcryst2[i],ffost2[i]]
      tushar_model = modeltwograin_old(model_x, data2)
      ;plot, model_x,tushar_model
      ;stop
      break
    ENDIF
  ENDFOR
ENDIF 



IF (fit_name eq 'disk') THEN BEGIN
  
  out_model = modelsinglespectrum(transpose(model_x), link, /disk)

  if keyword_set(separate) then begin
    out_model_separate = modelsinglespectrum(transpose(model_x), link, /disk, /separate)
    n_models = 0
  endif

  fmt = 'a,f,f,f,f,f,f,f,f,f,f,f'
  readcol, 'disk_new.csv',F=fmt,db_name,chisq,rin,rout,rlaw,amin,amax,alaw,diskmass,fcryst,foliv,ffost,/silent
  
  fmt = 'a,f'
  readcol, 'r_star.csv',F=fmt,r_star_name,c_r_star,/silent
  
  ; Load correct r_star
  for k=0, n_elements(r_star_name)-1 do begin
    if (r_star_name[k] eq object_name) then begin
      r_star = c_r_star[k]
      break
    endif
  endfor
  
  FOR i=0,(n_elements(db_name)-1) DO BEGIN
    IF (object_name eq db_name[i]) THEN BEGIN
      
      plot_old = 1
      data2 = [alog10(rin[i]/(r_sun*r_star)),alog(rout[i]/rin[i]),rlaw[i],amin[i],alog(amax[i])/amin[i],alaw[i],alog10((diskmass[i])*m_moon),foliv[i],fcryst[i],ffost[i]]
      tushar_model = modelsinglespectrum_old(transpose(model_x), data2 )
      
      break
    ENDIF
  ENDFOR
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
device, filename ='plots/PACS_'+object_name+'_'+fit_name+sep_label+'.ps',/COLOR,/HELVETICA,XSIZE=15,YSIZE=12.5 & !p.font =0
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

; Overlay data with different colors and markers
oplot,IRS_wave,IRS_spec,color=IRS_color,psym=IRS_psym
oplot,PACS_wave,PACS_spec,color=PACS_color,psym=PACS_psym

; Add error bars
oploterr,wave_irs,fl_diff,uncer_irs,0;,psym=1;,color=0

; Plot the models
oplot,model_x,out_model,color=3, thick=5,linestyle=lines[1]

legend_names = ['IRS','PACS','New Model','Old Model']
legend_psyms = [2,6,0,0]
legend_colors=[1,3,3,2]
legend_linestyle=[0,0,2,3]
legend_textcolors=[0,0,0,0]

if keyword_set(separate) then begin
  for j=0,n_models do begin
    for i=0,4 do begin
      oplot,model_x,out_model_separate[*,(j*5+i)],color=(i),linestyle=(1+j)
      legend_psyms = [legend_psyms,0]
      legend_colors = [legend_colors,i]
      legend_linestyle = [legend_linestyle,1+j]
      legend_textcolors = [legend_textcolors,0]
    endfor
    legend_names = [legend_names,'Olivine','Pyroxene','Forsterite','Enstatite','Water-ice']
  endfor
endif

if (plot_old eq 1) then begin  
  oplot,model_x,tushar_model,color=2, thick=5,linestyle=lines[2]
  print, "Old model plotted for: "+object_name
endif


; Create legend
legend,[legend_names],psym=legend_psyms,colors=legend_colors,linestyle=legend_linestyle,textcolors=legend_textcolors;,corners=[150.0,0.8,200.0,1.6]

xyouts,35,0.1,cggreek('chi')+'!E2!N!X / d.o.f. : '+strtrim(string(-2.0*chisq_best,format='(f18.2)'),1),/data


device,/close
set_plot,'x'

end