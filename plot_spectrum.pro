; +
; NAME:
;  plot_spectrum.pro
;
; PURPOSE:
;  Plot a spectrum of the given object and save as post script to /plots
;
; INPUTS:
;   OBJ_NAME: Name of object
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
;  Modified by EC (7/2/2014)
; -
; *************************************************** ;

pro plot_spectrum, object_name

  restore,'savfiles_MIPS_SED_corrected/'+object_name+'.sav'
  
  wave_irs = final_wave ;[final_wave,71.42]
  fl_diff = final_spec ;[final_spec,mips70_val]
  uncer_irs = final_specerr ;[final_specerr,mips70_error]
  
  restore, 'old_savfiles_MCMC/'+object_name+'.sav'
  
  ; *************************************************** ;
  ; Separate out data types and set colors/markers
  
  IRS_wave = wave_irs[where(strmatch(final_source, 'SpitzerIRS') EQ 1)]
  IRS_spec = fl_diff[where(strmatch(final_source, 'SpitzerIRS') EQ 1)]
  IRS_color = 1 ; red
  IRS_psym = 1 ; plus
  
  MIPS70_wave = [71.42,71.42]
  MIPS70_spec = [MIPS70_VAL,MIPS70_VAL]
  MIPS70_color = 2 ; red
  MIPS70_psym = 6 ; square
  
  MIPS_SED_wave = wave_irs[where(strmatch(final_source, 'SpitzerMIPS_SED') EQ 1)]
  MIPS_SED_spec = fl_diff[where(strmatch(final_source, 'SpitzerMIPS_SED') EQ 1)]
  MIPS_SED_color = 3 ; blue
  MIPS_SED_psym = 5 ; triangle
  
  ; *************************************************** ;
  ; Begin plotting
  
  ; Set up device
  set_plot,'PS'
  device, filename ='plots/'+object_name+'_spectrum_corrected.ps',/COLOR,/HELVETICA,XSIZE=15,YSIZE=12.5 & !p.font =0
  loadct,39
  !p.background=16777215
  
  ;Make IRS Spectrum Plot w/ Errorbars
  TVLCT,[0,255,0,0],[0,0,255,0],[0,0,0,255]
  
  ; Plot data points
  plot,wave_irs,fl_diff,title=object_name+' Spectrum', $
    ystyle=1,psym=4,xstyle=1,xtitle='Wavelength ('+cggreek('mu')+'m)', $
    ytitle='F'+cggreek('nu')+' (Jy)',charthick=1, thick=1, $
    xthick=2, ythick=2, charsize=1,color=0
    
  ; Connect the points
  oplot,wave_irs,fl_diff,color=0
  
  ; Overlay data with different colors and markers
  oplot,IRS_wave,IRS_spec,color=IRS_color,psym=IRS_psym
  oplot,MIPS70_wave,MIPS70_spec,color=MIPS70_color,psym=MIPS70_psym
  oplot,MIPS_SED_wave,MIPS_SED_spec,color=MIPS_SED_color,psym=MIPS_SED_psym
  
  ; Add error bars
  oploterr,wave_irs,fl_diff,uncer_irs;,psym=1;,color=0
  
  ; Create legend
  ;LEGEND,TEXTOIDL('\chi^{2}/d.o.f. : ')+  strtrim(string(tmp1,format='(f18.2)'),1), /top, /left,color=1,linestyle=lines[1]
  
  device,/close
  set_plot,'x'
  
END