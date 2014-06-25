; +
; NAME:
;  plot_spectrum.pro
;
; PURPOSE:
;  Plot a spectrum of the given object and save as post script to /plots
;
; INPUTS:
;   NAME: Name of object
;   WAVE: Array of wavelength values (x-axis)
;   SPEC: Array of fluxes for the given wavelengths (y-axis)
;   SPECERR: Array of errors associated with the array of fluxes
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
;  Modified by EC (6/25/2014)
; -
; *************************************************** ;
pro plot_spectrum,name,wave,spec,specerr,source

; Set device parameters
set_plot,'PS'
device, filename ='plots/'+name+'_spectrum.ps',/COLOR,/HELVETICA,XSIZE=15,YSIZE=12.5 & !p.font =0

;Define colors
TVLCT,[0,255,0,0],[0,0,255,0],[0,0,0,255]

; *************************************************** ;
; Separate out data types and set colors/markers

IRS_wave = wave[where(strmatch(source, 'SpitzerIRS') EQ 1)]
IRS_spec = spec[where(strmatch(source, 'SpitzerIRS') EQ 1)]
IRS_color = 1 ; red
IRS_psym = 1 ; plus

MIPS70_wave = wave[where(strmatch(source, 'Spitzer_MIPS_Phot') EQ 1)]
MIPS70_spec = spec[where(strmatch(source, 'Spitzer_MIPS_Phot') EQ 1)]
MIPS70_color = 2 ; red
MIPS70_psym = 6 ; square

MIPS_SED_wave = wave[where(strmatch(source, 'SpitzerMIPS_SED') EQ 1)]
MIPS_SED_spec = spec[where(strmatch(source, 'SpitzerMIPS_SED') EQ 1)]
MIPS_SED_color = 3 ; blue
MIPS_SED_psym = 5 ; triangle


; *************************************************** ;
; Plot    
plot,wave,spec,title=name,ystyle=1,xstyle=1,xtitle=('Wavelength (\mum)'),ytitle=('F_\nu (Jy)'),charthick=1, thick=1, xthick=2, ythick=2,psym=0, charsize=1,color=0

; Overlay data with different colors and markers
oplot,IRS_wave,IRS_spec,color=IRS_color,psym=IRS_psym
oplot,MIPS70_wave,MIPS70_spec,color=MIPS70_color,psym=MIPS70_psym
oplot,MIPS_SED_wave,MIPS_SED_spec,color=MIPS_SED_color,psym=MIPS_SED_psym
oploterr,wave,spec,specerr,0


;LEGEND, /best, color=1, linestyle=lines[1]

; *************************************************** ;
; End
device,/close
set_plot,'x'

end