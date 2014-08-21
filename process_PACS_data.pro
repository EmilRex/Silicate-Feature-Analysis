; +
; NAME:
;  process_PACS_data
;
; PURPOSE:
;  Join Spitzer IRS data with Herschel PACS data. Scales Herschel PACS
;  data through comparison of integrated emission with MIPS 70 point.
;  Also weights the Herschel PACS data by looking at data/wavelength density.
;
; INPUTS:
;   NONE
;
; KEYWORDS
;   PLOT: Plot object spectrum
;
; OUTPUTS:
;   *NAME*.sav - IDL save file with updated data file sorted by wavelength
;
; AUTHORS:
;  Emil Christensen - chris2er@dukes.jmu.edu
;  Christine Chen - cchen@stsci.edu
;
; DISCLAIMER
;  This software is provided as is without any warranty whatsoever.
;  Permission to use, copy, modify, and distribute modified or
;  unmodified copies is granted, provided this disclaimer
;  is included unchanged.
;
; MODIFICATION HISTORY:
;  Written from process_MIPS_data by EC (8/18/14)
; -
; *************************************************** ;
pro process_PACS_data, plot=plot

name = 'HD181327'

; Read in response curve
fmt='f,f,f,f,f,f'
readcol,'MIPS_SED_response.txt',F=fmt,v1,v2,resp_wave,response,v5,v6, /SILENT

; Load object data from IDL savefile
restore, 'old_savfiles_mcmc/'+name+'.sav'

; *************************************************** ;
; Read PACS data

file = strarr(4)
fmt='f,f,f,f,f'

file[0] = 'hipe11_0_2454_calver49_OBSID_1342231318_B2A_2_TB_CS_PSF_PO.csv'
file[1] = 'hipe11_0_2454_calver49_OBSID_1342231723_B2B_2_TB_CS_PSF_PO.csv'
file[2] = 'hipe11_0_2454_calver49_OBSID_1342231318_R1_102_TB_CS_PSF_PO.csv'
file[3] = 'hipe11_0_2454_calver49_OBSID_1342231723_R1_102_TB_CS_PSF_PO.csv'

pacs_wave = list(length=4)
pacs_flux = list(length=4)

FOR i=0,3 DO BEGIN
  
  ; Read data
  readcol, 'HD181327_PACS_data/'+file[i],F=fmt,pacs_wave_in,v2,pacs_flux_in,v4,v5,/silent
  
  ; Store the data
  pacs_wave[i] = pacs_wave_in
  pacs_flux[i] = pacs_flux_in

ENDFOR

; *************************************************** ;
; Align PACS data

; Trim first 4 points ... flux seems way to low.
;pacs_wave[0]  = pacs_wave[0,4:n_elements(pacs_wave[0])-1]

; Overlap of region 0 and 1
; Scale to first sequesnce
r0_max = max(pacs_wave[0])
r1_min = min(pacs_wave[1])
r0_overlap = (pacs_flux[0])[where(pacs_wave[0] ge r1_min)]
r1_overlap = (pacs_flux[1])[where(pacs_wave[1] le r0_max)]
r0_med = median(r0_overlap)
r1_med = median(r1_overlap)
scale1 = r0_med/r1_med
pacs_flux[1] = scale1*pacs_flux[1]
pacs_wave[1] = (pacs_wave[1])[where(pacs_wave[1] gt r0_max)]

; Overlap of region 1 and 2
; Scale to first sequesnce
r1_max = max(pacs_wave[1])
r2_min = min(pacs_wave[2])
r1_overlap = (pacs_flux[1])[where(pacs_wave[1] ge r2_min)]
r2_overlap = (pacs_flux[2])[where(pacs_wave[2] le r1_max)]
r1_med = median(r1_overlap)
r2_med = median(r2_overlap)
scale2 = r1_med/r2_med
pacs_flux[2] = scale2*pacs_flux[2]
pacs_wave[1] = (pacs_wave[1])[where(pacs_wave[1] lt r2_min)]

; Overlap of region 2 and 3
; Scale to first sequesnce
r2_max = max(pacs_wave[2])
r3_min = min(pacs_wave[3])
r2_overlap = (pacs_flux[2])[where(pacs_wave[2] ge r3_min)]
r3_overlap = (pacs_flux[3])[where(pacs_wave[3] le r2_max)]
r2_med = median(r2_overlap)
r3_med = median(r3_overlap)
scale3 = r2_med/r3_med
pacs_flux[3] = scale3*pacs_flux[3]
pacs_wave[3] = (pacs_wave[3])[where(pacs_wave[3] gt r2_max)]

plot, pacs_wave[0],pacs_flux[0],xrange=[1,220],yrange=[0,5.0],psym=1
oplot, pacs_wave[1],pacs_flux[1],psym=3
oplot, pacs_wave[2],pacs_flux[2],psym=2
oplot, pacs_wave[3],pacs_flux[3],psym=4
oplot, final_wave,final_spec

pacs_wave_tot = [pacs_wave[0],pacs_wave[1],pacs_wave[2],pacs_wave[3]]
pacs_flux_tot = [pacs_flux[0],pacs_flux[1],pacs_flux[2],pacs_flux[3]]

; *************************************************** ;
; Bin PACS data
  
n_pacs = n_elements(pacs_wave_tot)
min_pacs = min(pacs_wave_tot)
max_pacs = max(pacs_wave_tot)
range_pacs = max_pacs - min_pacs

n_new = round(n_pacs/50.0)
bin_width = range_pacs/n_new

wave_new = dblarr(n_new)
flux_new = dblarr(n_new)
pacs_err = dblarr(n_new)

;plot, pacs_wave_tot,pacs_flux_tot,xrange=[30,220],yrange=[0,5.0],psym=2
;stop

; Integrate over each bin
FOR j=0,(n_new-1) DO BEGIN
  start_val = min_pacs + j*bin_width
  stop_val = min_pacs + (j+1)*bin_width
  wave_new[j] = start_val + 0.5*bin_width
  
  wave_int = pacs_wave_tot[where( (pacs_wave_tot ge start_val) and (pacs_wave_tot le stop_val) )]
  flux_int = pacs_flux_tot[where( (pacs_wave_tot ge start_val) and (pacs_wave_tot le stop_val) )]
  
  
  flux_new[j] = int_tabulated(wave_int,flux_int)/(max(wave_int)-min(wave_int))
  pacs_err[j] = stddev(flux_int)/sqrt(n_elements(flux_int))
ENDFOR

; *************************************************** ;
; Trim the first three points

wave_new = wave_new[3:(n_elements(wave_new)-1)]
flux_new = flux_new[3:(n_elements(flux_new)-1)]
pacs_err = pacs_err[3:(n_elements(flux_new)-1)]

; Trim the unusually high point
wave_new = wave_new[where( (flux_new lt 4.0) and (flux_new gt 0.05) )]
flux_new = flux_new[where( (flux_new lt 4.0) and (flux_new gt 0.05) )]
pacs_err = pacs_err[where( (flux_new lt 4.0) and (flux_new gt 0.05) )]

; Remove outlier near 145 microns
wave_new = wave_new[where( flux_new ne max(flux_new[where(wave_new gt 140.0 and wave_new  lt 150.0)]) )]
flux_new = flux_new[where( flux_new ne max(flux_new[where(wave_new gt 140.0 and wave_new  lt 150.0)]) )]
pacs_err = pacs_err[where( flux_new ne max(flux_new[where(wave_new gt 140.0 and wave_new  lt 150.0)]) )]

; Store the data
pacs_wave = wave_new
pacs_flux = flux_new

; Test output
;plot, wave_new, flux_new,xrange=[1,220],yrange=[0,5.0],psym=1
;oplot,final_wave,final_spec
;oplot,[71.42],[MIPS70_VAL],psym=6
;stop

; *************************************************** ;
; Normalize MIPS SED data by comparing with MIPS70

; Define constants
c = 3.0e10 ; cm/s
to_cm = 1.0e-4

; Interpolate response curve at mips sed wavelengths
response_int = 10^(interpol(alog10(response),alog10(resp_wave*to_cm),alog10(pacs_wave*to_cm)))

; Define integrand
int1 = (c*1.0e-23*pacs_flux)/(pacs_wave*to_cm)^2 ; convert units
int2 = int1*response_int ; multiply flux and response

; Find integral of filter
int_filter = int_tabulated(pacs_wave*to_cm,response_int)

; Integrate the MIPS SED over the bandpass
synth_f70 = int_tabulated(pacs_wave*to_cm,int2)/int_filter

; Convert MIPS 70 to proper units
MIPS_as_flambda = (MIPS70_VAL*c*1.0e-23)/(71.42*to_cm)^2

; Find the normalization constant
norm = MIPS_as_flambda/synth_f70

; Normalize
pacs_flux = norm*pacs_flux

; Test output
;plot, pacs_wave, pacs_flux,xrange=[1,220],yrange=[0,5.0],psym=1
;oplot,final_wave,final_spec
;oplot,[71.42],[MIPS70_VAL],psym=6
;stop

; *************************************************** ;
; Subtract off the photosphere

phot_int = 10^interpol(alog10(final_phot_fnu),alog10(final_phot_wave*to_cm),alog10(pacs_wave*to_cm))

pacs_flux = pacs_flux - phot_int

; Test output
;plot, pacs_wave, pacs_flux,xrange=[1,220],yrange=[0,5.0],psym=1,/xlog,/ylog
;oplot,final_wave,final_spec
;oplot,final_phot_wave,final_phot_fnu
;stop

; *************************************************** ;
; Concatenate and save
pacs_source = MAKE_ARRAY(n_elements(pacs_wave), 1, /STRING, VALUE = 'Herschel_PACS')
FINAL_SOURCE = MAKE_ARRAY(n_elements(FINAL_WAVE), 1, /STRING, VALUE = 'SpitzerIRS')

; Find weighting for data points by looking at their spread
;IRS_density = n_elements(FINAL_WAVE)/(max(FINAL_WAVE)-min(FINAL_WAVE))
;mips_sed_density = n_elements(pacs_wave)/(max(pacs_wave)-min(pacs_wave))
;weight = round(IRS_density/mips_sed_density)/2
;print,string(name)+' Weight = '+strcompress(string(weight))
;print,string(name)+' MIPS70 Value = '+strcompress(string(MIPS70_VAL))+' Jy'
;print,string(name)+' MIPS70 Error = '+strcompress(string(MIPS70_ERROR))+' Jy'
;print,""
weight = 1
;pacs_error = MAKE_ARRAY(n_elements(pacs_wave), 1, /FLOAT, VALUE = stddev(pacs_flux))

; Append data structure recursively to account for data weight
FOR k = 0, weight-1 DO BEGIN
  FINAL_WAVE = [FINAL_WAVE, pacs_wave]
  FINAL_SPEC = [FINAL_SPEC, pacs_flux]
  FINAL_SPECERR = [FINAL_SPECERR, pacs_err]
  FINAL_SOURCE = [FINAL_SOURCE, pacs_source]
ENDFOR

plot,FINAL_WAVE,FINAL_SPEC,xrange=[1,220],yrange=[0,5.0],psym=1
oploterr,FINAL_WAVE,FINAL_SPEC,FINAL_SPECERR
stop

; Sort arrays by wavelength
order = SORT(FINAL_WAVE)
FINAL_WAVE = FINAL_WAVE[order]
FINAL_SPEC = FINAL_SPEC[order]
FINAL_SPECERR = FINAL_SPECERR[order]
FINAL_SOURCE = FINAL_SOURCE[order]

;stop

; Save to new data file
SAVE,$
  FINAL_PHOT_FNU,FINAL_PHOT_WAVE, FINAL_WAVE, FINAL_SPEC, FINAL_SPECERR, FINAL_SOURCE,$
  filename='savfiles_PACS/'+name+'.sav'

;stop

; Make a plot of the new spectrum
;IF KEYWORD_SET(plot) THEN BEGIN
;  plot_spectrum,name
;ENDIF

END