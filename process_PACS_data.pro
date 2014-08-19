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
; Read and bin PACS data

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
  readcol, 'HD181327_PACS_data/'+file[i],F=fmt,pacs_wave_in,v2,pacs_flux_in,v4,v5
  
  ; Prep bins
  n_pacs = n_elements(pacs_wave_in)
  min_pacs = min(pacs_wave_in)
  max_pacs = max(pacs_wave_in)
  range_pacs = max_pacs - min_pacs
  
  n_new = round(n_pacs/20.0)
  bin_width = range_pacs/n_new
  
  wave_new = dblarr(n_new)
  flux_new = dblarr(n_new)

  ;plot, pacs_wave_in,pacs_flux_in
  ;stop
  
  ; Integrate over each bin
  FOR j=0,(n_new-1) DO BEGIN
    start_val = min_pacs + j*bin_width
    stop_val = min_pacs + (j+1)*bin_width
    wave_new[j] = start_val + 0.5*bin_width
    
    wave_int = pacs_wave_in[where( (pacs_wave_in ge start_val) and (pacs_wave_in le stop_val) )]
    flux_int = pacs_flux_in[where( (pacs_wave_in ge start_val) and (pacs_wave_in le stop_val) )]
    
    flux_new[j] = int_tabulated(wave_int,flux_int)/(max(wave_int)-min(wave_int))
  ENDFOR
  
  ;plot, wave_new, flux_new
  ;stop
  
  ; Store the data
  pacs_wave[i] = wave_new;pacs_wave_in;wave_new
  pacs_flux[i] = flux_new;pacs_flux_in;flux_new

ENDFOR

plot, pacs_wave[0],pacs_flux[0],xrange=[1,220],yrange=[0,5.0],psym=1,/xlog,/ylog
oplot, pacs_wave[1],pacs_flux[1],psym=3
oplot, pacs_wave[2],pacs_flux[2],psym=2
oplot, pacs_wave[3],pacs_flux[3],psym=4

oplot, final_wave,final_spec
oplot, final_phot_wave,final_phot_fnu

stop
; *************************************************** ;
; Normalize MIPS SED data by comparing with MIPS70

; Define constants
c = 3.0e10 ; cm/s
to_cm = 1.0e-4

; Interpolate response curve at mips sed wavelengths
response_int = 10^(interpol(alog10(response),alog10(resp_wave*to_cm),alog10(mips_sed_wave*to_cm)))

; Define integrand
int1 = (c*1.0e-23*mips_sed_spec)/(mips_sed_wave*to_cm)^2 ; convert units
int2 = int1*response_int ; multiply flux and response

; Find integral of filter
int_filter = int_tabulated(mips_sed_wave*to_cm,response_int)

; Integrate the MIPS SED over the bandpass
synth_f70 = int_tabulated(mips_sed_wave*to_cm,int2)/int_filter

; Convert MIPS 70 to proper units
MIPS_as_flambda = (MIPS70_VAL*c*1.0e-23)/(71.42*to_cm)^2

; Find the normalization constant
norm = MIPS_as_flambda/synth_f70

; Normalize
mips_sed_spec = norm*mips_sed_spec

; *************************************************** ;
; Concatenate and save
mips_sed_source = MAKE_ARRAY(n_elements(mips_sed_wave), 1, /STRING, VALUE = 'SpitzerMIPS_SED')
FINAL_SOURCE = MAKE_ARRAY(n_elements(FINAL_WAVE), 1, /STRING, VALUE = 'SpitzerIRS')

; Find weighting for data points by looking at their spread
IRS_density = n_elements(FINAL_WAVE)/(max(FINAL_WAVE)-min(FINAL_WAVE))
mips_sed_density = n_elements(mips_sed_wave)/(max(mips_sed_wave)-min(mips_sed_wave))
weight = round(IRS_density/mips_sed_density)/2
print,string(name)+' Weight = '+strcompress(string(weight))
print,string(name)+' MIPS70 Value = '+strcompress(string(MIPS70_VAL))+' Jy'
print,string(name)+' MIPS70 Error = '+strcompress(string(MIPS70_ERROR))+' Jy'
print,""

; Append data structure recursively to account for data weight
FOR k = 0, weight-1 DO BEGIN
  FINAL_WAVE = [FINAL_WAVE, mips_sed_wave]
  FINAL_SPEC = [FINAL_SPEC, mips_sed_spec]
  FINAL_SPECERR = [FINAL_SPECERR, mips_sed_specerr]
  FINAL_SOURCE = [FINAL_SOURCE, mips_sed_source]
ENDFOR

; Sort arrays by wavelength
order = SORT(FINAL_WAVE)
FINAL_WAVE = FINAL_WAVE[order]
FINAL_SPEC = FINAL_SPEC[order]
FINAL_SPECERR = FINAL_SPECERR[order]
FINAL_SOURCE = FINAL_SOURCE[order]

; Save to new data file
SAVE,$
  FINAL_PHOT_FNU,FINAL_PHOT_WAVE, FINAL_WAVE, FINAL_SPEC, FINAL_SPECERR, FINAL_SOURCE,$
  filename='savfiles_MIPS_SED_corrected/'+name+'.sav'
  
; Mark that one object has been modified
count = count + 1

; Make a plot of the new spectrum
IF KEYWORD_SET(plot) THEN BEGIN
  plot_spectrum,name
ENDIF



print, 'Data updated for: '+string(count)+' objects.'
END