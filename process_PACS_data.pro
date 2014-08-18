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
; Read, concatenate and bin PACS data
fmt='f,f,f,f,f'

file = strarr(4)


file[0] = 'hipe11_0_2454_calver49_OBSID_1342231318_B2A_2_TB_CS_PSF_PO.txt'
file[1] = 'hipe11_0_2454_calver49_OBSID_1342231723_B2B_2_TB_CS_PSF_PO.txt'
file[2] = 'hipe11_0_2454_calver49_OBSID_1342231318_R1_102_TB_CS_PSF_PO.txt'
file[3] = 'hipe11_0_2454_calver49_OBSID_1342231723_R1_102_TB_CS_PSF_PO.txt'

readcol, 'HD181327_PACS_data/'+name+'.dat',F=fmt,mips_sed_wave,mips_sed_spec,Mips_sed_specerr,mips_sed_sn;, /SILENT



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