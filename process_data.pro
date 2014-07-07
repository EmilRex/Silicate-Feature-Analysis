; +
; NAME:
;  process_data
;
; PURPOSE:
;  Join Spitzer IRS data with Spitzer MIPS SED data. Scales MIPS SED
;  data through comparison of integrated emission with MIPS 70 point.
;  Also weights the MIPS SED data by looking at data/wavelength density.
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
;  Written by EC (7/2/2014)
; -
; *************************************************** ;
  
FUNCTION strip_ext, dir, ext
  ; ext - sav,dat,etc
  ; Strip all object names in directory
  
  CD, dir
  spawn, 'ls *.'+ext, names ; Get list of file names
  CD, '..'
  
  len = intarr(n_elements(names))
  ext_len = strlen(ext)
  
  FOR i=0, (n_elements(names)-1) DO BEGIN
    len[i] = strlen(names[i])
    ; Strip extension from file names
    names[i] = strmid(names[i],0,(len[i]-ext_len-1))
  ENDFOR
  
  return, names
END


; *************************************************** ;
pro process_data, plot=plot

; Set iterators
count = 0

; Get names of input objects
name_list = strip_ext('new_MIPS_data','dat')

; Get list of potential matches
matches = strip_ext('old_savfiles_mcmc','sav')

; Read in response curve
fmt='f,f,f,f,f,f'
readcol,'MIPS_SED_response.txt',F=fmt,v1,v2,resp_wave,response,v5,v6, /SILENT


; *************************************************** ;
; Loop for each name
FOR i=0, (n_elements(name_list)-1) DO BEGIN
  FOR j=0, (n_elements(matches)-1) DO BEGIN
    
    IF (name_list[i] eq matches[j]) THEN BEGIN
  
      ; Read in MIPS SED data
      fmt='f,f,f,f,f'
      readcol, 'new_MIPS_data/'+name_list[i]+'.dat',F=fmt,mips_sed_wave,mips_sed_spec,Mips_sed_specerr,mips_sed_sn, /SILENT
      
      ; Load object data from IDL savefile
      restore, 'old_savfiles_mcmc/'+name_list[i]+'.sav'
    
    
      ; *************************************************** ;
      ; Test input data
    
      ; Test the assumption that len(wave)=len(spec)=len(err)
      IF (n_elements(FINAL_WAVE) ne n_elements(FINAL_SPEC)) OR $ 
         (n_elements(FINAL_SPEC) ne n_elements(FINAL_SPECERR)) OR $ 
         (n_elements(FINAL_SPECERR) ne n_elements(FINAL_WAVE)) THEN BEGIN
        print, 'Error: Check diminsionality of '+name_list[i]+' data arrays!'
        return
      ENDIF
    
      ; Check to see arrays are consistent, then create source column
      IF (n_elements(mips_sed_wave) ne n_elements(mips_sed_spec)) OR $ 
         (n_elements(mips_sed_spec) ne n_elements(mips_sed_specerr)) OR $ 
         (n_elements(mips_sed_specerr) ne n_elements(mips_sed_wave)) THEN BEGIN
        print, 'Error: Check diminsionality of '+name_list[i]+' data arrays!'
        return
      ENDIF
      
      IF (n_elements(MIPS70_VAL) ne 1) AND (n_elements(MIPS70_ERROR) ne 1) OR $
         (MIPS70_VAL le 0.0) OR (MIPS70_ERROR le 0.0)  THEN BEGIN
        print, 'No MIPS photometry data found for '+name_list[i]
        return
      ENDIF
      
      ; *************************************************** ;
      ; Normalize MIPS SED data by comparing with MIPS70
      
;      c = 3.0e10 ; cm/s
;      to_cm = 1.0e-4
      
      ; Find delta lambda
;      dlambda = (max(mips_sed_wave)-min(mips_sed_wave))*to_cm
;      nu_to_lambda = c/(mips_sed_wave*to_cm)^2

      ; Interpolate responses at MIPS SED wavelengths
;      interpolated_response = interpol(response,resp_wave,mips_sed_wave)

      ; Adjust
;      synthetic_F70 = interpolated_response*mips_sed_spec*nu_to_lambda
;      tot_expected_flux = INT_TABULATED(mips_sed_wave,synthetic_F70)/dlambda

;      MIPS70_lambda = float((MIPS70_VAL*(c/(71.42*to_cm)^2)))

      ; Find scaling factor
;      scale_factor = MIPS70_lambda/tot_expected_flux

      intflux=10^(interpol(alog10(mips_sed_spec),alog10(mips_sed_wave*1.0e-4),alog10(resp_wave)))
      flux70=int_tabulated(resp_wave,intflux*response)/int_tabulated(resp_wave,response)*1.0e26

      scale_factor = MIPS70_VAL/flux70
      

stop
      ; Scale
      mips_sed_spec = scale_factor*mips_sed_spec
        
      ; *************************************************** ;
      ; Concatenate and save
      mips_sed_source = MAKE_ARRAY(n_elements(mips_sed_wave), 1, /STRING, VALUE = 'SpitzerMIPS_SED')
      FINAL_SOURCE = MAKE_ARRAY(n_elements(FINAL_WAVE), 1, /STRING, VALUE = 'SpitzerIRS')
      
      ; Find weighting for data points by looking at their spread
      IRS_density = n_elements(FINAL_WAVE)/(max(FINAL_WAVE)-min(FINAL_WAVE))
      mips_sed_density = n_elements(mips_sed_wave)/(max(mips_sed_wave)-min(mips_sed_wave))
      weight = FLOOR(IRS_density/mips_sed_density)
      
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
        filename='savfiles_MIPS_SED_corrected/'+name_list[i]+'.sav'
        
      ; Mark that one object has been modified
      count = count + 1
      
      ; Make a plot of the new spectrum
      IF KEYWORD_SET(plot) THEN BEGIN
        plot_spectrum,name_list[i]
      ENDIF
      
    ENDIF 
  ENDFOR
ENDFOR


print, 'Data updated for: '+string(count)+' objects.'
END