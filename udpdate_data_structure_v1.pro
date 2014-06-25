; +
; NAME:
;  udpdate_data_structure_v1
;
; PURPOSE:
;  Add new data to old data savefiles
;
; INPUTS:
;   SOURCE: Data source identifier
;           - SpitzerIRS_[SL,LL,SH,LH]
;           - SpitzerMIPS_[PHOT,SED]
;           - Herschel_PACS_[PHOT,SPEC]
;           - Herschel_SPIRE
;   DIR_IN: Directory containing input files
;   DIR_OUT: Directory for output files
;
; KEYWORDS
;   INITIAL: Take in old data and add MIPS_70 value to data array
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
;  Written by EC (6/24/2014)
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

pro udpdate_data_structure_v1,$
  dir_in=dir_in, dir_out=dir_out, new_dir=new_dir,$
  initial=initial, new_data=new_data,$
  source_name=source_name

; Set iterators
loop_one = 1
initial_count = 0
new_count = 0

; *************************************************** ;
; Check for directory
test = FILE_TEST(dir_out, /DIRECTORY)
; If directory does not exist, print warning
IF (test eq 0) THEN BEGIN; Directory does not exist
  print, 'Error! The directory does not exist.'
  print, 'Please create the desired directory: /'+dir_out
  return
ENDIF

; Get names of input objects
name_list = strip_ext(dir_in,'sav')

; *************************************************** ;
; Loop for each name
FOR i=0, (n_elements(name_list)-1) DO BEGIN

  ; Load object data from IDL savefile
  restore, dir_in+'/'+name_list[i]+'.sav'

; *************************************************** ;
; Call Example:
 ;udpdate_data_structure_v1, dir_in='old_savfiles_mcmc',dir_out='savfiles',/initial
  ; If first time, destroy old data structure
  IF KEYWORD_SET(initial) THEN BEGIN

    ; Test the assumption that len(wave)=len(spec)=len(err)
    IF (n_elements(FINAL_WAVE) ne n_elements(FINAL_SPEC)) OR (n_elements(FINAL_SPEC) ne n_elements(FINAL_SPECERR)) OR (n_elements(FINAL_SPECERR) ne n_elements(FINAL_WAVE)) THEN BEGIN
      print, 'Error: Check diminsionality of '+name_list[i]+' data arrays!'
      return
    ENDIF
      
    ; Put source into FINAL_SOURCE
    FINAL_SOURCE = MAKE_ARRAY(n_elements(FINAL_WAVE), 1, /STRING, VALUE = 'SpitzerIRS')

    ; Check for MIPS_70 value
    IF (n_elements(MIPS70_VAL) eq 1) AND (n_elements(MIPS70_ERROR) eq 1) THEN BEGIN
    
      ; Append FINAL_[WAVE,SPEC,SPECERR,SOURCE]
      FINAL_WAVE = [FINAL_WAVE, 71.4200] ; 71.42 = average lambda according to instrument handbook
      FINAL_SPEC = [FINAL_SPEC, MIPS70_VAL]
      FINAL_SPECERR = [FINAL_SPECERR, MIPS70_ERROR]
      FINAL_SOURCE = [FINAL_SOURCE, 'Spitzer_MIPS_Phot']
    ENDIF ELSE BEGIN
      print, 'No MIPS photometry data found for '+name_list[i]
    ENDELSE
    
    ; Save to new data file
    SAVE,$
    FINAL_PHOT_FNU,FINAL_PHOT_WAVE, FINAL_WAVE, FINAL_SPEC, FINAL_SPECERR, FINAL_SOURCE,$
    filename=dir_out+'/'+name_list[i]+'.sav'
      
    initial_count = initial_count+1
  ENDIF 
  
; *************************************************** ;
  ; If not first time, merely update data 
  ; Call Example:
  ;udpdate_data_structure_v1, dir_in='savfiles',new_dir='new_MIPS_data',dir_out='savfiles_MIPS_SED', /new_data,source_name = 'SpitzerMIPS_SED'
  
  IF KEYWORD_SET(new_data) THEN BEGIN
    
    ; Strip names from directory on first loop
    IF (loop_one eq 1) THEN BEGIN
      new_name = strip_ext(new_dir,'dat')
      loop_one = 0 ; Identify first loop as done
    ENDIF
        
    ; Find match to current object - see temp template
    FOR j = 0, n_elements(new_name)-1 DO BEGIN
      IF (name_list[i] eq new_name[j]) THEN BEGIN
        
        ; Read in new data
        fmt='f,f,f,f,f'
        readcol, new_dir+'/'+new_name[j]+'.dat',F=fmt,NEW_WAVE,NEW_SPEC,NEW_SPECERR,NEW_SN, /SILENT
        
        ; Check to see arrays are consistent, then create source column
        IF (n_elements(NEW_WAVE) ne n_elements(NEW_SPEC)) OR (n_elements(NEW_SPEC) ne n_elements(NEW_SPECERR)) OR (n_elements(NEW_SPECERR) ne n_elements(NEW_WAVE)) THEN BEGIN
          print, 'Error: Check diminsionality of '+new_name[j]+' data arrays!'
          return
        ENDIF ELSE BEGIN
          NEW_SOURCE = MAKE_ARRAY(n_elements(NEW_WAVE), 1, /STRING, VALUE = source_name)
        ENDELSE
        
        ; Append data structure
        FINAL_WAVE = [FINAL_WAVE, NEW_WAVE]
        FINAL_SPEC = [FINAL_SPEC, NEW_SPEC]
        FINAL_SPECERR = [FINAL_SPECERR, NEW_SPECERR]
        FINAL_SOURCE = [FINAL_SOURCE, NEW_SOURCE]
        
        ; Sort arrays by wavelength
        order = SORT(FINAL_WAVE)
        FINAL_WAVE = FINAL_WAVE[order]
        FINAL_SPEC = FINAL_SPEC[order]
        FINAL_SPECERR = FINAL_SPECERR[order]
        FINAL_SOURCE = FINAL_SOURCE[order]
        
        ; Save to new data file
        SAVE,$
        FINAL_PHOT_FNU,FINAL_PHOT_WAVE, FINAL_WAVE, FINAL_SPEC, FINAL_SPECERR, FINAL_SOURCE,$
        filename=dir_out+'/'+name_list[i]+'.sav'
        
        ; Mark that one object has been modified
        new_count = new_count + 1

        plot_spectrum,new_name[j],FINAL_WAVE,FINAL_SPEC,FINAL_SPECERR
        BREAK
        
      ENDIF
    ENDFOR    

  ENDIF
  
ENDFOR

; *************************************************** ;
; Print objects affected
IF KEYWORD_SET(initial) THEN BEGIN
  print, 'Initial - Data updated for: '+string(initial_count)+' objects.'
ENDIF
IF KEYWORD_SET(new_data) THEN BEGIN
  print, 'New - Data updated for: '+string(new_count)+' objetcs.'
ENDIF

; End program
END