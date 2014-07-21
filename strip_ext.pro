; Separated from process_data.pro
; By EC on 7/21/14
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