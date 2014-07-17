pro test_fit, single=single

seed = !Null

; Generate a random distribution

if keyword_set(single) then begin
  
  n = 6; number of fit parameters
  
  param_bnd = dblarr(2,6)
  param_bnd[*,0]= [10.0, 1000.0] ; Temperatures limited to reasonable values
  param_bnd[*,1]= [0.1, 30.0] ; .1 as lower bound
  param_bnd[*,2]= [16.5, 23.5] ; positive scale factors
  param_bnd[*,3]= [0, 1.0] ; olivine/pyroxene composition
  param_bnd[*,4]= [0, 1.0] ; crystalline fraction
  param_bnd[*,5]= [0, 1.0] ; forsterite/enstatite composition
  
endif

; Initialize Randomness
seed = !NULL
params = dblarr(n)

; Get random parameter values
for i=0,(n-1) do begin
  range = param_bnd[1,i] - param_bnd[0,i]
  params[i] = param_bnd[0,i] + range*randomu(seed)
endfor

; Load data for 
restore, 'savfiles_MIPS_SED_corrected/HD146897.sav'


if keyword_set(single) then begin
  ; Generate spectrum 
  FINAL_SPEC = (modelonegrain(FINAL_WAVE,params))*1e17
endif


save, FINAL_PHOT_FNU, FINAL_PHOT_WAVE, FINAL_SOURCE, $
      FINAL_SPEC, FINAL_SPECERR, FINAL_WAVE, $
      filename='savfiles_MIPS_SED_corrected/Test1.sav'

; Check the result
plot,FINAL_WAVE,FINAL_SPEC
FINAL_SPECERR[*,*] = 0
stop
; *************************************************** ;
; Perform setup operations

; Home directory
home_dir = '~/Summer2014/Silicate_Feature_Analysis/'
CD, home_dir

; Path
!PATH=!PATH+':'+Expand_Path('+'+home_dir)

; Functional directories
COMMON file_path, in_dir, out_dir, fit_name
in_dir = 'savfiles_MIPS_SED_corrected'
out_dir = '../Silicate_Feature_Analysis_output'

fit_name = 'single'
name = 'Test1'

; *************************************************** ;

fits_v1, name=name, fittype=fit_name
;plot_result, name



end