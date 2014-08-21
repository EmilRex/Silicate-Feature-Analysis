; ****************************************************************************************************** ;
pro main
; ****************************************************************************************************** ;
; Perform setup operations

; Set to directory containing repository
home_dir = '~/'
CD, home_dir

; Expand IDL path
!PATH=!PATH+':'+Expand_Path('+'+home_dir)

; Functional directories
COMMON file_path, in_dir, out_dir, fit_name, object_name
; Set to directory containing data
; Default set to savfiles
in_dir = 'savfiles'

; Set directory location for output. 
; Location should be relative to your home directory
; Default set to output
out_dir = 'output'


; *************************************************** ;
;RUN
; *************************************************** ;
; Set name of object(s) to be run
names = ['']

; Set fits to run.
; Default is all
fit_names = ['single','multi','disk']

FOREACH name, names DO BEGIN
  FOREACH fit_name, fit_names DO BEGIN

    object_name = name
    print, object_name

    fits, name=name, fittype=fit_name
    
    plot_result

    null = display_results(name)

  ENDFOREACH
ENDFOREACH


end
