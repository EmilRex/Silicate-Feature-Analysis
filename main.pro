; ****************************************************************************************************** ;
pro main
; ****************************************************************************************************** ;
; Perform setup operations

; Home directory
home_dir = '~/Summer2014/Silicate_Feature_Analysis/'
CD, home_dir

; Path
!PATH=!PATH+':'+Expand_Path('+'+home_dir)

; Functional directories
COMMON file_path, in_dir, out_dir, fit_name, object_name
in_dir = 'savfiles_MIPS_SED_corrected'
;in_dir = 'savfiles_PACS'

out_dir = '../Silicate_Feature_Analysis_output'
;out_dir = '../Science3_output'


; *************************************************** ;
;RUN MULTIPLE
; *************************************************** ;

;names = ['HD95086','HD106906','HD108257','HD110058','HD111520','HD113556','HD113766','HD114082','HD115600','HD117214','HD145560','HD146181','HD146897']
names = ['HD95086','HD110058','HD113556','HD114082','HD115600','HD117214','HD145560','HD146181','HD146897']
;names = ['HD117214'];,'HD146897']
;names = ['HD181327']
fit_names = ['single','multi','disk']
;fit_names = ['disk'];['multi','single'];,'disk']

FOREACH fit_name, fit_names DO BEGIN
  FOREACH name, names DO BEGIN

    object_name = name
    fits_v1, name=name, fittype=fit_name
    ;print, object_name
    
    ; !!! Before plotting, make sure pro/legend.pro is compiled
    ;plot_result_PACS,/separate
    ;plot_result_MIPS,/separate
    
    ;null = display_historic(name)
    ;null = display_results(name)
    ;stop
  ENDFOREACH
ENDFOREACH


end
