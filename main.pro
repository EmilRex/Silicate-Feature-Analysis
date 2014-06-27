pro main

; Set up directory
home_dir = '~/Summer2014/Silicate_Feature_Analysis/'
CD, home_dir

; Set up path
!PATH=!PATH+':'+Expand_Path('+'+home_dir)

; Run program
name = 'HD105857'
fittype = 'multi_mips'
fits_v1, name=name, fittype=fittype

; Get Data
; Gives error for cggreek - variable undefined - when run from command line, but not GUI
;udpdate_data_structure_v1, dir_in='savfiles',new_dir='new_MIPS_data',dir_out='savfiles_MIPS_SED', /new_data,source_name = 'SpitzerMIPS_SED'


end