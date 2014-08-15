; +
; NAME:
;  printmodels_v1
;
; PURPOSE:
;  Restore save data for specified object and pass to modelfit_v1
;  via 'params' structure.
;
; INPUTS:
;  NAME1:  Name of object to be analyzed. Main identifier.
;  FITTYPE: Model to be fitted (in order of appearance below)
;          1) 'single' - Single grain model
;          2) 'multi' - Two discrete grain populations, (T1,a1), (T2,a2)
;          3) 'disk' - Continuous radial distribution of particles 
;  PS:
;  NSTEP:
;  THIN_VAL:
;  TEFF: 
;  DIST:
;  AMIN: Minimum grain size (units?)
;  SEED: 
;  CNT:
;  SCALE_VAL:
;  NUM_CHAINS: Number of MCMC chains to run
;
; OUTPUTS:
;   NONE
;
; AUTHORS:
;  Tushar Mittal - mittal.tushar22@gmail.com
;  Christine Chen - cchen@stsci.edu
;  Emil Christensen - chris2er@dukes.jmu.edu
;
; DISCLAIMER
;  This software is provided as is without any warranty whatsoever.
;  Permission to use, copy, modify, and distribute modified or
;  unmodified copies is granted, provided this disclaimer
;  is included unchanged.
;
; MODIFICATION HISTORY:
;  Written by TM (June 2013) as printmodels_new.pro
;  Organized and commented by EC (6/24/2014)
; -
; *************************************************** ;
pro printmodels_v1, name1, fittype=fittype,ps=ps,nstep=nstep,$
  thin_val=thin_val,teff=teff,dist=dist,amin=amin,seed=seed,$
  cnt=cnt,scale_val =scale_val,num_chains=num_chains

COMMON file_path, in_dir, out_dir

; Assign default if fit type not set
if not keyword_set(fittype) then begin
  print, 'Error: fit type not set. Exiting from printmodels'
  return
endif

; Load object data from IDL savefile
restore,in_dir+'/'+name1+'.sav' ; Changed to new directory name for the time being

; Fill 'params' structure with data passed by fits
; and with data from object savefile
params={distval:dist,$
        teff_val:teff,$
        amin_val:amin,$
        lambdafit:final_wave,$
        diskfit:final_spec,$
        errfit:final_specerr,$
        name_obj:name1,$
        seed : seed,$
        nstep :nstep,$
        thin_val:thin_val,$
        ps :  ps,$
        cnt:cnt,$
        scale_val:scale_val,$
        num_chains: num_chains}
        
; Print fit type to console
print,fittype

; Pass 'params' to modelfit
modelfit_v1,param_struc=params,fittype=fittype

; End program
return
end








