; +
; NAME:
;  fits_v1
;
; PURPOSE:
;  Script called from main.pro. Runs appropriate pso and mcmc versions.
;
; INPUTS:
;  NAME:  Name of object to be analyzed. Main identifier.
;  FITTYPE: Model to be fitted (in order of appearance below)
;          1) 'single' - Single grain model
;          2) 'multi' - Two discrete grain populations, (T1,a1), (T2,a2)
;          3) 'disk' - Continuous radial distribution of particles 
;
; KEYWORDS:
;   NONE
;
; OUTPUTS:
;   NONE
;
; AUTHORS:
;  Tushar Mittal
;  Christine Chen
;  Emil Christensen - chris2er@dukes.jmu.edu
;
; DISCLAIMER
;  This software is provided as is without any warranty whatsoever.
;  Permission to use, copy, modify, and distribute modified or
;  unmodified copies is granted, provided this disclaimer
;  is included unchanged.
;
; MODIFICATION HISTORY:
;  Written by TM (June 2013) as fits_new.pro
;  Organized and commented by EC (6/23/2014)
; -
; *************************************************** ;
pro fits_v1, name=name,fittype=fittype

COMMON file_path, in_dir, out_dir, fit_name, object_name
COMMON star_params, t_star, dist_to_star, star_lambda, star_fnu

; Load text file with names, Teff, amin and dist
fmt='a,f,f,f'
readcol,'input_files/input_param_file.txt',F=fmt,catalog_name,c_teff,c_amin,c_dist_val,/silent

; Match loaded data names with name of object being analyzed 
amin = .1
for i = 0, size(catalog_name,/n_elements)-1 do begin
   if (catalog_name[i] eq object_name) then begin
       Teff=c_teff[i]
       amin=c_amin[i]
       dist_val=c_dist_val[i]
    endif
endfor

; Load stellar photosphere model data
restore,in_dir+'/'+object_name+'.sav'

t_star = Teff
dist_to_star = dist_val
star_lambda = final_phot_wave
star_fnu = final_phot_fnu

; *************************************************** ;
; Set simulation best values

best_val1 = [99999,99999] ; for 'single'
best_val2 = [99999,99999] ; for 'multi'
best_val3 = [99999,99999] ; for 'disk'

; *************************************************** ;
; Begin fittype cases

case fittype of
; ****************************************************************************************************** ;
     'single': begin
; ****************************************************************************************************** ;      

          ; *************************************************** ;
          ; Short MCMC run for prior data
          val=dblarr(8,7) ; Define aray to hold fit values
          go = 1.0 ; Proceed to pso 
          i = 0.0 ; Begin iteration
          yes = 0.0 ; Do not procedd to MCMC simulation yet

          ; If prior best fit loaded, run a short MCMC to check for convergence
          if (best_val1[0] lt 99999 ) then begin

             ; Call MCMC routine with 100 steps of 10 chains
             printmodels_v1,name,fittype=[1],ps='X',nstep=100,thin_val=1,teff=Teff,$
             dist=dist_val,amin=amin,seed=best_val1[0:6],cnt=1,scale_val=0,num_chains=10

             ; Read in and store data from printmodels_v1 called just above
             data = readfits(out_dir+'/'+name+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=1)

             ; ???
             if ( min(-2.*data(7,99,*)) le 5.0) then begin
                yes = 1.0 ; Proceed to MCMC routine
                go = -1.0 ; Skip PSO routine
             endif

             ; ???
             if ( min(-2.*data(7,99,*)) gt 100.0) then begin
                yes=0.0 ; Forgo MCMC routine and exit to main
                go = -1.0 ; Skip PSO routine
             endif
          endif
          
          ; *************************************************** ;
          ; PSO Routine
          while (go gt 0.0) and ( i lt 7.0) do begin
             pso_min_example,name=name,teff=Teff,amin=amin,output= value,sequence = i ; Call PSO
             val[*,i]=value ; Save PSO's optimized values for the run

             ; If ??? proceed to MCMC             
             if(val[7,i] le 1.00) and ( i le 2.0) then begin
                 go = -1.0 ; Exit PSO routine
                 yes =1.0 ; Proceed to MCMC routine
                 best_val1 = value ; Update best value
             endif

             ; If ??? proceed to MCMC
             if(val[7,i] le 1.50) and ( i gt 2.0) then begin
                 go = -1.0 ; Exit PSO routine
                 yes =1.0 ; Proceed to MCMC routine
                 best_val1 = value ; Update best value
              endif

              ; If ??? exit to main
             if(val[7,0] ge 300.0) and  (val[7,1] ge 300.0) and (val[7,2] ge 300.0)  and (i eq 2) then begin
                 go = -1.0 ; Exit PSO routine
                 yes =0.0 ; Forgo MCMC routine and exit to main
             endif

             ; If ??? proceed to exit to main
             if (val[7,0] ge 100.0) and  (val[7,1] ge 100.0) and (val[7,2] ge 100.0) and (val[7,3] ge 100.0) and  (val[7,4] ge 100.0) and (val[7,5] ge 100.0) and (i eq 5) then begin
                     go = -1.0 ; Exit PSO routine
                     yes = 0.0 ; Forgo MCMC routine and exit to main
              endif
                 
              ; If at last run of PSO, proceed to MCMC
              if ( i eq 6) then begin
                     abc = min(val[7,*],indx)
                     best_val1=val[*,indx] ; Set best value to best of all 6 runs
                     yes = 1.0 ; Proceed to MCMC routine
                     go  = -1.0 ; Exit PSO routine
               endif
               i++ ; iterate
            endwhile

          ; *************************************************** ;
          ; MCMC Routine
          if (yes eq 1.0) then begin
              ; Call MCMC routine with 6000 steps of 10 chains
              printmodels_v1,name,fittype=[1],ps='X',nstep=5000,thin_val=1,teff=Teff,$
              dist=dist_val,amin=amin,seed=best_val1[0:6],cnt=1,scale_val=0,num_chains=10
          endif

          ; Uncomment to save output and make available for re-run
          ;openw, 1,'output_fin/'+name+'_best_fits_pso.txt',/append
          ;printf,1,'single:',val
          ;close,1

      end


; ****************************************************************************************************** ;
      'multi': begin
; ****************************************************************************************************** ;

        ; For more detailed comments see case 'single' above
        ; *************************************************** ;
        ; Re-define variables
          val=dblarr(15,9)
          go =1.0
          i =0.0
          yes=0.0

          if (best_val2[0] lt 99999 ) then begin 

              printmodels_v1,name,fittype=[2],ps='X',nstep=100,thin_val=1,teff=Teff,$
              dist=dist_val,amin=amin,seed=best_val2[0:13],cnt=1,scale_val=0,num_chains=20

              data = readfits(out_dir+'/'+name+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=1)

              if ( min(-2.*data(14,99,*)) le 5.0) then begin
                yes=1.0
                go = -1.0
              endif

              if ( min(-2.*data(14,99,*)) gt 100.0) then begin
                yes=0.0
                go = -1.0
              endif
          endif

          ; *************************************************** ;
          ; PSO Routine
          while (go gt 0.0) and ( i lt 9.0) do begin
             pso_min_example,name=name,teff=Teff,amin=amin,output= value,sequence = i
             val[*,i]=value
             
             if(( min(val[14,0:3]) le 1.00) and ( i eq 3)) then begin
                 go = -1.0
                 yes =1.0
                 best_val2 = value
             endif

             if(val[14,i] le 1.50) and ( i gt 3) then begin
                 go = -1.0
                 yes =1.0 
                 best_val2 = value
             endif

             if(val[14,0] ge 300.0) and  (val[14,1] ge 300.0) and (val[14,2] ge 300.0) and (val[14,3] ge 300.0) and (i eq 3) then begin
                 go = -1.0
                 yes =0.0 
             endif

             if (val[14,0] ge 100.0) and  (val[14,1] ge 100.0) and (val[14,2] ge 100.0) and (val[14,3] ge 100.0) and  (val[14,4] ge 100.0) and (val[14,5] ge 100.0)  and (i eq 5) then begin
                 go = -1.0
                 yes = 0.0
             endif
                 
             if ( i eq 8) then begin
                 abc = min(val[14,*],indx)
                 best_val2=val[*,indx]
                 yes=1.0
                 go = -1.0 
            endif
            i++
        endwhile

        ; *************************************************** ;
        ; MCMC Routine
          if (yes eq 1.0) then begin
              printmodels_v1,name,fittype=[2],ps='X',nstep=5000,thin_val=1,teff=Teff,$
              dist=dist_val,amin=amin,seed=best_val2[0:13],cnt=1,scale_val=0,num_chains=20
          endif

      end
      

; ****************************************************************************************************** ;
      'disk': begin
; ****************************************************************************************************** ;
        
        ; For more detailed comments see case 'single' above
        ; *************************************************** ;
        ; Re-define variables
        val = dblarr(12,6)
        go = 1.0
        i = 0.0
        yes=0.0
        
        if (best_val3[0] lt 99999 ) then begin
        
          printmodels_v1,name,fittype=[3],ps='X',nstep=100,thin_val=1,teff=Teff,$
            dist=dist_val,amin=amin,seed=best_val3[0:10],cnt=1,scale_val=0,num_chains=15
            
          data = readfits(out_dir+'/'+name+'_chn_mcmc_'+fit_name+'_part.fits',EXTEN_NO=1)
          
          if ( min(-2.*data(11,99,*)) le 5.0) then begin
            yes = 1.0
            go = -1.0
          endif
          
          if ( min(-2.*data(11,99,*)) gt 100.0) then begin
            yes = 0.0
            go = -1.0
          endif
          
        endif
        
        ; *************************************************** ;
        ; PSO Routine
        while (go gt 0.0) and ( i lt 6.0) do begin
          pso_min_example,name=name,teff=Teff,amin=amin,output= value,sequence = i

          val[*,i]=value

          if(val[11,i] le 1.00) and ( i le 2.0) then begin
            go = -1.0
            yes = 1.0
            best_val3 = value
          endif
          
          if(val[11,i] le 1.50) and ( i gt 2.0) then begin
            go = -1.0
            yes =1.0
            best_val3 = value
          endif
          
          if(val[11,0] ge 300.0) and  (val[11,1] ge 300.0) and (val[11,2] ge 300.0)  and (i eq 2) then begin
            go = -1.0
            yes =0.0
          endif
          
          if (val[11,0] ge 100.0) and  (val[11,1] ge 100.0) and (val[11,2] ge 100.0) and (val[11,3] ge 100.0) and  (val[11,4] ge 100.0) and (i eq 4) then begin
            go = -1.0
            yes = 0.0
          endif
          
          if ( i eq 5) then begin
            abc = min(val[11,*],indx)
            best_val3=val[*,indx]
            yes=1.0
            go = -1.0
          endif
          i++
        endwhile
        
        ; *************************************************** ;
        ; MCMC Routine
        if (yes eq 1.0) then begin
          printmodels_v1,name,fittype=[3],ps='X',nstep=5000,thin_val=1,teff=Teff,$
            dist=dist_val,amin=amin,seed=best_val3[0:10],cnt=1,scale_val=0,num_chains=15
        endif
        
        ; Uncomment to save output and make available for re-run
        ;openw, 1, 'output_fin/'+name+'_best_fits_pso.txt',/append
        ;printf,1,'disk:',val
        ;close,1
        
      end


endcase
return
end