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
;          3) 'single' - Single grain model
;          0) 'disk' - Continuous radial distribution of particles 
;          7) 'disk_mips' -  
;          4) 'ring' - Gaussian ring 
;          1) 'multi' - Two discrete grain populations, (T1,a1), (T2,a2)
;          6) 'multi_mips' - 
;
; KEYWORDS:
;   DO_OLD: Set to load old best fit values as starting point.
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
pro fits_v1, name=name,fittype=fittype,DO_OLD = do_old

; Load text file with names, Teff, amin and dist
fmt='a,f,f,f'
readcol,'input_files/input_param_file.txt',F=fmt,catalog_name,c_teff,c_amin,c_dist_val

; Match loaded data names with name of object being analyzed 
amin = .1
for i = 0, size(catalog_name,/n_elements)-1 do begin
   if (catalog_name[i] eq name) then begin
       Teff=c_teff[i]
       amin=c_amin[i]
       dist_val=c_dist_val[i]
    endif
endfor

; *************************************************** ;
; If name set as below, run test data

if (name eq 'Test1') or (name eq 'Test2') or (name eq 'Test3') or (name eq 'Test4') or (name eq 'Test5')  then begin
    Teff=8750
    amin =1.049
    dist_val=114.5
 endif

if (name eq 'Test6') or (name eq 'Test7')  then begin
    Teff=8750
    amin =1.049
    dist_val=114.5
 endif

; *************************************************** ;
; Set simulation best values

best_val1 = [99999,99999] ; for 'single'
best_val2 = [99999,99999] ; for 'disk', 'disk_mips'
best_val3 = [99999,99999] ; for 'ring'
best_val4 = [99999,99999] ; for 'multi', 'multi_mips'

; *************************************************** ;
; If keyword 'do_old' set, load saved best fit values

IF KEYWORD_SET(do_old) THEN BEGIN

fmt='a,f,f,f,f,f,f'
readcol,'input_files/input_param_file_single.txt',F=fmt,catalog_name,v1,v2,v3,v4,v5,v6

for i = 0, size(catalog_name,/n_elements)-1 do begin
   if (catalog_name[i] eq name) then begin
       best_val1=[v1[i],v2[i],v3[i],v4[i],v5[i],v6[i]]
   endif
endfor

fmt='a,f,f,f,f,f,f,f,f,f,f,f,f'
readcol,'input_files/input_param_file_multi.txt',F=fmt,catalog_name,v1,v2,v3,v4,v5,v6,va1,va2,va3,va4,va5,va6

for i = 0, size(catalog_name,/n_elements)-1 do begin
   if (catalog_name[i] eq name) then begin
       best_val4=[v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],va1[i],va2[i],va3[i],va4[i],va5[i],va6[i]]
   endif
endfor

fmt='a,f,f,f,f,f,f,f,f,f,f'
readcol,'input_files/input_param_file_disk.txt',F=fmt,catalog_name,v1,v2,v3,v4,v5,v6,va1,va2,va3,va4
for i = 0, size(catalog_name,/n_elements)-1 do begin
   if (catalog_name[i] eq name) then begin
       best_val2=[v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],va1[i],va2[i],va3[i],va4[i]]
   endif;
endfor

fmt='a,f,f,f,f,f,f,f,f,f'
readcol,'input_files/input_param_file_ring.txt',F=fmt,catalog_name,v1,v2,v3,v4,v5,v6,va1,va2,va3
for i = 0, size(catalog_name,/n_elements)-1 do begin
   if (catalog_name[i] eq name) then begin
       best_val3=[v1[i],v2[i],v3[i],v4[i],v5[i],v6[i],va1[i],va2[i],va3[i]]
 endif;
endfor
endif

; *************************************************** ;
; Begin fittype cases

; ****************************************************************************************************** ;
case fittype of
     'single': begin
; ****************************************************************************************************** ;      

          ; *************************************************** ;
          ; Short MCMC run for prior data
          val=dblarr(7,7) ; Define aray to hold fit values
          go = 1.0 ; Proceed to pso 
          i = 0.0 ; Begin iteration
          yes = 0.0 ; Do not procedd to MCMC simulation yet

          ; If prior best fit loaded, run a short MCMC to check for convergence
          if (best_val1[0] lt 99999 ) then begin

             ; Call MCMC routine with 100 steps of 10 chains
             printmodels_v1,name,fittype=[3],ps='X',nstep=100,thin_val=1,teff=Teff,$
             dist=dist_val,amin=amin,seed=best_val1[0:5],cnt=1,scale_val=0,num_chains=10

             ; Read in and store data from printmodels_v1 called just above
             data = readfits('output_v1/'+name+'_chn_mcmc_single_part.fits',EXTEN_NO=1)

             ; ???
             if ( min(-2.*data(6,99,*)) le 5.0) then begin
                yes = 1.0 ; Proceed to MCMC routine
                go = -1.0 ; Skip PSO routine
             endif

             ; ???
             if ( min(-2.*data(6,99,*)) gt 100.0) then begin
                yes=0.0 ; Forgo MCMC routine and exit to main
                go = -1.0 ; Skip PSO routine
             endif
          endif
          
          ; *************************************************** ;
          ; PSO Routine
          while (go gt 0.0) and ( i lt 7.0) do begin
             pso_min_example_s,name=name,teff=Teff,amin=amin,output= value,sequence = i ; Call PSO
             val[*,i]=value ; Save PSO's optimized values for the run

             ; If ??? proceed to MCMC             
             if(val[6,i] le 1.00) and ( i le 2.0) then begin
                 go = -1.0 ; Exit PSO routine
                 yes =1.0 ; Proceed to MCMC routine
                 best_val1 = value ; Update best value
             endif

             ; If ??? proceed to MCMC
             if(val[6,i] le 1.50) and ( i gt 2.0) then begin
                 go = -1.0 ; Exit PSO routine
                 yes =1.0 ; Proceed to MCMC routine
                 best_val1 = value ; Update best value
              endif

              ; If ??? exit to main
             if(val[6,0] ge 300.0) and  (val[6,1] ge 300.0) and (val[6,2] ge 300.0)  and (i eq 2) then begin
                 go = -1.0 ; Exit PSO routine
                 yes =0.0 ; Forgo MCMC routine and exit to main
             endif

             ; If ??? proceed to exit to main
             if (val[6,0] ge 100.0) and  (val[6,1] ge 100.0) and (val[6,2] ge 100.0) and (val[6,3] ge 100.0) and  (val[6,4] ge 100.0) and (val[6,5] ge 100.0) and (i eq 5) then begin
                     go = -1.0 ; Exit PSO routine
                     yes = 0.0 ; Forgo MCMC routine and exit to main
              endif
                 
              ; If at last run of PSO, proceed to MCMC
              if ( i eq 6) then begin
                     abc = min(val[6,*],indx)
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
              printmodels_v1,name,fittype=[3],ps='X',nstep=6000,thin_val=1,teff=Teff,$
              dist=dist_val,amin=amin,seed=best_val1[0:5],cnt=1,scale_val=0,num_chains=10
          endif

          ; Uncomment to save output and make available for re-run
          ;openw, 1,'output_fin/'+name+'_best_fits_pso.txt',/append
          ;printf,1,'single:',val
          ;close,1

      end

; ****************************************************************************************************** ;
     'disk': begin
; ****************************************************************************************************** ;
       
          ; *************************************************** ;
          ; Short MCMC run for prior data
          val=dblarr(11,6) ; Define aray to hold fit values
          go = 1.0 ; Proceed to pso
          i = 0.0 ; Begin iteration
          yes = 0.0 ; Do not procedd to MCMC simulation yet

          ; If prior best fit loaded, run a short MCMC to check for convergence
          if (best_val2[0] lt 99999 ) then begin 

            ; Call MCMC routine with 100 steps of 15 chains
            printmodels_v1,name,fittype=[0],ps='X',nstep=100,thin_val=1,teff=Teff,$
            dist=dist_val,amin=amin,seed=best_val2[0:9],cnt=1,scale_val=0,num_chains=15

            ; Read in and store data from printmodels_v1 called just above
            data = readfits('output_v1/'+name+'_chn_mcmc_disk_part.fits',EXTEN_NO=1)

            ; ???
             if ( min(-2.*data(10,99,*)) le 5.0) then begin
                yes=1.0 ; Proceed to MCMC routine
                go = -1.0 ; Skip PSO routine
             endif

             ; ???
             if ( min(-2.*data(10,99,*)) gt 100.0) then begin
                yes=0.0 ; Forgo MCMC routine and exit to main
                go = -1.0 ; Skip PSO routine
             endif
          endif

          ; *************************************************** ;
          ; PSO Routine
          while (go gt 0.0) and ( i lt 6.0) do begin
             pso_min_example_sd,name=name,teff=Teff,amin=amin,output= value,sequence = i ; Call PSO
             val[*,i] = value ; Save PSO's optimized values for the run
             
             ; If ??? proceed to MCMC
             if(val[10,i] le 1.00) and ( i le 2.0) then begin
                 go = -1.0 ; Exit PSO routine
                 yes = 1.0 ; Proceed to MCMC routine
                 best_val2 = value ; Update best value
             endif

             ; If ??? proceed to MCMC
             if(val[10,i] le 1.50) and ( i gt 2.0) then begin
                 go = -1.0 ; Exit PSO routine
                 yes =1.0 ; Proceed to MCMC routine 
                 best_val2 = value ; Update best value
             endif

                 ; If ??? exit to main
             if(val[10,0] ge 300.0) and  (val[10,1] ge 300.0) and (val[10,2] ge 300.0)  and (i eq 2) then begin
                 go = -1.0 ; Exit PSO routine
                 yes = 0.0 ; Forgo MCMC routine and exit to main 
             endif

             ; If ??? exit to main
             if (val[10,0] ge 100.0) and  (val[10,1] ge 100.0) and (val[10,2] ge 100.0) and (val[10,3] ge 100.0) and  (val[10,4] ge 100.0) and (i eq 4) then begin
                 go = -1.0 ; Exit PSO routine
                 yes = 0.0 ; Forgo MCMC routine and exit to main
             endif
                 
             ; If at last run of PSO, proceed to MCMC
             if ( i eq 5) then begin
                 abc = min(val[10,*],indx)
                 best_val2 = val[*,indx] ; Set best value to best of all 5 runs
                 yes = 1.0 ; Proceed to MCMC routine
                 go = -1.0 ; Exit PSO routine
             endif
             i++ ; iterate
         endwhile

         ; *************************************************** ;
         ; MCMC Routine
         if (yes eq 1.0) then begin
             ; Call MCMC routine with 5000 steps of 15 chains
             printmodels_v1,name,fittype=[0],ps='X',nstep=5000,thin_val=1,teff=Teff,$
             dist=dist_val,amin=amin,seed=best_val2[0:9],cnt=1,scale_val=0,num_chains=15
         endif

         ; Uncomment to save output and make available for re-run
         ;openw, 1, 'output_fin/'+name+'_best_fits_pso.txt',/append
         ;printf,1,'disk:',val
         ;close,1
          
      end

; ****************************************************************************************************** ;
      'disk_mips': begin
; ****************************************************************************************************** ;
         
        ; For more detailed comments see case 'single' above
        ; *************************************************** ;
        ; Re-define variables
          val = dblarr(11,6) 
          go = 1.0 
          i = 0.0 
          yes=0.0 

          if (best_val2[0] lt 99999 ) then begin 

              printmodels_v1,name,fittype=[7],ps='X',nstep=100,thin_val=1,teff=Teff,$
              dist=dist_val,amin=amin,seed=best_val2[0:9],cnt=1,scale_val=0,num_chains=15

              data = readfits('output_v1/'+name+'_chn_mcmc_disk_part.fits',EXTEN_NO=1)

              if ( min(-2.*data(10,99,*)) le 5.0) then begin
                yes = 1.0
                go = -1.0
              endif

              if ( min(-2.*data(10,99,*)) gt 100.0) then begin
                yes = 0.0
                go = -1.0
              endif

          endif

          ; *************************************************** ;
          ; PSO Routine
          while (go gt 0.0) and ( i lt 6.0) do begin
             pso_min_example_sd,name=name,teff=Teff,amin=amin,output= value,sequence = i
             val[*,i]=value
             
             if(val[10,i] le 1.00) and ( i le 2.0) then begin
                 go = -1.0
                 yes = 1.0 
                 best_val2 = value
             endif

             if(val[10,i] le 1.50) and ( i gt 2.0) then begin
                 go = -1.0
                 yes =1.0 
                 best_val2 = value
             endif

             if(val[10,0] ge 300.0) and  (val[10,1] ge 300.0) and (val[10,2] ge 300.0)  and (i eq 2) then begin
                 go = -1.0
                 yes =0.0 
             endif

             if (val[10,0] ge 100.0) and  (val[10,1] ge 100.0) and (val[10,2] ge 100.0) and (val[10,3] ge 100.0) and  (val[10,4] ge 100.0) and (i eq 4) then begin
                 go = -1.0
                 yes = 0.0
             endif
                 
             if ( i eq 5) then begin
                 abc = min(val[10,*],indx)
                 best_val2=val[*,indx]
                 yes=1.0
                 go = -1.0
             endif
             i++
         endwhile

         ; *************************************************** ;
         ; MCMC Routine            
         if (yes eq 1.0) then begin
             printmodels_v1,name,fittype=[7],ps='X',nstep=5000,thin_val=1,teff=Teff,$
             dist=dist_val,amin=amin,seed=best_val2[0:9],cnt=1,scale_val=0,num_chains=15
         endif

         ; Uncomment to save output and make available for re-run
         ;openw, 1, 'output_fin/'+name+'_best_fits_pso.txt',/append
         ;printf,1,'disk:',val
         ;close,1
          
      end

; ****************************************************************************************************** ;
      'ring': begin
; ****************************************************************************************************** ;
        
        ; For more detailed comments see case 'single' above
        ; *************************************************** ;
        ; Re-define variables
          val=dblarr(10,6)
          go =1.0
          i =0.0
          yes=0.0

          if (best_val3[0] lt 99999 ) then begin 

              printmodels_v1,name,fittype=[4],ps='X',nstep=100,thin_val=1,teff=Teff,$
              dist=dist_val,amin=amin,seed=best_val3[0:8],cnt=1,scale_val=0,num_chains=15

              data = readfits('output_v1/'+name+'_chn_mcmc_ring_part.fits',EXTEN_NO=1)

              if ( min(-2.*data(9,99,*)) le 5.0) then begin
                yes=1.0
                go = -1.0
             endif

             if ( min(-2.*data(9,99,*)) gt 100.0) then begin
                yes=0.0
                go = -1.0
             endif
          endif

          ; *************************************************** ;
          ; PSO Routine
          while (go gt 0.0) and ( i lt 6.0) do begin
             pso_min_example_cl,name=name,teff=Teff,amin=amin,output= value,sequence = i
             val[*,i]=value
             
             if(val[9,i] le 1.00) and ( i le 2.0) then begin
                 go = -1.0
                 yes =1.0 
                 best_val3 = value
             endif

             if(val[9,i] le 1.50) and ( i gt 2.0) then begin
                 go = -1.0
                 yes =1.0 
                 best_val3 = value
             endif

             if(val[9,0] ge 300.0) and  (val[9,1] ge 300.0) and (val[9,2] ge 300.0)  and (i eq 2.0) then begin
                 go = -1.0
                 yes =0.0 
             endif

             if (val[9,0] ge 100.0) and  (val[9,1] ge 100.0) and (val[9,2] ge 100.0) and (val[9,3] ge 100.0) and  (val[9,4] ge 100.0) and (i eq 4.0) then begin
                 go = -1.0
                 yes = 0.0
             endif
                 
             if ( i eq 5) then begin
                 abc = min(val[9,*],indx)
                 best_val3=val[*,indx]
                 go = -1.0
                 yes=1.0
             endif

             i++

          endwhile

          ; *************************************************** ;
          ; MCMC Routine
          if (yes eq 1.0) then begin
              printmodels_v1,name,fittype=[4],ps='X',nstep=5000,thin_val=1,teff=Teff,$
              dist=dist_val,amin=amin,seed=best_val3[0:8],cnt=1,scale_val=0,num_chains=15
          endif

          ; Uncomment to save output and make available for re-run
          ;openw, 1, 'output_fin/'+name+'_best_fits_pso.txt',/append
          ;printf,1,'ring:',val
          ;close,1

      end

; ****************************************************************************************************** ;
      'multi': begin
; ****************************************************************************************************** ;        

        ; For more detailed comments see case 'single' above
        ; *************************************************** ;
        ; Re-define variables
          val=dblarr(13,9)
          go =1.0
          i =0.0
          yes=0.0

          if (best_val4[0] lt 99999 ) then begin 

              printmodels_v1,name,fittype=[1],ps='X',nstep=100,thin_val=1,teff=Teff,$
              dist=dist_val,amin=amin,seed=best_val4[0:11],cnt=1,scale_val=0,num_chains=20

              data = readfits('output_v1/'+name+'_chn_mcmc_multi_part.fits',EXTEN_NO=1)

             if ( min(-2.*data(12,99,*)) le 5.0) then begin
                yes=1.0
                go = -1.0
             endif

             if ( min(-2.*data(12,99,*)) gt 100.0) then begin
                yes=0.0
                go = -1.0
             endif

          endif

          ; *************************************************** ;
          ; PSO Routine
          while (go gt 0.0) and ( i lt 9.0) do begin
             pso_min_example_m,name=name,teff=Teff,amin=amin,output= value,sequence = i
             val[*,i]=value
             
             if(( min(val[12,0:3]) le 1.00) and ( i eq 3)) then begin
                 go = -1.0
                 yes =1.0
                 best_val4 = value
              endif

 
             if(val[12,i] le 1.50) and ( i gt 3) then begin
                 go = -1.0
                 yes =1.0 
                 best_val4 = value
             endif

             if(val[12,0] ge 300.0) and  (val[12,1] ge 300.0) and (val[12,2] ge 300.0) and (val[12,3] ge 300.0) and (i eq 3) then begin
                 go = -1.0
                 yes =0.0 
             endif

             if (val[12,0] ge 100.0) and  (val[12,1] ge 100.0) and (val[12,2] ge 100.0) and (val[12,3] ge 100.0) and  (val[12,4] ge 100.0) and (val[12,5] ge 100.0)  and (i eq 5) then begin
                 go = -1.0
                 yes = 0.0
             endif
                 
             if ( i eq 8) then begin
                 abc = min(val[12,*],indx)
                 best_val4=val[*,indx]
                 yes=1.0
                 go = -1.0 
             endif
         i++
         endwhile

         ; *************************************************** ;
         ; MCMC Routine
         if (yes eq 1.0) then begin
             printmodels_v1,name,fittype=[1],ps='X',nstep=5000,thin_val=1,teff=Teff,$
             dist=dist_val,amin=amin,seed=best_val4[0:11],cnt=1,scale_val=0,num_chains=20
         endif

             
;;;;;;;;;; Figure out the second point for the mcmc ;;;;;;;;;;;;
;; First find non zero values of val[12,*], then find ones which are
;; less than 15 (chisq) , then find the ones with temp diff from best
;; val of at least 40 K in at least one temp. - Then choose one with
;; lowest chisq
;          d_val=where(val[12,*] gt 0 and val[12,*] le 15 and (( abs(min(val[0,*],val[6,*]) - min(best_val[0,*],best_val[6,*])) ge 40) or ( abs(max(val[0,*],val[6,*]) - max(best_val[0,*],best_val[6,*])) ge 40) ))
;          if (yes eq 1.0) and (d_val[0] ne -1) then begin
;               min2_val =  min(val[12,d_val],indx) 
;               best_val2 = val[*,d_val[indx]];
;
;               printmodels_v1,name,fittype=[1],ps='X',nstep=60000,thin_val=1,teff=Teff,$
;               dist=dist_val,amin=amin,seed=best_val2[0:11],cnt=2,scale_val=0
;          endif


         ; Uncomment to save output and make available for re-run
         ;openw, 1, 'output_fin/'+name+'_best_fits_pso.txt',/append
         ;printf,1,'multi:',val
         ;close,1
      end

; ****************************************************************************************************** ;
      'multi_mips': begin
; ****************************************************************************************************** ;

        ; For more detailed comments see case 'single' above
        ; *************************************************** ;
        ; Re-define variables
          val=dblarr(13,9)
          go =1.0
          i =0.0
          yes=0.0

          if (best_val4[0] lt 99999 ) then begin 

              printmodels_v1,name,fittype=[6],ps='X',nstep=100,thin_val=1,teff=Teff,$
              dist=dist_val,amin=amin,seed=best_val4[0:11],cnt=1,scale_val=0,num_chains=20

              data = readfits('output_v1/'+name+'_chn_mcmc_multi_part.fits',EXTEN_NO=1)

              if ( min(-2.*data(12,99,*)) le 5.0) then begin
                yes=1.0
                go = -1.0
              endif

              if ( min(-2.*data(12,99,*)) gt 100.0) then begin
                yes=0.0
                go = -1.0
              endif
          endif

          ; *************************************************** ;
          ; PSO Routine
          while (go gt 0.0) and ( i lt 9.0) do begin
             pso_min_example_m,name=name,teff=Teff,amin=amin,output= value,sequence = i
             val[*,i]=value
             
             if(( min(val[12,0:3]) le 1.00) and ( i eq 3)) then begin
                 go = -1.0
                 yes =1.0
                 best_val4 = value
             endif

             if(val[12,i] le 1.50) and ( i gt 3) then begin
                 go = -1.0
                 yes =1.0 
                 best_val4 = value
             endif

             if(val[12,0] ge 300.0) and  (val[12,1] ge 300.0) and (val[12,2] ge 300.0) and (val[12,3] ge 300.0) and (i eq 3) then begin
                 go = -1.0
                 yes =0.0 
             endif

             if (val[12,0] ge 100.0) and  (val[12,1] ge 100.0) and (val[12,2] ge 100.0) and (val[12,3] ge 100.0) and  (val[12,4] ge 100.0) and (val[12,5] ge 100.0)  and (i eq 5) then begin
                 go = -1.0
                 yes = 0.0
             endif
                 
             if ( i eq 8) then begin
                 abc = min(val[12,*],indx)
                 best_val4=val[*,indx]
                 yes=1.0
                 go = -1.0 
            endif
            i++
        endwhile

        ; *************************************************** ;
        ; MCMC Routine
          if (yes eq 1.0) then begin
              printmodels_v1,name,fittype=[6],ps='X',nstep=5000,thin_val=1,teff=Teff,$
              dist=dist_val,amin=amin,seed=best_val4[0:11],cnt=1,scale_val=0,num_chains=20
          endif

             
;;;;;;;;;; Figure out the second point for the mcmc ;;;;;;;;;;;;
;; First find non zero values of val[12,*], then find ones which are
;; less than 15 (chisq) , then find the ones with temp diff from best
;; val of at least 40 K in at least one temp. - Then choose one with
;; lowest chisq
;
;          d_val=where(val[12,*] gt 0 and val[12,*] le 15 and (( abs(min(val[0,*],val[6,*]) - min(best_val[0,*],best_val[6,*])) ge 40) or ( abs(max(val[0,*],val[6,*]) - max(best_val[0,*],best_val[6,*])) ge 40) ))
;          if (yes eq 1.0) and (d_val[0] ne -1) then begin
;               min2_val =  min(val[12,d_val],indx) 
;               best_val2 = val[*,d_val[indx]];
;
;               printmodels_v1,name,fittype=[1],ps='X',nstep=60000,thin_val=1,teff=Teff,$
;               dist=dist_val,amin=amin,seed=best_val2[0:11],cnt=2,scale_val=0
;
;          endif

          ; Uncomment to save output and make available for re-run
          ;openw, 1, 'output_fin/'+name+'_best_fits_pso.txt',/append
          ;printf,1,'multi:',val
          ;close,1
      end
endcase

; End Program
return
end