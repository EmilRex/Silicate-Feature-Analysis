function run_fit_prg

  compile_opt idl2

COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs
COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, effectiveTempArray, stellar_emit
 
restore, 'graintempdata.sav'
restore, 'qtables_withcrys2.sav' ; qastrosil, qolivine, qpyroxene                                                                     
              
   common mcmc_common, seed1  ;- save the random seed value
 common input_var,params,param_bnd,mcmc_res,num_parameter
  common trial_slct,current,curr_all,numb_shft
 common fnc_name,mips70_val,mips70_error,dof1,use_mips70      ;- save the function_name value   

  num_chains = params.num_chains
  data_base=[transpose(params.lambdafit),transpose(params.diskfit),transpose(params.errfit)]


 fxhmake,header1,data_base,/date
 header_lines = strarr(2)
 header_lines[0] = '/  MCMC Fit - Dust Model - Jang Condell et al. 2013, Mittal et al. 2013'
 header_lines[1] = '/  Spitzer IRS spectrum, Chen et al. 2013'
 sxaddhist, header_lines,header1
 file='output_fin_new/'+params.name_obj+'_chn_mcmc_disk_part.fits'
 FITS_WRITE,file,data_base,header1

  mcmc_res[*,0,*] = init_mcmc(num_chains,data_base)
  nsuccess = 0 & nfail = 0

  min_val_glb = max(mcmc_res[num_parameter,0,*],ind)
  elements = mcmc_res[*,0,ind]

  numb_shft =3.
  t = 1
  Lm = [0.0,0.0,0.0]
  p1 = 1./3.
  p2 = 1./3.
  p3 = 1./3.
  delm= [0.0,0.0,0.0]
  k = 1
  count = 0

  for i = 1, params.nstep - 1, 1 do begin

     min_val = max(mcmc_res[num_parameter,k-1,*],ind)           ; The min_val here is actaully log(exp(-chisq/2/dof))

     if (((-1.)*min_val) le ((-1.)*min_val_glb)) then begin
        elements = mcmc_res[*,k-1,ind]
        min_val_glb = min_val
     endif

;stop
       curr_all = mcmc_res[*,k-1,*]     
     for chn =0,num_chains-1,1 do begin

;        data= (mcmc_res[num_parameter,0:i-1,chn])
;        minVal = min(data)
;        maxVal = max(data)
;        medianVal = median(data,/even)

                                ; Find the quartiles.
 ;       Q1 = Median(data[Where(data LE medianVal, countlowerhalf)])
 ;       Q3= Median(data[Where(data GT medianVal, countupperhalf)])

           current = mcmc_res[0:num_parameter-1,k-1,chn]
           currentValue = mcmc_res[num_parameter,k-1,chn]

;;;;;;;;;;; After 500 steps, if a chain is not within 100*min_chi_sq,
;;;;;;;;;;; then it is replaced by the best chain

           if i>500 then begin
              if (-1.*max(mcmc_res[num_parameter,0:k-1,chn] ) > 100.*(-1.0)*min_val_glb)  then begin
                 current = elements[0:num_parameter-1]
                 currentValue = min_val_glb
        endif 



     endif

        prob = [p1,p2,p3]
        tmp  = max(multinom(1,prob,1,seed=seed19),in2)
        CR = 1.0  - prob[in2]
        Lm[in2] = Lm[in2]  +1

       u = randomu(seed1)
   
       ech_val = randomu(seed6,num_parameter)*(randomn(seed8,num_parameter,BINOMIAL=[1,.5]) -1.)
       ech_val  = (ech_val + 1.)

       tst_sub =  randomu(seed4,num_parameter)
       fac_num = ((CR-tst_sub) gt 0)*1D + ((CR-tst_sub) le 0)*0D

       while (n_elements(fac_num) lt 2) do begin
          tst_sub =  randomu(seed4,num_parameter)
          fac_num = ((CR-tst_sub) gt 0)*1D + ((CR-tst_sub) le 0)*0D
       endwhile


       diff_chng =  selectTrial(num_chains,fac_num,chn)       
       dprime = total(fac_num)

       gamd = 2.38/sqrt(2.*numb_shft*dprime) 
      if (i mod 5) le 0 then gamd =1.0

       chng_val = ech_val*gamd*fac_num*diff_chng + randomn(seed3,num_parameter,1)*param_bnd[1,*]*.0075
       trial = current + chng_val
;stop
       counts=where(trial gt param_bnd[1,*],number)
       if (number gt 0) then begin
          trial[counts]= trial[counts] - 2.*(trial[counts] - param_bnd[1,counts]) 
       endif
   
       counts=where(trial lt param_bnd[0,*],number)
       if (number gt 0) then begin
          trial[counts]= trial[counts] + 2.*(-trial[counts] + param_bnd[0,counts]) 
       endif

       counts=where(trial  gt param_bnd[1,*] or  trial  lt param_bnd[0,*],number)

       while(number ge 1) do begin
          ech_val[counts] = ech_val[counts]*.5
          chng_val = ech_val*gamd*fac_num*diff_chng + randomn(seed3,num_parameter,1)*param_bnd[1,*]*.0075
          trial[counts] = current[counts] + chng_val[counts]
          counts=where(trial gt param_bnd[1,*] or  trial lt param_bnd[0,*],number)
       end

          if(min(finite(trial))) eq 0 then begin
             trial = current
          end

       newValue = logTargetDistribution(trial,data_base)
      ;- determine acceptance probability via MH algorithm
      ;     alpha = exp(newValue - currentValue) * transitionRatio
     alpha = exp(newValue - currentValue)
  
     if u  lt alpha then begin  ;- new trial accepted
;        current = trial
;        currentValue = newValue
       mcmc_res[0:num_parameter-1,k,chn] = trial
       mcmc_res[num_parameter,k,chn] = newValue
        nsuccess++
    endif else begin 
       mcmc_res[0:num_parameter-1,k,chn] = current
       mcmc_res[num_parameter,k,chn] = currentValue
     nfail++                    ;- new trial rejected
   endelse

;    std_delm = stdev(mcmc_res[0:num_parameter-1,0:i,chn])
;    delm[in2] = delm[in2] + total(( mcmc_res[0:num_parameter-1,i,chn] -  mcmc_res[0:num_parameter-1,i-1,chn])^2/std_delm/std_delm)
;stop
 endfor


;     p1= t*num_chains*(delm[0]/Lm[0])/total(delm)
;     p2= t*num_chains*(delm[1]/Lm[1])/total(delm)
;     p3= t*num_chains*(delm[2]/Lm[2])/total(delm)
;     t = t+1
;       print,'Number: ',i
;       p1 =p1/(p1+p2+p3)
;       p2 =p2/(p1+p2+p3)
;    p3 =p3/(p1+p2+p3)
;stop

     k++
;;; Every 100th chain link, the mcmc chains (all chains) are stored in
;;; the fits file as a separate extension. Note that in each save set,
;;; the first row is the same as the last row of the preceding record
;;; - This is because the last record is needed to calculate whether a
;;;   new chain link should be accepted or not - Consequently, the
;;;                                              effective number of
;;;                                              links in a chain are .99*num_step
 
     if (k mod 100) le 0 then begin
;     save, file='output_fin_new/'+params.name_obj+'_chn_mcmc_multi_part.sav',i,nsuccess,nfail,num_chains,param_bnd,data_base,mcmc_res 
;       print,'Number: ',i
        count++
        file='output_fin_new/'+params.name_obj+'_chn_mcmc_disk_part.fits'
        FITS_OPEN,file,fcb,/append
        fxhmake,header1,mcmc_res,/date
        fxaddpar,header1,'total_num_chains',num_chains
        fxaddpar,header1,'total_num_steps',params.nstep
        fxaddpar,header1,'num_step_current',i
        fxaddpar,header1,'nfail',nfail
        fxaddpar,header1,'nsuccess',nsuccess
        fxaddpar,header1,'dof',dof1
        fxaddpar,header1,'global_min_val_red_chisq',min_val_glb*(-2.)

        nm = 'chain_part_' + strtrim(count,1)
        FITS_WRITE,fcb,mcmc_res,extname=nm,extver=1
        FITS_CLOSE,FCB
        mcmc_res[*,0,*] = mcmc_res[*,99,*]

        k = 1

     endif

  endfor

     min_val = max(mcmc_res[num_parameter,99,*],ind)           ; The min_val here is actaully log(exp(-chisq/2/dof))

     if (((-1.)*min_val) le ((-1.)*min_val_glb)) then begin
        elements = mcmc_res[*,99,ind]
        min_val_glb = min_val
     endif

        file='output_fin_new/'+params.name_obj+'_chn_mcmc_disk_part.fits'
        FITS_OPEN,file,fcb,/append
        fxhmake,header1,elements,/date
        fxaddpar,header1,'total_num_chains',num_chains
        fxaddpar,header1,'total_num_steps',params.nstep
        fxaddpar,header1,'num_step_current',i
        fxaddpar,header1,'nfail',nfail
        fxaddpar,header1,'nsuccess',nsuccess
        fxaddpar,header1,'dof',dof1 
        fxaddpar,header1,'global_min_val_red_chisq',min_val_glb*(-2.)

        nm = 'final_output'
        FITS_WRITE,fcb,elements,extname=nm,extver=1
        FITS_CLOSE,FCB
 
;     save, file='output_fin_new/'+params.name_obj+'_chn_mcmc_multi_fin.sav',i,nsuccess,nfail,num_chains,param_bnd,data_base,mcmc_res 

return,1
end

function selectTrial,num_chains,fac_num,index_no

  common trial_slct,current,curr_all,numb_shft

indx = round(randomu(seed_val,2*numb_shft,[n_elements(current)],/double)*(num_chains-3))
outpt = dblarr(n_elements(current))

for i =0,n_elements(current)-1,1 do begin
   if (fac_num[i] ge 1) then begin
      indx[*,i]  = indx[*,i] + ((index_no-indx[*,i]) eq 0)*1D
      if (indx[0,i] eq indx[1,i]) then  indx[1,i] = indx[1,i]+1
      if (indx[2,i] eq indx[3,i]) then  indx[3,i] = indx[3,i]+1
      if (indx[4,i] eq  indx[5,i]) then  indx[4,i] = indx[5,i]+1
      indx[*,i]  = indx[*,i] + ((index_no-indx[*,i]) eq 0)*1D

      outpt[i] =  curr_all[i,0,indx[0]]+curr_all[i,0,indx[2]]+curr_all[i,0,indx[4]]-curr_all[i,0,indx[1]]-curr_all[i,0,indx[3]]-curr_all[i,0,indx[5]]

   endif

endfor


return, outpt
end


function init_mcmc,num_chains,data_base

common input_var,params,param_bnd,mcmc_res,num_parameter

     outpt = dblarr(num_parameter+1,1,num_chains)
     u = randomu(sd2,num_parameter,1,num_chains)

     for i=0,num_chains-1,1 do begin
        par = params.seed*(1.+ randomu(sd3,num_parameter)*.005 )
        outpt[0:num_parameter-1,0,i] = par
        outpt[num_parameter,0,i] = logTargetDistribution(par,data_base)
     endfor

return,outpt
end

function logTargetDistribution,link,result

common fnc_name,mips70_val,mips70_error,dof1,use_mips70 

COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs
COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, effectiveTempArray, stellar_emit

;COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs
;COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, effectiveTempArray, stellar_emit

;stop
;'disk model': begin
;link=[113.575,5.38780,19.2318,0.263435,0.0392570,0.299372,333.999,26.3491,18.2289,0.862035,0.00286000,0.419083]
 
                                                                                                                                 
          link[0]=10^(link[0])
          link[6]=10^(link[6])
        
          if (use_mips70 gt 0) then begin

             spectra1 = modelsinglespectrum([[result[0,*]],[71.42]], link) ;                                                                         
             mips70 = spectra1[n_elements(result[0,*])]
             spectra= spectra1[0:n_elements(result[0,*])-1]

             chisq = TOTAL ( ((result[1,*]-spectra)^2.0)/((.05*result[2,*])^2.0+(result[2,*])^2.0))


                                ; The assumption here is that since
                                ; there can exist an additional dust
                                ; component besides the one modeled
                                ; here contributing principally only
                                ; to the mips70 flux, all models which
                                ; either fit mips70 flux or
                                ; under-represent the flux are
                                ; considered equally likely

                chisq2 = ((mips70_val-mips70)^2.0)/((mips70_error)^2.0)
                 like_func  = -(chisq+chisq2)/(2.0*dof1) 
 ;  Log of likelihood func - where likelihood funct is exp(-chisq/2 ) 
 ;This likelihood is valid under the assumption of gaussian errors
          endif else begin

             spectra =  modelsinglespectrum(result[0,*], link)
             chisq = TOTAL ( ((result[1,*]-spectra)^2.0)/((.05*result[2,*])^2.0+(result[2,*])^2.0))
              like_func = -(chisq)/(2.0*dof1)

           endelse
   
         link[6]=alog10(link[6])
          link[0]=alog10(link[0])


;  print,'like_finc ',-2.*like_func
  return, like_func
end


pro mcmc_sd_mips,par_bound=param_bnd_inp,parameter=params_inp

  compile_opt idl2

common fnc_name,mips70_val,mips70_error,dof1,use_mips70      ;- save the function_name value   
common input_var,params,param_bnd,mcmc_res,num_parameter

params = params_inp
param_bnd = param_bnd_inp
num_parameter = n_elements(params.seed)
 
err_chk=where(params.errfit le .01*params.diskfit)
if (max(err_chk) ge 0) then params.errfit[err_chk]=.01*params.diskfit[err_chk]


mcmc_res = dblarr(num_parameter+1,100.,params.num_chains) 

dof1 = n_elements(params.lambdafit) - num_parameter

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
;; Checking to see if mips70 value is available and if yes, adding a
;; degree of freedom (dof1) 

use_mips70=0
if(params.mips70_val gt 0.0) then begin 
   dof1 = dof1 + 1.0 ; Add 1 for mips 70 value 
   use_mips70 = 1
endif

mips70_val=params.mips70_val ; This is the photosphere subtracted MIPS70 value
mips70_error=params.mips70_error

done_yes =  run_fit_prg() ; Call the MCMC fitting routine 

end




