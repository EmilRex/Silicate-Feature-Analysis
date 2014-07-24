; +
; NAME
;  mcmc_general
;
; PURPOSE
;  Perform MCMC routine for the given model
;
; INPUTS
;  PARAM_BND_INP: Structure containing parameter bounds
;  PARAMS_INP: Structure containing parameters and data
;
; KEYWORDS
;   NONE
;
; OUTPUTS
;   FITS files containing steps and results of MCMC simulation
;
; AUTHORS
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
; MODIFICATION HISTORY
;  Written by TM (June 2013) as mcmc_m_mips.pro
;  Organized and commented by EC (6/26/2014)
;  Generalized and renamed by EC (7/14/14)
; -
; *************************************************** ;


; ****************************************************************************************************** ;
function run_fit_prg
; ****************************************************************************************************** ;

compile_opt idl2

COMMON file_path, in_dir, out_dir, fit_name, object_name

; Create global variables relating to silicate features
COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs
COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, effectiveTempArray, stellar_emit

; Load data into some of the global parameters above
restore, 'graintempdata.sav'
restore, 'qtables_withcrys2.sav' ; qastrosil, qolivine, qpyroxene                                                                     

; Create more global variables
COMMON mcmc_common, seed1  ;- save the random seed value
COMMON input_var,params,param_bnd,mcmc_res,num_parameter
COMMON trial_slct,current,curr_all,numb_shft
COMMON fnc_name,dof1      ;- save the function_name value   

; *************************************************** ;
; Extract data from params and initialize result storage
num_chains = params.num_chains
data_base=[transpose(params.lambdafit),transpose(params.diskfit),transpose(params.errfit)]

; Create FITS file for MCMC chain storage
fxhmake,header1,data_base,/date
header_lines = strarr(2)
header_lines[0] = '/  MCMC Fit - Dust Model - Jang Condell et al. 2013, Mittal et al. 2013'
header_lines[1] = '/  Spitzer IRS spectrum, Chen et al. 2013'
sxaddhist, header_lines,header1
file=out_dir+'/'+params.name_obj+'_chn_mcmc_'+fit_name+'_part.fits'
FITS_WRITE,file,data_base,header1

; Initialize MCMC chain storage
mcmc_res[*,0,*] = init_mcmc(num_chains,data_base)

; Reset success/fail counters
nsuccess = 0 & nfail = 0

; Set global min based on initial
min_val_glb = max(mcmc_res[num_parameter,0,*],ind)
elements = mcmc_res[*,0,ind]

; Initialize MCMC parameters
numb_shft =3.
t = 1
Lm = [0.0,0.0,0.0]
p1 = 1./3.
p2 = 1./3.
p3 = 1./3.
delm= [0.0,0.0,0.0]
k = 1
count = 0

; *************************************************** ;
; Begin main loop
; data stored in mcmc_res: [parameters+fcn_val,step,chain]


; *************************************************** ;
; For each step
for i = 1, params.nstep - 1, 1 do begin ; iterate through nsteps
; *************************************************** ;

  min_val = max(mcmc_res[num_parameter,k-1,*],ind) ; The min_val here is actaully log(exp(-chisq/2/dof))

  ; If new min is better than global min, replace
  if (((-1.)*min_val) le ((-1.)*min_val_glb)) then begin
    elements = mcmc_res[*,k-1,ind]
    min_val_glb = min_val
  endif

  ; Pick current step
  curr_all = mcmc_res[*,k-1,*]  
  
  ; *************************************************** ;
  ; For each chain   
  for chn =0,num_chains-1,1 do begin ; iterate through chains
  ; *************************************************** ;

    ; Statistics
    ;data= (mcmc_res[num_parameter,0:i-1,chn])
    ;minVal = min(data)
    ;maxVal = max(data)
    ;medianVal = median(data,/even)
    ;Find the quartiles.
    ;Q1 = Median(data[Where(data LE medianVal, countlowerhalf)])
    ;Q3= Median(data[Where(data GT medianVal, countupperhalf)])
  
    ; Select data for given step and chain
    current = mcmc_res[0:num_parameter-1,k-1,chn]
    currentValue = mcmc_res[num_parameter,k-1,chn]
  
    ; After 500 steps, if a chain is not within 100*min_chi_sq, then it is replaced by the best chain
    if i>500 then begin
      if (-1.*max(mcmc_res[num_parameter,0:k-1,chn] ) > 100.*(-1.0)*min_val_glb)  then begin
        current = elements[0:num_parameter-1]
        currentValue = min_val_glb
      endif 
    endif

    ; *************************************************** ;
    ; Random walk... take a stroll
    
    ; Simulate probabilities with multinom in /pro
    prob = [p1,p2,p3]
    tmp  = max(multinom(1,prob,1,seed=seed19),in2) ; SIMULATE MULTINOMIAL RANDOM VARIABLES
    CR = 1.0  - prob[in2]
    Lm[in2] = Lm[in2]  +1

    ; 
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
  
    ; *************************************************** ;
    ; Evaluate and store new position
    ;print, 'Iteration: '+strcompress(string(k))
    
    newValue = logTargetDistribution(trial,data_base)
    
    ; Determine acceptance probability via MH algorithm
    ;alpha = exp(newValue - currentValue) * transitionRatio
    alpha = exp(newValue - currentValue) 
    
    if u lt alpha then begin ; New trial accepted
      ;current = trial
      ;currentValue = newValue
      mcmc_res[0:num_parameter-1,k,chn] = trial
      mcmc_res[num_parameter,k,chn] = newValue
      nsuccess++
    endif else begin ; New trial rejected
      mcmc_res[0:num_parameter-1,k,chn] = current
      mcmc_res[num_parameter,k,chn] = currentValue
      nfail++ 
    endelse
  
    ;std_delm = stdev(mcmc_res[0:num_parameter-1,0:i,chn])
    ;delm[in2] = delm[in2] + total(( mcmc_res[0:num_parameter-1,i,chn] -  mcmc_res[0:num_parameter-1,i-1,chn])^2/std_delm/std_delm)

  endfor

  ;p1= t*num_chains*(delm[0]/Lm[0])/total(delm)
  ;p2= t*num_chains*(delm[1]/Lm[1])/total(delm)
  ;p3= t*num_chains*(delm[2]/Lm[2])/total(delm)
  ;t = t+1
  ;print,'Number: ',i
  ;p1 =p1/(p1+p2+p3)
  ;p2 =p2/(p1+p2+p3)
  ;p3 =p3/(p1+p2+p3)

  k++ ; iterate 100th step counter

  ; *************************************************** ;  
  ; Every 100th chain link, the mcmc chains (all chains) are stored in
  ; the fits file as a separate extension. Note that in each save set,
  ; the first row is the same as the last row of the preceding record.
  ; 
  ; This is because the last record is needed to calculate whether a
  ; new chain link should be accepted or not. 
  ; 
  ; Consequently, the effective number of links in a chain are .99*num_step.
 
  if (k mod 100) le 0 then begin
  ;save, file='output_fin_new/'+params.name_obj+'_chn_mcmc_multi_part.sav',i,nsuccess,nfail,num_chains,param_bnd,data_base,mcmc_res 
  ;print,'Number: ',i
  count++
  file=out_dir+'/'+params.name_obj+'_chn_mcmc_'+fit_name+'_part.fits'
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

  k = 1 ; Reset 100th step counter

  endif

endfor

; *************************************************** ;
; End main loop

; Edited two lines below!!! 99 was k-1 (which is then 0)
;k =1 ; Reset k
min_val = max(mcmc_res[num_parameter,99,*],ind) ; The min_val here is actaully log(exp(-chisq/2/dof))

; Update best value if needed
if (((-1.)*min_val) le ((-1.)*min_val_glb)) then begin
  elements = mcmc_res[*,99,ind]
  min_val_glb = min_val
endif

; Write final best value as last extension of associated FITS file
file=out_dir+'/'+params.name_obj+'_chn_mcmc_'+fit_name+'_part.fits'
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
 
;save, file='output_fin_new/'+params.name_obj+'_chn_mcmc_multi_fin.sav',i,nsuccess,nfail,num_chains,param_bnd,data_base,mcmc_res 

; End of main script -> exit out
return,1
end


; ****************************************************************************************************** ;
function selectTrial,num_chains,fac_num,index_no
; ****************************************************************************************************** ;
  

COMMON trial_slct,current,curr_all,numb_shft

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


; ****************************************************************************************************** ;
function init_mcmc,num_chains,data_base
; ****************************************************************************************************** ;


COMMON input_var,params,param_bnd,mcmc_res,num_parameter

outpt = dblarr(num_parameter+1,1,num_chains)
u = randomu(sd2,num_parameter,1,num_chains)

for i=0,num_chains-1,1 do begin
  par = params.seed*(1.+ randomu(sd3,num_parameter)*.005 )
  outpt[0:num_parameter-1,0,i] = par
  outpt[num_parameter,0,i] = logTargetDistribution(par,data_base)
endfor

return,outpt
end


; ****************************************************************************************************** ;
function logTargetDistribution,link,result
; ****************************************************************************************************** ;
; Take current 'link' in current 'chain', model the spectrum and compare to actual 'result'

COMMON fnc_name,dof1 
COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs
COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, effectiveTempArray, stellar_emit
COMMON file_path, in_dir, out_dir, fit_name, object_name


; *************************************************** ;
; ONE GRAIN MODEL
;
;temp1, grain1, scale1, olivine/pyroxene, crystalline, enstatite/fosterite,

IF (fit_name eq 'single') THEN BEGIN

  link[2]=10^(link[2])
  
  spectra = modelonegrain(result[0,*], link)
  chisq = TOTAL ( ((result[1,*]-spectra)^2.0)/((.05*result[2,*])^2.0+(result[2,*])^2.0))
  like_func = -(chisq)/(2.0*dof1)
  
  link[2]=alog10(link[2])
  
ENDIF

; *************************************************** ;
; TWO GRAIN MODEL
;
;temp1, grain1, scale1, olivine/pyroxene, crystalline, enstatite/fosterite,  
;temp2, grain2, scale2, olivine/pyroxene2, crystalline2, enstatite/fosterite2 

IF (fit_name eq 'multi_mips') THEN BEGIN

; un'log' value
  link[2]=10^(link[2])
  link[8]=10^(link[8])

  spectra = modeltwograin(result[0,*], link)
  chisq = TOTAL ( ((result[1,*]-spectra)^2.0)/((.05*result[2,*])^2.0+(result[2,*])^2.0))
  like_func = -(chisq)/(2.0*dof1)

; re'log' value
  link[2]=alog10(link[2])
  link[8]=alog10(link[8])

ENDIF

; *************************************************** ;
; CONTINUOUS DISK MODEL
;
;   rin  = 10^fitparam[0]
;   rout = exp(fitparam[1])*rin
;   sigmalaw = fitparam[2]
;   amin = fitparam[3]
;   amax = exp(fitparam[4])*amin
;   grainlaw = fitparam[5]
;   diskmass = 10^(fitparam[6]) * massscale
;   fo = fitparam[7]
;   fc = fitparam[8]
;   ff = fitparam[9]

IF (fit_name eq 'disk_mips') THEN BEGIN
  
  link[0]=10^(link[0])
  link[6]=10^(link[6])

  spectra =  modelsinglespectrum(result[0,*], link)
  chisq = TOTAL ( ((result[1,*]-spectra)^2.0)/((.05*result[2,*])^2.0+(result[2,*])^2.0))
  like_func = -(chisq)/(2.0*dof1)
  
  link[6]=alog10(link[6])
  link[0]=alog10(link[0])
  
ENDIF

; return the likelyhood value (i.e. fitness of link)
return, like_func
end


; ****************************************************************************************************** ;
pro mcmc_general,par_bound=param_bnd_inp,parameter=params_inp
; ****************************************************************************************************** ;


compile_opt idl2

; Create global variables relating to input data
common fnc_name,dof1      ;- save the function_name value   
common input_var,params,param_bnd,mcmc_res,num_parameter

; Rename data structure and extract number of variables
params = params_inp
param_bnd = param_bnd_inp
num_parameter = n_elements(params.seed)

; Check errors - ???
err_chk=where(params.errfit le .01*params.diskfit)
if (max(err_chk) ge 0) then params.errfit[err_chk]=.01*params.diskfit[err_chk]

; Define chain storage array; variables stored = free parameters plus function value (chisq)
mcmc_res = dblarr(num_parameter+1,100.,params.num_chains) 

; Define Degrees of freedom -> (# of data points(varies >10^2) - # of parameters(12) )
dof1 = n_elements(params.lambdafit) - num_parameter

; *************************************************** ;
; Checking to see if mips70 value is available and if yes, adding a
; degree of freedom (dof1) 

; -------------------------------------------------- ;
; Below determines whether to include mips... delete when changing code
;use_mips70=0
;if(params.mips70_val gt 0.0) then begin 
;   dof1 = dof1 + 1.0 ; Add 1 for mips 70 value 
;   use_mips70 = 1
;endif
;mips70_val=params.mips70_val ; This is the photosphere subtracted MIPS70 value
;mips70_error=params.mips70_error
; -------------------------------------------------- ;

; Setup complete -> run program
done_yes =  run_fit_prg() ; Call the MCMC fitting routine 

; End program
end