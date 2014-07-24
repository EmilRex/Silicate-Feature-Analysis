; +
; NAME:
;  modelfit_v1
;
; PURPOSE:
;  Declare parameter bounds for each model and pass to MCMC routine
;
; INPUTS:
;  PARAMS:
;  FITTYPE: Model to be fitted (in order of appearance below)
;          not set) 'disk' - Continuous radial distribution of particles 
;          1) 'multi' - Two discrete grain populations, (T1,a1), (T2,a2)
;          3) 'single' - Single grain model
;          4) 'ring' - Gaussian ring
;          5) Two ring
;          6) 'multi_mips' -
;          7) 'disk_mips' -  
;          
; KEYWORDS:
;   NOCRYS = Do not use crystals in model
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
;  Written by TM (June 2013) as modelfit_new.pro
;  Organized and commented by EC (6/26/2014)
; -
; *************************************************** ;
pro modelfit_v1,param_struc=params,fittype=fittype

; Create global variables relating to silicate features
COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs
COMMON stellarprops, temptable, folivine, effectiveTemp, lambdastar, fluxstar
COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, effectiveTempArray, stellar_emit

; Load data into some of the global parameters above
restore, 'graintempdata.sav'
restore, 'qtables_withcrys2.sav'

;Set physical parameters: ASSUMED BULK DENSITY OF GRAINS (SUBJECT TO INTERPRETATION)
rho_s = 3.3 
AU_in_cm = 1.496e13

; Fill parameters from params structure
dist=params.distval
Teff=params.teff_val
print, 'T_star=', Teff
effectiveTemp = Teff
amin = params.amin_val


; *************************************************** ;
; Find the right grain model for calculating temperatures

; Get file names and store in grainfiles
cmd = 'ls modelgrids/Teff*grains.sav'
spawn, cmd, grainfiles

; Strip temps from names
strbeg = strpos(grainfiles, 'Teff')+3
strend = strpos(grainfiles, 'grains')
tarray = fltarr(n_elements(grainfiles))

for i=0,n_elements(tarray)-1 do begin
   tarray[i] = float(strmid(grainfiles[i], strbeg[i]+1, strend[i]-strbeg[i]-1))
endfor

; reorder arrays
ii = sort(tarray)
grainfiles = grainfiles[ii]
tarray = tarray[ii]

; Find and load correct value
kuruczindex = interpol(findgen(n_elements(tarray)),tarray,effectiveTemp)
ki = round(kuruczindex) < (n_elements(tarray)-1) > 0
; next command will retrieve temptable, folivine, Teff
; as generated in generateallgrid.pro
; folivine = [0.0, 1.0]
restore, grainfiles[ki]
   
; *************************************************** ;
; Begin model scenarios

if not keyword_set(fittype) then begin

; Work in progress
print, 'Fit type not currenly available...'
print, 'Please check back later.'
return

; CASE 1: no fittype selected
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

  ; Define physically acceptible parameter bounds
    param_bnd = dblarr(3,10)
    param_bnd[*,0]= [1, 5.0,0] ; rin
    param_bnd[*,1]= [0.0,10.0,0] ; r range
    param_bnd[*,2]= [-5,5 ,0] ; no limits on alpha
    param_bnd[*,3]= [amin, 30.0,0]; amin - blowout radius as lower bound
    param_bnd[*,4]= [0.047, 7.0,0] ; a range
    param_bnd[*,5]= [1.5, 5.5,0] ; no limits on pindex
    param_bnd[*,6]= [16.5, 23.5,0] ; no limits on mass
    param_bnd[*,7]= [0.0, 1.0,0] ; amorph composition
    param_bnd[*,8]= [0.0, 1.0,0] ; crystalline fraction
	
	  if keyword_set(nocrys) then begin
		  param_bnd[*,8]= [0.0, 1.0,1]
      initparams[8] = 0.0
  	endif
  	
	  param_bnd[*,9]= [0.0, 1.0,0] ; crystalline composition

;  scale_val=[[.030,.030,.0220,0.80,.030,.0220, .0220, .0150, .015, .0150],$
;  [.050,.050,.0350,1.60,.050,.0350, .0420, .0350, .0325, .0150],$
;  [.080,.080,.0350,1.80,.080,.0350, .0620, .0350, .0325, .0150]]

  ; Run MCMC simulation with given bounds
  mcmc_sd,par_bound=param_bnd,parameter=params
  ; mcmc__define_sd,'modelsinglespectrum',par_bound=param_bnd,parameter=params,scale_val=scale_val[*,params.scale_val]

endif else begin
  case fittype of 
  
  ; *************************************************** ;
  ; TWO COMPONENT GRAIN MODEL
  ; params = [T1, a1, scale1, fo1, fc1, ff1,
  ;           T2, a2, scale2, fo2, fc2, ff2]
  1: begin 
  
    ; Work in progress
    print, 'Fit type not currenly available...'
    print, 'Please check back later.'
    return  
  
  	param_bnd = dblarr(3,12)
		param_bnd[*,0]= [30.0, 300.0,0] ; Temperatures limited to reasonable values
		param_bnd[*,1]= [amin, 100.0,0] ; blowout radius as lower bound
		param_bnd[*,2]= [16.5, 23.5,0] ; positive scale factors
		param_bnd[*,3]= [0, 1.0,0] ; olivine/pyroxene composition 
		param_bnd[*,4]= [0, 1.0,0] ; crystalline fraction
		param_bnd[*,5]= [0, 1.0,0] ; forsterite/enstatite composition 
 		param_bnd[*,6]= [100.0, 1000.0,0] ; Temperatures limited to reasonable values
		param_bnd[*,7]= [amin, 100.0,0] ; blowout radius as lower bound
		param_bnd[*,8]= [16.5, 23.5,0] ; positive scale factors
		param_bnd[*,9]= [0, 1.0,0] ; olivine/pyroxene composition 
		param_bnd[*,10]= [0, 1.0,0] ; crystalline fraction
		param_bnd[*,11]= [0, 1.0,0] ; forsterite/enstatite composition 

    if keyword_set(nocrys) then begin
	   	param_bnd[*,4]= [0, 1.0,1]
      param_bnd[*,5]= [0, 1.0,1]
      param_bnd[*,10]= [0, 1.0,1]
      param_bnd[*,11]= [0, 1.0,1]
      initparams[4] = 0.0
	    initparams[10] = 0.0
    endif

;   scale_val=[1.010,.70, .0220, .0150, .0150, .0150,  1.50, .530,  .0120, .0150, .0150, .0150]
;   mcmc__define_m,'modeltwograin',par_bound=param_bnd,parameter=params,scale_val=scale_val
    mcmc_m,par_bound=param_bnd,parameter=params

    end

  ; *************************************************** ;
  ; SINGLE GRAIN MODEL
  ; params = [Temperature, agrain, folivine, fcrystal, scale]
  3: begin
    
    print, 'Fitting single grain model...'
    
	  param_bnd = dblarr(3,6)
 		;param_bnd[*,0]= [10.0, 1000.0,0] ; Temperatures limited to reasonable values
		param_bnd[*,0]= [0.1, 500.0,0] ; New distance bounds
		param_bnd[*,1]= [amin, 30.0,0] ; .1 as lower bound
		param_bnd[*,2]= [16.5, 23.5,0] ; positive scale factors
		param_bnd[*,3]= [0, 1.0,0] ; olivine/pyroxene composition 
		param_bnd[*,4]= [0, 1.0,0] ; crystalline fraction
		param_bnd[*,5]= [0, 1.0,0] ; forsterite/enstatite composition 


    mcmc_general,par_bound=param_bnd,parameter=params
         
  end

  ; *************************************************** ;
  ; RING MODEL
  ;        radius = fitparam[0]
  ;        width = fitparam[1]*radius
  ;        amin = fitparam[2]
  ;        amax = exp(fitparam[3])*amin
  ;        grainlaw = fitparam[4]
  ;        fo = fitparam[6]
  ;        fc = fitparam[7]
  ;        ff = fitparam[8]
  ; radius, width/radius, amin, log(amax/amin), p, mass, fo, fc, f_fost/enst
  4: begin

    ; Work in progress
    print, 'Fit type not currenly available...'
    print, 'Please check back later.'
    return
  
 	  param_bnd = dblarr(3,9)
 		param_bnd[*,0]= [1.0, 5.0,0] ; radius
 		param_bnd[*,1]= [0, 2.0,0] ; width: as a fraction of radius
 		param_bnd[*,2]= [amin, 30.0,0] ; amin
 		param_bnd[*,3]= [0.047,7.0,0] ; a range
 		param_bnd[*,4]= [1.5,4.5,0] ; p

   ; composition
		param_bnd[*,5]= [16.5, 23.5,0] ; positive scale factors  - MASS
		param_bnd[*,6]= [0, 1.0,0] ; olivine/pyroxene composition 
		param_bnd[*,7]= [0, 1.0,0] ; crystalline fraction
		param_bnd[*,8]= [0, 1.0,0] ; forsterite/enstatite composition 
    
    ;scale_val=[[.030,.0120,0.50,.010,.0120, .0120, .0250, .0225, .0225],$
    ;[.050,.0350,1.60,.040,.0350, .0420, .0350, .0325, .0150],$
    ;[.0800,.0350,1.80,.060,.0350, .0620, .0350, .0325, .0150]]
 
    ;mcmc__define_cl,'modelring',par_bound=param_bnd,parameter=params,scale_val=scale_val[*,params.scale_val]
    mcmc_cl,par_bound=param_bnd,parameter=params

  end
    
  ; *************************************************** ;
  ; TWO INDEPENDENT RING MODEL
  ;        radius = fitparam[0]
  ;        width = fitparam[1]*radius
  ;        amin = fitparam[2]
  ;        amax = exp(fitparam[3])*amin
  ;        grainlaw = fitparam[4]
  ;        fo = fitparam[6]
  ;        fc = fitparam[7]
  ;        ff = fitparam[8]
  ; radius, width/radius, amin, log(amax/amin), p, mass, fo, fc, f_fost/enst
  5: begin
 
    ; Work in progress
    print, 'Fit type not currenly available...'
    print, 'Please check back later.'
    return
    
    param_bnd = dblarr(3,18)
 		param_bnd[*,0]= [1.0, 5.0,0] ; radius
		param_bnd[*,1]= [0, 2.0,0] ; width: as a fraction of radius
 		param_bnd[*,2]= [amin, 30.0,0] ; amin
 		param_bnd[*,3]= [0.047,7.0,0] ; a range 
 		param_bnd[*,4]= [1.5,4.5,0] ; p

    ; composition
		param_bnd[*,5]= [16.5, 23.5,0] ; positive scale factors  - MASS
		param_bnd[*,6]= [0, 1.0,0] ; olivine/pyroxene composition 
		param_bnd[*,7]= [0, 1.0,0] ; crystalline fraction
		param_bnd[*,8]= [0, 1.0,0] ; forsterite/enstatite composition 
 		param_bnd[*,9]= [1.0, 5.0,0] ; radius
 		param_bnd[*,10]= [0, 2.0,0] ; width: as a fraction of radius
 		param_bnd[*,11]= [amin, 30.0,0] ; amin
 		param_bnd[*,12]= [0.047,7.0,0] ; a range
 		param_bnd[*,13]= [1.5,4.5,0] ; p

    ; composition
		param_bnd[*,14]= [16.5, 23.5,0] ; positive scale factors  - MASS
		param_bnd[*,15]= [0, 1.0,0] ; olivine/pyroxene composition 
		param_bnd[*,16]= [0, 1.0,0] ; crystalline fraction
		param_bnd[*,17]= [0, 1.0,0] ; forsterite/enstatite composition 

    mcmc_cl2,par_bound=param_bnd,parameter=params

  end

  ; *************************************************** ;
  ; TWO COMPONENT GRAIN MODEL - Fixed MIPS
  ; params = [T1, a1, scale1, fo1, fc1, ff1,
  ;           T2, a2, scale2, fo2, fc2, ff2]
  6: begin

    print, 'Fitting two grain model...'

	  param_bnd = dblarr(3,12)
 		;param_bnd[*,0]= [30.0, 300.0,0] ; Temp limited to reasonable values
		param_bnd[*,0]= [0.1, 500.0,0] ; New distance bounds
		param_bnd[*,1]= [amin, 100.0,0] ; blowout radius as lower bound
		param_bnd[*,2]= [16.5, 23.5,0] ; positive scale factors
		param_bnd[*,3]= [0, 1.0,0] ; olivine/pyroxene composition 
		param_bnd[*,4]= [0, 1.0,0] ; crystalline fraction
		param_bnd[*,5]= [0, 1.0,0] ; forsterite/enstatite composition 
 		;param_bnd[*,6]= [100.0, 1000.0,0] ; Temp limited to reasonable values
		param_bnd[*,6]= [0.1, 500.0,0] ; New distance bounds
		param_bnd[*,7]= [amin, 100.0,0] ; blowout radius as lower bound
		param_bnd[*,8]= [16.5, 23.5,0] ; positive scale factors
		param_bnd[*,9]= [0, 1.0,0] ; olivine/pyroxene composition 
		param_bnd[*,10]= [0, 1.0,0] ; crystalline fraction
		param_bnd[*,11]= [0, 1.0,0] ; forsterite/enstatite composition 

    if keyword_set(nocrys) then begin
		  param_bnd[*,4]= [0, 1.0,1]
      param_bnd[*,5]= [0, 1.0,1]
      param_bnd[*,10]= [0, 1.0,1]
      param_bnd[*,11]= [0, 1.0,1]
      initparams[4] = 0.0
	    initparams[10] = 0.0
    endif

    ;scale_val=[1.010,.70, .0220, .0150, .0150, .0150,  1.50, .530,  .0120, .0150, .0150, .0150]

     mcmc_general,par_bound=param_bnd,parameter=params

  end

  ; *************************************************** ;
  ; SINGLE CONTINUOUS DISK - Fixed MIPS
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
	7: begin
 
	  param_bnd = dblarr(3,10)
		param_bnd[*,0]= [1, 5.0,0] ; rin
		param_bnd[*,1]= [0.0,10.0,0] ; r range
		param_bnd[*,2]= [-5,5 ,0] ; no limits on alpha
  	param_bnd[*,3]= [amin, 30.0,0] ; amin - blowout radius as lower bound
		param_bnd[*,4]= [0.047, 7.0,0] ; a range
		param_bnd[*,5]= [1.5, 5.5,0] ; no limits on pindex
		param_bnd[*,6]= [16.5, 23.5,0] ; no limits on mass
		param_bnd[*,7]= [0.0, 1.0,0] ; amorph composition
		param_bnd[*,8]= [0.0, 1.0,0] ; crystalline fraction

	  if keyword_set(nocrys) then begin
			param_bnd[*,8]= [0.0, 1.0,1]
      initparams[8] = 0.0
  	endif
		
		; crystalline composition
		param_bnd[*,9]= [0.0, 1.0,0]

    ;scale_val=[[.030,.030,.0220,0.80,.030,.0220, .0220, .0150, .015, .0150],$
    ;[.050,.050,.0350,1.60,.050,.0350, .0420, .0350, .0325, .0150],$
    ;[.080,.080,.0350,1.80,.080,.0350, .0620, .0350, .0325, .0150]]

    mcmc_general,par_bound=param_bnd,parameter=params
    ;mcmc__define_sd,'modelsinglespectrum',par_bound=param_bnd,parameter=params,scale_val=scale_val[*,params.scale_val]

  end
endcase
endelse

; End Program
return
end