; +
; NAME
;  modelonegrain
;
; PURPOSE
;  Models the spectrum of a dust disk 
;
; INPUTS
;   LAMBDA: Wavelengths for which to model the flux
;   PARAMS: array of disk properties used in model
;
; KEYWORDS
;   MIE: Set to use MIE theory
;
; OUTPUTS
;   FLUX: array of fluxes corresponding to given wavelengths for the modeled spectrum
;
; AUTHORS
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
; MODIFICATION HISTORY
;  Written by TM (June 2013) as modelonegrain.pro
;  Organized and commented by EC (6/27/2014)
;  Modified by EC (7/21/14)
;   Takes dist as given and computes temp for each graintype
;   instead of visa versa.
;   Incorporates elements of CC's dist.pro
; -
; *************************************************** ;

function modelonegrain, lambda, params, mie=mie, temp=temp
; 2-population single-component grain model                                                                                                          

; Carry global variables relating to silicate features
COMMON stellarprops, temptable, folivine, effectiveTemp, lambdastar, fluxstar
COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs
COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, $ 
       enstatite_emit, effectiveTempArray, stellar_emit
COMMON file_path, in_dir, out_dir, fit_name, object_name

COMMON star_params, t_star, dist_to_star, star_lambda, star_fnu

; Test
;object_name = 'HD109573'

;restore, 'graintempdata.sav'
;restore, 'qtables_withcrys2.sav' ; qastrosil, qolivine, qpyroxene

; Define constants
rho_s = 3.3
pc_in_AU = 2.0626e5
AU_in_cm = 1.496e13
mm_to_cm = 1.0e-4
c = 2.998d10 ; cm/s
h = 6.626d-27 ; erg s
k = 1.381d-16 ; erg/K

; Define paramters passed                                                                                                                                                                                                 
dist     = params[0] ; remember to update bounds
agr      = params[1]
dustmass = params[2]
foliv    = params[3]
fcrys    = params[4]
ffors    = params[5]

; *************************************************** ;
; Look up relevant data
; 
; Lookup qabs (not it will return a [1,n(lambda),4] array)
;   use keyword /separate to return qabs for each graintype individually
;   Order: olivine, pyroxene, forsterite, enstatite
qlookup, [agr], lambda, foliv, fcrys, ffors, qabs, /separate

; *************************************************** ;
; Begin calculation

; Define arrays
flux = dblarr(n_elements(lambda),4)
tot_flux = dblarr(n_elements(lambda))
temp = dblarr(4)

; Define scale
scale = dustmass*0.75/(rho_s*agr*1e-4)/AU_in_cm^2

; Convert distance in parsecs to au
dist_val_AU = dist_to_star*pc_in_AU

; Get q_abs over phot lambdas 
qlookup, [agr], star_lambda, foliv, fcrys, ffors, qabs_phot, /separate

; Calculate LHS integration constant
int_const = (c*(dist_val_AU^2))/(((star_lambda*mm_to_cm)^2)*(dist^2)) ; mm->mum

FOR i = 0, 3 DO BEGIN
  
  ; ******************************** ;
  ; Calculate LHS: Heating from star
  
  ; Put photosphere model in terms of data lambda
  ;phot_spec = 10^interpol(alog10(star_fnu*1.0e-23),alog10(star_lambda*mm_to_cm),alog10(uniq_lambda*mm_to_cm))

  lhs = INT_TABULATED(star_lambda*mm_to_cm,qabs_phot[*,*,i]*int_const*star_fnu*1.0e-23)
  
  ; ******************************** ;
  ; Calculate RHS: Cooling from dust
  
  ; Define temperature range
  t = 3.0*(findgen(1000) + 1.0)
  
  blambda = dblarr(n_elements(t),n_elements(star_lambda))
  rhs = dblarr(n_elements(t))
  
  ; Iterate for each temperature
  FOR j=0,(n_elements(t)-1) DO BEGIN
    blambda[j,*] = ( (2.0*h*(c^2))/((star_lambda*mm_to_cm)^5) )/( exp( (h*c)/((star_lambda*mm_to_cm)*k*t[j]) ) -1.0 )
    rhs[j] = 4.0*INT_TABULATED((star_lambda*mm_to_cm),qabs_phot[*,*,i]*blambda[j,*])
  ENDFOR
  
  ; Calculate Temperature
  temp[i] = 10^interpol(alog10(t),alog10(rhs),alog10(lhs))

  ; Compute spectrum
  ; Try multiplying "scale" by the "grain fraction" - will make sure each graintype only contributes
  ; its fraction of the flux. Although, the fractions are already implemented in qabs...
  ; qabs is really qabs*grain_mass_fraction.
  flux[*,i] = !pi*blackbody(lambda,temp[i])*(reform(qabs[*,*,i],n_elements(lambda))*scale)
  
ENDFOR


;print, temp
; Sum the fluxs over graintypes
FOR j=0,n_elements(lambda)-1 DO BEGIN
  tot_flux[j] = total(flux[j,*]) 
ENDFOR

; test case HR4796, HD109573
;plot,final_wave,final_spec
;oplot,lambda,flux,psym=4

;stop
return, tot_flux
end
