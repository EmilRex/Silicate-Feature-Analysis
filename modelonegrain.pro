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

function modelonegrain, lambda, params, mie=mie
; 2-population single-component grain model                                                                                                          

; Carry global variables relating to silicate features
COMMON stellarprops, temptable, folivine, effectiveTemp, lambdastar, fluxstar
COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs
COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, $ 
       enstatite_emit, effectiveTempArray, stellar_emit
;restore, 'graintempdata.sav'
;restore, 'qtables_withcrys2.sav' ; qastrosil, qolivine, qpyroxene

; Define constants
rho_s = 3.3
AU_in_cm = 1.496e13
c = 2.998d10
h = 6.626d-27
k = 1.381d-16

; Define paramters passed                                                                                                                                                                                                 
temp     = params[0] ; set to dist
agr      = params[1]
dustmass = params[2]
foliv    = params[3]
fcrys    = params[4]
ffors    = params[5]


; Lookup qabs (not it will return a [1,n(lambda),4] array)
;   use keyword /separate to return qabs for each graintype individually
qlookup, [agr], lambda, foliv, fcrys, ffors, qabs, /separate

; Define scales
scale = dustmass*0.75/(rho_s*agr*1e-4)/AU_in_cm^2


FOR i = 0, 3 DO BEGIN
  
  ; Calculate Temperature
  temp = f(dist)
  
  ; Compute spectrum
  flux[*,i] = !pi*blackbody(lambda,temp)*(reform(qabs[*,*,i],n_elements(lambda))*scale)
  
ENDFOR

; Sum the fluxs over graintypes
tot_flux = dblarr(n_elements(lambda))

FOR j=0,n_elements(lambda) DO BEGIN
  tot_flux[i] = total(flux[i,*]) 
ENDFOR

; Return the flux and exit
return, tot_flux
end
