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
;restore, 'qtables_withcrys.sav' ; qastrosil, qolivine, qpyroxene                                                                                                                                                                   
 
; Define paramters passed = [T1, a1, f1, scale1]                                                                                                                                                                                                    
temp1     = params[0]
agr1      = params[1]
dustmass1 = params[2]
foliv1    = params[3]
fcrys1    = params[4]
ffors1    = params[5]

; Define constants
rho_s = 3.3
AU_in_cm = 1.496e13

; Define scales
scale1 = dustmass1*0.75/(rho_s*agr1*1e-4)/AU_in_cm^2
scale2 = 0.0

; Lookup relevant data
;if keyword_set(mie) then begin
  qlookup, [agr1], lambda, foliv1, fcrys1, ffors1, qabs1
  abscoeff = replicate(0,n_elements(lambda))
;endif else begin
;  qlookup, [agr1], lambda, foliv1, 0.0, 0.0, qabs1
;  abscoeff = crystallineopacities(lambda, temp1, ffors1, /grid)
;  stop
;  scale1 = scale1 * (1.0-fcrys1)
;  scale2 = dustmass1*fcrys1/AU_in_cm^2
;endelse

; Compute spectrum
spect1 = !pi*blackbody(lambda,temp1)* $
         ( reform(qabs1,n_elements(lambda))*scale1 + $
           reform(abscoeff,n_elements(lambda))*scale2 )

; Return the flux and exit
flux = spect1
return, flux
end
