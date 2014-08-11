; +
; NAME
;  modeltwograin
;
; PURPOSE
;  Models the spectrum of a disk by passing parameters to modelonegrain twice
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
;  Written by TM (June 2013) as modeltwograin.pro
;  Organized and commented by EC (6/27/2014)
; -
; *************************************************** ;

function modeltwograin_old, lambda, params, spect1=spect1, spect2=spect2, mie=mie
; 2-population single-component grain model                                                 
; params = [T1, a1, fo1, fc1, scale1,                                                        
;           T2, a2, fo2, fc2, scale2]                                                        

; Create global variables relating to silicate features
COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs
COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, effectiveTempArray, stellar_emit

; Call modelone grain for each set of parameters
spect1 = modelonegrain_old(lambda, params[0:5], mie=mie)
spect2 = modelonegrain_old(lambda, params[6:11], mie=mie)

; Combine the two spectra 
flux = spect1+spect2

; Return the flux and exit
return, flux
end