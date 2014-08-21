; +
; NAME
;  diskspectrum
;
; PURPOSE
;  Computes the spectrum of a circumstellar dust disk 
;
; INPUTS
;   RIN: Disk inner radius (stellar radii)
;   ROUT: Disk outer radius (stellar radii)
;   AMIN: Minimum grain size (microns)
;   AMAX: Maximum grain size (microns)
;   TEFF: Effective temperature of star (Kelvin)
;   DISKMASS: Normalized mass of dust disk
;   FOLIVINE = olivine fraction - folivine*(1-fcrystalline)
;   FCRYSTALLINE: crystalline fraction
;   FFORST fraction - fforst*fcrystalline
;     enstatite fraction = (1-fforst)*fcrystalline
;     pyroxene fraction = (1-folivine)*(1-fcrystalline)
;   RLAW: Disk radius size power law
;   ALAW: Grain size power law
;   LAMBDA: Wavelengths for which to model the flux
;
; KEYWORDS
;   VERBOSE: Set to run additional diagnostics
;   SINGLE: Calculate for one grain size at one location
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
; Written by HJ-C ???
;  Adapted by TM (June 2013) as diskspectrum.pro
;  Organized and commented by EC (7/18/2014)
;  Added comprehensive computation time diagnostics
;   see test_disk_benchmarks.pro for analysis
;  Removed benchmarking scheme and commented (7/30/14)
;  Modified to also include single model by EC on (7/31/14)
; -
; *************************************************** ;
pro diskspectrum, rin, rout, amin, amax, Teff, diskmass, $ 
                  folivine, fcrystalline, fforst, fwaterice, $
                  rlaw=rlaw, alaw=alaw, verbose=verbose, $
                  lambda, flux, mie=mie, single=single

;COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, Qwaterice, crystallineabs

; Set constants
rho_s = 3.3 ;bulk density of grains in g/cm^3
AU_in_cm = 1.496e13
c  = 2.998d10
h = 6.626d-27
k = 1.381d-16

if not keyword_set(rlaw) then rlaw=1.5
if not keyword_set(alaw) then alaw=3.5

if keyword_set(single) then begin
  
  ; Set up a grain parameters
  NA = 1
  dla = 1.0
  aall = [amin]
  
  ; Set up radius parameters
  NR = 1
  dlr = 1.0
  rall = [rin]
  sigma = 1.0/(2.0*!dpi*rall^2)
  
endif else begin
  
  ; Set up a grain parameters
  dloga = 0.02d
  NA = ceil(alog10(amax/amin)/dloga)+1
  dloga = alog10(amax/amin)/(NA-1.0d)
  aall = amin*10.0D^(dindgen(NA)*dloga)
  dla = replicate(dloga*alog(10.0D), NA)
  dla[0] = 0.5*dla[0]
  dla[NA-1] = dla[0]
  
  if keyword_set(ring) then begin
    radius = rin
    width = rout
    rin = (radius - 4.0*width) > 5.0
    rout = radius + 4.0*width
  endif
  
  ; Set up radius parameters
  if (abs(rout/rin-1.0) le 1e-3) then begin ; single radius ring, no width
    NR = 1
    dlr = 1.0
    rall = [rin]
    sigma = 1.0/(2.0*!dpi*rall^2)
  endif else begin
    dlogr = 0.02d
    NR = ceil(alog10(rout/rin)/dlogr) + 1
    dlogr = alog10(rout/rin)/(NR-1.0d)
    rall = rin*10.0D^(dindgen(NR)*dlogr)
    dlr = replicate(dlogr*alog(10.0D), NR)
    dlr[0] = 0.5*dlr[0]
    dlr[NR-1] = dlr[0]
    if keyword_set(ring) then $
      sigma = exp(-0.5*(rall-radius)^2/width^2) $
    else $
      sigma = rall^(-rlaw)
  endelse
  
endelse


; Get temperatures - don't include crystalline silicates
equiltemplookup, Teff, aall, rall, tempr, folivine, fcrystalline, fforst, fwaterice

; Get Qabs, from Mie calculation
qlookup, aall, lambda, folivine, fcrystalline, fforst, fwaterice, qabsall, /separate

; *************************************************** ;
; Calculate spectrum
flux = dblarr(n_elements(lambda),5)
NL = n_elements(lambda)
if nl le 1 then print, 'OH WOE!'

rintegrand = 2.0*!dpi * rall * sigma*dlr*rall ; multiply by  R_sun^2
crosssection = !dpi*(aall*1e-4)^2 * aall^(-alaw)*dla*aall ; Multiply by 1e-4 for cm

; Integrate to get normalization 
totsigma = total(rintegrand)
totmass = 4.0/3.0 * rho_s * total(aall*1e-4 * crosssection) * totsigma
amorphnorm = diskmass/totmass
crysnorm = diskmass/totsigma/total(crosssection)

; *************************************************** ;
; Calculate flux for each grain type
for g=0,4 do begin
  ; Calculate Blackbody
  temp_t1 =reform(transpose(tempr[*,*,g]),NA*NR) 
  wave_new1= (REBIN(lambda,NA*NR,NL))*1.d-4
  temp_t2= (REBIN(temp_t1,NA*NR,NL))
  
  
  hnukT = ((h*c/k)/wave_new1)/temp_t2
  tmp_p1 = (2.0*h*c)/(wave_new1)^3
  brightness = (tmp_p1*(1.d23/206265.0^2))/(exp(hnukT)-1.0)
  subl = where(temp_t2 gt 1000.0) ; include dust sublimation
  
  if subl[0] ne -1 then brightness[subl] = 0.0
  ;  tmp_p1 = REBIN(rintegrand,1,n_elements(rintegrand)*n_elements(tempr[*,0]))
  ;  rintegrand2 = REBIN(tmp_p1,n_elements(tmp_p1),n_elements(lambda))
  ;  brightness2 = reform(brightness,NR,NA,NL)
  
  brightness = transpose(brightness)
  ;brightness2 = (brightness)*abscoeff
  tmp_p1 = FLTARR(NA,NL) 
  
  
  for i=0,NA-1 do begin
      intFdr = (matrix_multiply(brightness[*,0+NR*(i):NR-1+NR*(i)], rintegrand)*qabsall[i,*,g])*(crosssection[i]*amorphnorm)
      ;intFdr2 = matrix_multiply(brightness2[*,0+NR*(i):NR-1+NR*(i)], rintegrand)*(crosssection[i]*crysnorm)
      tmp_p1[i,*] = intFdr;*(1.0-fcrystalline)+intFdr2*fcrystalline
  endfor

  ; Return the flux and end
  flux[*,g] = total(tmp_p1,1)
  flux[*,g] = flux[*,g]/AU_in_cm^2

endfor

return
end
