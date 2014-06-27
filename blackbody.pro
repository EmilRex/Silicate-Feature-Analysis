; +
; NAME
;  blackbody
;
; PURPOSE
; Given temperature (T), return blackbody spectrum in units of 
; Jansky/asec^2 vs. wavelength in microns   
;
; INPUTS
;   WAVELENGTH: Wavelengths at which to calculate fluxes
;   T = Temperature of blackbody
;
; KEYWORDS
;   GRID: 
;
; OUTPUTS
;   B_NU: Array of computed fluxes for given wavelengths
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
;  Written by TM (June 2013) as blackbody.pro
;  Organized and commented by EC (6/27/2014)
; -
; *************************************************** ;

function blackbody, wavelength, T, grid=grid
                                                                                                         
; Define physical constants
c = 2.998d10
h = 6.626d-27
k = 1.381d-16

; -------------------------------------------------- ;
;if keyword_set(grid) then begin;
;  T2 =reform(T,n_elements(T))
;  tst1=rebin(wavelength,n_elements(T2),n_elements(wavelength))
;  tst2=rebin(T2,n_elements(T2),n_elements(wavelength))

;  hnu = rebin(h*nu, n_elements(wavelength), n_elements(T))
;  kT  = transpose(rebin(k*T, n_elements(T), n_elements(wavelength)))
;  B_nu = 2.0*hnu/(exp(hnu/kT)-1.0)/(rebin(wavelength*1e-4, n_elements(wavelength), n_elements(T)))^2

;  nu = c/(tst1*1.d-4)
;  hnukT = (h*nu/k)/tst2
;  p1 = 2.0*h*nu/(tst1*1e-4)^2
;  B_nu = p1/(exp(hnukT)-1.0)
;endif else begin
; -------------------------------------------------- ;

  ; Calculate blackbody
  nu = c/(wavelength*1.d-4)
  hnukT = (h*nu/k)/T
  p1 = 2.0*h*nu/(wavelength*1e-4)^2
  B_nu = p1/(exp(hnukT)-1.0)

;endelse

; Convert from ergs/s/cm^2/Hz/ster to janskys/arcsec^2                                                                                                                  
B_nu = B_nu*1.d23/206265.0^2

return, B_nu
end