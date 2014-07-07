; +
; NAME
;  qlookup
;
; PURPOSE
;  Recover/calculate silicate features given basic grain properties 
;
; INPUTS
;   GRAINRAD: Radius of grains
;   WAVELENGTH: Array of wavelengths to be modeled
;   OLIVRATIO: Ratio of olivine, w.r.t. number density. 
;              Assume the rest is pyroxene.
;   CRYSFRAC:
;   FORSTFRAC:
;
; KEYWORDS
;   NONE
;
; OUTPUTS
;   QABS: Value returned for actual spectrum modeling
;   Q*** = tables of Q values
;   QEXT: 
;   QSCAT: 
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
;  Written by TM (June 2013) as qlookup.pro
;  Organized and commented by EC (6/27/2014)
; -
; *************************************************** ;

pro qlookup, grainrad, wavelength, olivratio, crysfrac, forstfrac, $
             qabs, qext=qext, qscat=qscat

; Carry global variables relating to silicate features
COMMON grainprops, qastrosil, qolivine, qpyroxene, qenstatite, qforsterite, crystallineabs
;COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs                                           

qextall = dblarr(n_elements(grainrad), n_elements(wavelength),4)
qscatall = dblarr(n_elements(grainrad), n_elements(wavelength),4)

; agrain indices
ii = interpol(dindgen(n_elements(qastrosil.agrain)), $
              alog(qastrosil.agrain), alog(grainrad) ) $
              < n_elements(qastrosil.agrain)-1.0 


shortl = where(wavelength lt 0.2, nshort)
if nshort gt 0 then begin
; short wavelengths: use astronomical silcates 
  jj = interpol(dindgen(n_elements(qastrosil.lambda)), $
                alog(qastrosil.lambda), alog(wavelength[shortl]) )
  
  qextall[*, shortl,*] = $
     rebin(exp(interpolate(alog(Qastrosil.Qext), ii, jj, /grid)), $
     n_elements(grainrad), nshort, 4)
           
  qscatall[*, shortl,*] = $
     rebin(exp(interpolate(alog(Qastrosil.Qscat), ii, jj, /grid)), $
     n_elements(grainrad), nshort, 4)
      
endif 


; *************************************************** ;
; Compute for each graintype [olivine,pyroxene,forsterite,enstatite]
; *************************************************** ;


; *************************************************** ;
graintype =0
qstruct = qolivine
;longl = where(wavelength ge min(qstruct.agrain), nlong)
;if nlong gt 0 then begin
  jj = interpol(dindgen(n_elements(qstruct.lambda)), alog(qstruct.lambda), alog(wavelength) )
  qextall[*,*,graintype] = exp(interpolate((Qstruct.logQext), ii, jj, /grid))
  qscatall[*,*,graintype] = exp(interpolate((Qstruct.logQscat), ii, jj, /grid))
;endif


; *************************************************** ;
graintype =1
qstruct = qpyroxene
;longl = where(wavelength ge min(qstruct.agrain), nlong)
;if nlong gt 0 then begin
  jj = interpol(dindgen(n_elements(qstruct.lambda)), alog(qstruct.lambda), alog(wavelength) )
  qextall[*,*,graintype] = exp(interpolate((Qstruct.logQext), ii, jj, /grid))
  qscatall[*,*,graintype] = exp(interpolate((Qstruct.logQscat), ii, jj, /grid))
;endif


; *************************************************** ;
graintype =2
qstruct = qforsterite
;longl = where(wavelength ge min(qstruct.agrain), nlong)
;if nlong gt 0 then begin
  jj = interpol(dindgen(n_elements(qstruct.lambda2)), alog(qstruct.lambda2), alog(wavelength) )
  qextall[*,*,graintype] = exp(interpolate((Qstruct.logQext2), ii, jj, /grid))
  qscatall[*,*,graintype] = exp(interpolate((Qstruct.logQscat2), ii, jj, /grid))
;endif


; *************************************************** ;
graintype =3
qstruct = qenstatite
;longl = where(wavelength ge min(qstruct.agrain), nlong)
;if nlong gt 0 then begin
  jj = interpol(dindgen(n_elements(qstruct.lambda)), alog(qstruct.lambda), alog(wavelength) )
  qextall[*,*,graintype] = exp(interpolate((Qstruct.logQext), ii, jj, /grid))
  qscatall[*,*,graintype] = exp(interpolate((Qstruct.logQscat), ii, jj, /grid))
;endif



; *************************************************** ;
; MIX
; olivratio = olivine/(pyroxene+olivine)
; crysfrac = crystalline/(crystalline+amorphous)
; forstfrac = forsterite/(enstatite+forsterite)

qext = reform( qextall[*,*,0]*olivratio*(1.0-crysfrac) + $
               qextall[*,*,1]*(1.0-olivratio)*(1.0-crysfrac) + $
               qextall[*,*,2]*forstfrac*crysfrac + $
               qextall[*,*,3]*(1.0-forstfrac)*crysfrac, $
               n_elements(grainrad), n_elements(wavelength) )

qscat = reform( qscatall[*,*,0]*olivratio*(1.0-crysfrac) + $
                qscatall[*,*,1]*(1.0-olivratio)*(1.0-crysfrac) + $
                qscatall[*,*,2]*forstfrac*crysfrac + $
                qscatall[*,*,3]*(1.0-forstfrac)*crysfrac, $
                n_elements(grainrad), n_elements(wavelength) )

qabs = qext-qscat

stop
return
end
