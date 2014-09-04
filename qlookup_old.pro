pro qlookup_old, grainrad, wavelength, olivratio, crysfrac, forstfrac, $
             qabs, qext=qext, qscat=qscat, separate=separate
; grainrad = a
; wavelength = lambda
; olivratio = ratio of olivine, w.r.t. number density . 
;             Assume the rest is pyroxene.
; Q*** = tables of Q values
; qext, qscat = return values of Q
COMMON grainprops, qastrosil, qolivine, qpyroxene, qenstatite, qforsterite, $
   crystallineabs

restore, 'qtables_withcrys2.sav'

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


;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;

graintype =0
         qstruct = qolivine
;   longl = where(wavelength ge min(qstruct.agrain), nlong)
;   t = systime(1)

;   if nlong gt 0 then begin
      jj = interpol(dindgen(n_elements(qstruct.lambda)), $
                    alog(qstruct.lambda), alog(wavelength) )
      qextall[*,*,graintype] = $
         exp(interpolate((Qstruct.logQext), ii, jj, /grid))
      qscatall[*,*,graintype] = $
         exp(interpolate((Qstruct.logQscat), ii, jj, /grid))
;   endif

;print,systime(1)-t,' sec1'
;t=systime(1)

   graintype =1
         qstruct = qpyroxene
;   longl = where(wavelength ge min(qstruct.agrain), nlong)
;   if nlong gt 0 then begin
      jj = interpol(dindgen(n_elements(qstruct.lambda)), $
                    alog(qstruct.lambda), alog(wavelength) )
      qextall[*,*,graintype] = $
         exp(interpolate((Qstruct.logQext), ii, jj, /grid))
      qscatall[*,*,graintype] = $
         exp(interpolate((Qstruct.logQscat), ii, jj, /grid))
;   endif


;print,systime(1)-t,' sec1'
;t=systime(1)

   graintype =2
         qstruct = qforsterite
;   longl = where(wavelength ge min(qstruct.agrain), nlong)
;   if nlong gt 0 then begin
      jj = interpol(dindgen(n_elements(qstruct.lambda2)), $
                    alog(qstruct.lambda2), alog(wavelength) )
      qextall[*,*,graintype] = $
         exp(interpolate((Qstruct.logQext2), ii, jj, /grid))
      qscatall[*,*,graintype] = $
         exp(interpolate((Qstruct.logQscat2), ii, jj, /grid))
;   endif


;stop

;print,systime(1)-t,' sec1'
;t=systime(1)
   graintype =3
         qstruct = qenstatite
;   longl = where(wavelength ge min(qstruct.agrain), nlong)
;   if nlong gt 0 then begin
      jj = interpol(dindgen(n_elements(qstruct.lambda)), $
                    alog(qstruct.lambda), alog(wavelength) )
      qextall[*,*,graintype] = $
         exp(interpolate((Qstruct.logQext), ii, jj, /grid))
      qscatall[*,*,graintype] = $
         exp(interpolate((Qstruct.logQscat), ii, jj, /grid))
;   endif

;print,systime(1)-t,' sec1'
;stop
;t=systime(1)


;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;

IF NOT KEYWORD_SET(separate) THEN BEGIN

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

ENDIF ELSE BEGIN

  qabs = dblarr(n_elements(grainrad),n_elements(wavelength),4)
  qext = dblarr(n_elements(grainrad),n_elements(wavelength),4)
  qscat = dblarr(n_elements(grainrad),n_elements(wavelength),4)

  qext[*,*,0] = qextall[*,*,0]*olivratio*(1.0-crysfrac)
  qext[*,*,1] = qextall[*,*,1]*(1.0-olivratio)*(1.0-crysfrac)
  qext[*,*,2] = qextall[*,*,2]*forstfrac*crysfrac
  qext[*,*,3] = qextall[*,*,3]*(1.0-forstfrac)*crysfrac
  
  qscat[*,*,0] = qscatall[*,*,0]*olivratio*(1.0-crysfrac)
  qscat[*,*,1] = qscatall[*,*,1]*(1.0-olivratio)*(1.0-crysfrac)
  qscat[*,*,2] = qscatall[*,*,2]*forstfrac*crysfrac
  qscat[*,*,3] = qscatall[*,*,3]*(1.0-forstfrac)*crysfrac
  
  
  FOR i=0,3 DO qabs[*,*,i] = qext[*,*,i] - qscat[*,*,i]
  
ENDELSE

;stop
return
end
