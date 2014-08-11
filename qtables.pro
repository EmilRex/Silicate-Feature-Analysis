pro qtables, Qastrosil, Qolivine, Qpyroxene
; read in files of optical constants
; generate Qext, Qabs values, as a function of grain size & wavelength
; save in qtables.sav
; three structures: Qolivine, Qpyroxine, Qastrosil
; each structure contains:
; - agrain: 0.1 to 1mm, 100 bins per decade
; - lambda: read in from files of refraction indices
; - Qext: 2D array, value for each agrain & lambda, index in that order.
; - Qscat: " "

; from 0.1 microns to 1 mm
agrain = 0.1D * 10.0D^(dindgen(401)*0.01d) 
NA = n_elements(agrain)
; maximum grain size parameter allowed to mie_single is 12000.0
maxx = 12000.0


; DRAINE & LEE astronomical silicates

; read in optical data
readcol, 'opacities/eps_Sil', lasil, e1, e2, nasil, kasil, $
         format='(F,F,F,F,F)'
nasil = nasil + 1.0
ii = where(lasil le 0.2)
lasil = lasil[reverse(ii)]
nasil = nasil[reverse(ii)]
kasil = kasil[reverse(ii)]

; when calling mie_single:
; 1st argument: x =2*pi*a/lambda - only one lambda at a time!
; 2nd argument: complex refractive index, expressed as an IDL dcomplex number
; 3rd argument: return vector for Qext
; 4th argument: return vector for Qscat

; calculate Qastrosil
nl = n_elements(lasil)
Qastrosil = {agrain:agrain, lambda:lasil, $
             Qext:dblarr(NA,nl), Qscat:dblarr(NA,nl)}

print, 'Astronomical silicates'
for i=0, nl-1 do begin
   mie_single, maxx, complex(nasil[i], -kasil[i]), qext_lg, qscat_lg
   Qastrosil.Qext [*,i] = qext_lg
   Qastrosil.Qscat[*,i] = qscat_lg

   x = 2.0*!dpi*agrain/lasil[i]
   smallgrains = where(x le maxx)
   mie_single, x[smallgrains], complex(nasil[i], -kasil[i]), qext, qscat
   Qastrosil.Qext [smallgrains,i] = qext
   Qastrosil.Qscat[smallgrains,i] = qscat
   print, 'lambda = ', lasil[i]
endfor


; OLIVINE

; read in optical data
readcol, 'opacities/olmg50.lnk', loliv, noliv, koliv, format='(F,F,F)'

; calculate Qolivine
nl = n_elements(loliv)
Qolivine = {agrain:agrain, lambda:loliv, $
             Qext:dblarr(NA,nl), Qscat:dblarr(NA,nl)}

print, 'Olivine'
for i=0, nl-1 do begin
   mie_single, maxx, complex(noliv[i], -koliv[i]), qext_lg, qscat_lg
   Qolivine.Qext [*,i] = qext_lg
   Qolivine.Qscat[*,i] = qscat_lg

   x = 2.0*!dpi*agrain/loliv[i]
   smallgrains = where(x le maxx)
   mie_single, x[smallgrains], complex(noliv[i], -koliv[i]), qext, qscat
   Qolivine.Qext [smallgrains,i] = qext
   Qolivine.Qscat[smallgrains,i] = qscat
   print, 'lambda = ', loliv[i]
endfor


; PYROXENE

; read in optical data
readcol, 'opacities/pyrmg50.lnk', lpyr, npyr, kpyr, format='(F,F,F)'

; calculate Qpyroxene
nl = n_elements(lpyr)
Qpyroxene = {agrain:agrain, lambda:lpyr, $
             Qext:dblarr(NA,nl), Qscat:dblarr(NA,nl)}

print, 'Pyroxene'
for i=0, nl-1 do begin
   mie_single, maxx, complex(npyr[i], -kpyr[i]), qext_lg, qscat_lg
   Qpyroxene.Qext [*,i] = qext_lg
   Qpyroxene.Qscat[*,i] = qscat_lg

   x = 2.0*!dpi*agrain/lpyr[i]
   smallgrains = where(x le maxx)
   mie_single, x[smallgrains], complex(npyr[i], -kpyr[i]), qext, qscat
   Qpyroxene.Qext [smallgrains,i] = qext
   Qpyroxene.Qscat[smallgrains,i] = qscat
   print, 'lambda = ', lpyr[i]
endfor

save, file='qtables.sav', Qastrosil, Qolivine, Qpyroxene

return
end
