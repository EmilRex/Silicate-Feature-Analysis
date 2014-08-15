pro tempdata
; read in optical data
COMMON grainprops, qastrosil, qolivine, qpyroxene, qenstatite, qforsterite, qwaterice, crystallineabs
RESTORE, 'qtables_withcrys2.sav'
RESTORE, 'qwaterice.sav'


NT = 751
tgrain = 2.0*10.0d^(3.0*dindgen(NT)/(NT-1.0))

agrain = qastrosil.agrain
NA = n_elements(agrain)

; from 2 K to 2000K
olivine_emit = dblarr(NA,NT)
pyroxene_emit = dblarr(NA,NT)
forsterite_emit = dblarr(NA,NT)
enstatite_emit = dblarr(NA,NT)
waterice_emit = dblarr(NA,NT)

alllambda = [qastrosil.lambda[where(qastrosil.lambda lt qolivine.lambda[0])], $
          qolivine.lambda]
nl = n_elements(alllambda)
NA = n_elements(agrain)

qlookup, agrain, alllambda, 1.0, 0.0, 0.0, 0.0, qabsolivine
qlookup, agrain, alllambda, 0.0, 0.0, 0.0, 0.0, qabspyroxene
qlookup, agrain, alllambda, 0.0, 1.0, 1.0, 0.0, qabsforsterite
qlookup, agrain, alllambda, 0.0, 1.0, 0.0, 0.0, qabsenstatite
qlookup, agrain, alllambda, 0.0, 0.0, 0.0, 1.0, qabswaterice

nu = 2.998e10/(alllambda*1e-4)
dnu = -[nu[1]-nu[0], 0.5*(nu[2:nl-1]-nu[0:nl-3]), nu[nl-1] - nu[nl-2]]

for i=0, NA-1 do begin
    print, 'grain size ', agrain[i]
    for j=0, NT-1 do begin
        olivine_emit[i,j] = total(qabsolivine[i,*]*dnu*blackbody(alllambda,tgrain[j])*206265.0^2*1e-23)
        pyroxene_emit[i,j] = total(qabspyroxene[i,*]*dnu*blackbody(alllambda,tgrain[j])*206265.0^2*1e-23)
        forsterite_emit[i,j] = total(qabsforsterite[i,*]*dnu*blackbody(alllambda,tgrain[j])*206265.0^2*1e-23)
        enstatite_emit[i,j] = total(qabsenstatite[i,*]*dnu*blackbody(alllambda,tgrain[j])*206265.0^2*1e-23)
        waterice_emit[i,j] = total(qabswaterice[i,*]*dnu*blackbody(alllambda,tgrain[j])*206265.0^2*1e-23)
    endfor
endfor

; calculate array of stellar irradiation versus effective temperature
; divide by distance^2 in units of stellar radius to get incident flux
; list all kurucz models 
cmd = 'ls ../kurucz/kp00*g40.txt'
spawn, cmd, photfiles
; read in temperatures 
strbeg = strpos(photfiles, '_')
strend = strpos(photfiles, 'g')
effectiveTempArray = fltarr(n_elements(photfiles))
for i=0,n_elements(photfiles)-1 do $
   effectiveTempArray[i] = float(strmid(photfiles[i], strbeg[i]+1, strend[i]-strbeg[i]-1))
ii = sort(effectiveTempArray)
photfiles = photfiles[ii]
effectiveTempArray = effectiveTempArray[ii]

print, 'calculating incident stellar fluxes'
nu = 2.998e10/(alllambda*1e-4)
dnu = -[nu[1]-nu[0], 0.5*(nu[2:nl-1]-nu[0:nl-3]), nu[nl-1] - nu[nl-2]]
stellar_emit = dblarr(n_elements(agrain), n_elements(photfiles), 5)

for i=0, n_elements(photfiles)-1 do begin
   readcol, photfiles[i], modlambda, modflux, format = '(F,F)', /silent
; convert flux units to ergs/cm^2/s/Hz/ster
   flux = modflux*modlambda*modlambda*3.336e-19/!dpi
; convert from angstroms to microns
   lambda = modlambda*1e-4
   fstar = interpol(flux, lambda, alllambda)

   for j=0, 4 do begin 
; each grain type, integrate over stellar spectrum
      case j of 
         0: qabs = qabsolivine
         1: qabs = qabspyroxene
         2: qabs = qabsforsterite
         3: qabs = qabsenstatite
         4: qabs = qabswaterice
      endcase

      absflux = matrix_multiply(qabs, fstar*dnu)
      stellar_emit[*,i,j] = absflux*0.25
      
   endfor
endfor

save, file='graintempdata.sav', tgrain, agrain, $
      olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, waterice_emit, $
      effectiveTempArray, stellar_emit
return
end


