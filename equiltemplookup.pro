pro equiltemplookup, teff, grainrad, distance, $
                     tdust, $
                     folivine, fcrystalline, fforst
  
;print, "Distance: ", distance                 
; 
; calculates equilibrium temperature of grains at the given grain size
; from table lookup
; INPUTS:
; teff = effective temperature of the star in K
; distance = distance from the star, in units of R_star. can be 1D array.
; grainrad = radius of grain, in microns.  can be 1D array.
; OUTPUT:
; tdust = the equilibrium temperature vs. grainrad, 
;         n_elements(grainrad) x n_elements(distance) array

; restore, file='graintempdata.sav'
; tgrain, agrain, olivine_emit, pyroxene_emit, 
;                 forsterite_emit, enstatite_emit
; tables of total emission versus temperature and grain size

COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, effectiveTempArray, stellar_emit
COMMON GRAINPROPS

alllambda = [qastrosil.lambda[where(qastrosil.lambda lt qolivine.lambda[0])], $
          qolivine.lambda]
nl = n_elements(alllambda)


; read in optical data
qlookup, grainrad, alllambda, folivine, fcrystalline, fforst, qabsall

; look up stellar fluxes
starfluxarr = dblarr(size(stellar_emit,/dimensions))

starfluxarr[*,*,0] = stellar_emit[*,*,0] * folivine * (1.0-fcrystalline)
starfluxarr[*,*,1] = stellar_emit[*,*,1] * (1.0-folivine) * (1.0-fcrystalline)
starfluxarr[*,*,2] = stellar_emit[*,*,2] * fforst * fcrystalline
starfluxarr[*,*,3] = stellar_emit[*,*,3] * (1.0-fforst) * fcrystalline


ii = interpol(dindgen(n_elements(agrain)), $
             alog(agrain), alog(grainrad))
jj = interpol(dindgen(n_elements(effectiveTempArray)), $
              effectiveTempArray, teff)

absflux = dblarr(n_elements(ii),4)
totabs = dblarr(n_elements(ii),n_elements(distance),4)

for g=0,3 do begin
  absflux[*,g] = exp(interpolate(alog(starfluxarr[*,*,g]), ii, jj, /grid)) 
  
  ; divide by distance to get actual flux
  totabs[*,*,g] = rebin(absflux[*,g],n_elements(grainrad),n_elements(distance)) $
           /transpose(rebin([distance]^2, n_elements(distance), n_elements(grainrad) ))
endfor

nte = size(olivine_emit,/dimensions)
thisemit = dblarr(nte[0],nte[1],4)
thisemit[*,*,0] = olivine_emit*folivine*(1.0-fcrystalline)
thisemit[*,*,1] = pyroxene_emit*(1.0-folivine)*(1.0-fcrystalline)
thisemit[*,*,2] = forsterite_emit*fforst*fcrystalline
thisemit[*,*,3] = enstatite_emit*(1.0-fforst)*fcrystalline

; look up temperature in table
tindex = dindgen(n_elements(tgrain))
aindex = interpol(dindgen(n_elements(agrain)), alog(agrain),alog(grainrad))

fluxtable = dblarr(n_elements(aindex),n_elements(tindex),4)
tdust = dblarr(n_elements(grainrad), n_elements(distance),4)

for g=0,3 do begin
  
  fluxtable[*,*,g] = interpolate(thisemit[*,*,g], aindex, tindex, /grid)
  for i=0,n_elements(grainrad)-1 do begin
      tdust[i,*,g] = exp(interpol(alog(tgrain), reform(alog(fluxtable[i,*,g])) $
                     , alog(totabs[i,*,g])))
  endfor

endfor

return
end