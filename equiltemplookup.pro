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

; RESTORE, 'qtables_withcrys.sav'

COMMON GRAINPROPS

alllambda = [qastrosil.lambda[where(qastrosil.lambda lt qolivine.lambda[0])], $
          qolivine.lambda]
nl = n_elements(alllambda)


; read in optical data
qlookup, grainrad, alllambda, folivine, fcrystalline, fforst, $
         qabsall

; look up stellar fluxes
starfluxarr = stellar_emit[*,*,0] * folivine * (1.0-fcrystalline) + $
              stellar_emit[*,*,1] * (1.0-folivine) * (1.0-fcrystalline) + $
              stellar_emit[*,*,2] * fforst * fcrystalline + $
              stellar_emit[*,*,3] * (1.0-fforst) * fcrystalline


ii = interpol(dindgen(n_elements(agrain)), $
             alog(agrain), alog(grainrad))
jj = interpol(dindgen(n_elements(effectiveTempArray)), $
              effectiveTempArray, teff)
absflux = exp(interpolate(alog(starfluxarr), ii, jj, /grid))

; divide by distance to get actual flux
totabs = rebin(absflux,n_elements(grainrad),n_elements(distance)) $
         /transpose(rebin([distance]^2, $
                          n_elements(distance), $
                          n_elements(grainrad) ))

; print, 'total absorbed:', totabs

thisemit = olivine_emit*folivine*(1.0-fcrystalline) + $
           pyroxene_emit*(1.0-folivine)*(1.0-fcrystalline) + $
           forsterite_emit*fforst*fcrystalline + $
           enstatite_emit*(1.0-fforst)*fcrystalline

; look up temperature in table
tindex = dindgen(n_elements(tgrain))
aindex = interpol(dindgen(n_elements(agrain)), alog(agrain),alog(grainrad))
fluxtable = interpolate(thisemit, aindex, tindex, /grid)

tdust = dblarr(n_elements(grainrad), n_elements(distance))
for i=0,n_elements(grainrad)-1 do begin
    tdust[i,*] = exp(interpol(alog(tgrain), reform(alog(fluxtable[i,*])) $
                            , alog(totabs[i,*])))
endfor

;print, "tdust: ", tdust

return
end


