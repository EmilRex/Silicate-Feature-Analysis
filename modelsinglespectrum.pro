; Modified by EC on 7/31/14

function modelsinglespectrum, lambda, params, single=single


COMMON stellarprops, temptable, folivine, effectiveTemp, lambdastar, fluxstar
;COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs

if keyword_set(single) then begin
  
  ;params = [dist, agr, diskmass, foliv, fcrys, fforst]
  rin  = 10^params[0]
  rout = 0.0
  rlaw = 1.0
  amin = params[1]
  amax = 0.0
  alaw = 1.0
  diskmass = params[2]
  foliv = params[3]
  fcrys = params[4]
  fforst = params[5]
  
  diskspectrum, rin, rout, amin, amax, effectiveTemp, diskmass, $
    foliv, fcrys, fforst, rlaw=rlaw, alaw=alaw, lambda, flux, /single
  
endif else begin
  
  ;params = [rin, log(rout/rin), rlaw, amin, log(amax/amin), alaw, diskmass, foliv, fcrys, fforst]
  rin  = params[0]
  rout = exp(params[1])*rin
  rlaw = params[2]
  amin = params[3]
  amax = exp(params[4])*amin
  alaw = params[5]
  diskmass = params[6]
  foliv = params[7]
  fcrys = params[8]
  fforst = params[9]
  
  diskspectrum, rin, rout, amin, amax, effectiveTemp, diskmass, $
    foliv, fcrys, fforst, rlaw=rlaw, alaw=alaw, lambda, flux
  
endelse                                                                                               

return, flux
end
