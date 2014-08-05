; Modified by EC on 7/31/14

function modelsinglespectrum, lambda, params, $
  single=single, multi=multi, disk=disk, separate=separate


COMMON stellarprops, temptable, folivine, effectiveTemp, lambdastar, fluxstar
;COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs

; *************************************************** ;
; Single
if keyword_set(single) then begin
  
  ;params = [dist, agr, diskmass, foliv, fcrys, fforst]
  rin  = 10^params[0]
  rout = 0.0
  rlaw = 1.0
  amin = params[1]
  amax = 0.0
  alaw = 1.0
  diskmass = 10^params[2]
  foliv = params[3]
  fcrys = params[4]
  fforst = params[5]
  
  diskspectrum, rin, rout, amin, amax, effectiveTemp, diskmass, $
    foliv, fcrys, fforst, rlaw=rlaw, alaw=alaw, lambda, flux, /single
  
endif 

; *************************************************** ;
; Multi
if keyword_set(multi) then begin

  ;params = [dist, agr, diskmass, foliv, fcrys, fforst]
  rin1  = 10^params[0]
  amin1 = params[1]
  diskmass1 = 10^params[2]
  foliv1 = params[3]
  fcrys1 = params[4]
  fforst1 = params[5]
  
  rin2  = 10^params[6]
  amin2 = params[7]
  diskmass2 = 10^params[8]
  foliv2 = params[9]
  fcrys2 = params[10]
  fforst2 = params[11]
  
  rout = 0.0
  rlaw = 1.0
  amax = 0.0
  alaw = 1.0
  
  diskspectrum, rin1, rout, amin1, amax, effectiveTemp, diskmass1, $
    foliv1, fcrys1, fforst1, rlaw=rlaw, alaw=alaw, lambda, flux1, /single

  diskspectrum, rin2, rout, amin2, amax, effectiveTemp, diskmass2, $
    foliv2, fcrys2, fforst2, rlaw=rlaw, alaw=alaw, lambda, flux2, /single

  flux = dblarr(n_elements(lambda),8)
  flux[*,0:3] = flux1
  flux[*,4:7] = flux2

endif

; *************************************************** ;
; Disk
if keyword_set(disk) then begin
  ;params = [rin, log(rout/rin), rlaw, amin, log(amax/amin), alaw, diskmass, foliv, fcrys, fforst]
  rin  = 10^params[0]
  rout = exp(params[1])*rin
  rlaw = params[2]
  amin = params[3]
  amax = exp(params[4])*amin
  alaw = params[5]
  diskmass = 10^params[6]
  foliv = params[7]
  fcrys = params[8]
  fforst = params[9]
  
  diskspectrum, rin, rout, amin, amax, effectiveTemp, diskmass, $
    foliv, fcrys, fforst, rlaw=rlaw, alaw=alaw, lambda, flux
endif

; *************************************************** ;
; Sum spectrum unless instructed otherwise
if not keyword_set(separate) then begin
  flux = total(flux,2)
endif

return, flux
end
