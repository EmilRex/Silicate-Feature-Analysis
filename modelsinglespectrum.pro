function modelsinglespectrum, lambda, params
; params = [rin, log(rout/rin), rlaw, amin, log(amax/amin), alaw, $
;           diskmass, foliv]
COMMON stellarprops, temptable, folivine, effectiveTemp, lambdastar, fluxstar
;COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs


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

; T = systime(1)

;print, 'param : ',params
diskspectrum, rin, rout, amin, amax, effectiveTemp, diskmass, $
              foliv, fcrys, fforst, $
              rlaw=rlaw, alaw=alaw, $
              lambda, flux

;PRINT, SYSTIME(1) - T,   ' Seconds - A1 '                                                                                                            

return, flux
end
