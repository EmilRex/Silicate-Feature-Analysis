function modeltwograin_old, lambda, params, spect1=spect1, spect2=spect2, mie=mie
; 2-population single-component grain model                                                 
; params = [T1, a1, fo1, fc1, scale1,                                                        
;           T2, a2, fo2, fc2, scale2]                                                        

COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs
COMMON GRAINTEMPDATA, tgrain, agrain, olivine_emit, pyroxene_emit, forsterite_emit, enstatite_emit, effectiveTempArray, stellar_emit

spect1 = modelonegrain_old(lambda, params[0:5], mie=mie)
spect2 = modelonegrain_old(lambda, params[6:11], mie=mie)


flux = spect1+spect2

return, flux
end
