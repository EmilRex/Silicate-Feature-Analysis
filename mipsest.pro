pro mipsest

;;read in stellar properties from mark pecaut's initial spreadsheet
;readcol,'younga_stellarprop.txt',hip,av,eav,hplx,ehplx,cplx,ecplx,dist,$
;  pdh,mdh,dc,pdc,mdh,dgrp,teff,logteff,loglstar,eloglstar,clogl,eclogl,mdlogl
;read in stellar properties from mark pecaut's final spreadsheet
fmt='d,a,a,d,d,i,d,d,d'
readcol,'younga_stellarprop2.txt',F=fmt,hip,spt,grp,av,eav,teff,logteff,$
  loglstar,eloglstar,dist

lstar = 10^loglstar*3.9d33
dist = dist*3.08d18
teffstr = strtrim(string(long(teff)),2)

openw,1,'younga_mipsest.txt'
for i = 0, n_elements(hip)-1 do begin

;scale kurucz model atmosphere model to luminosity and distance
    readcol,'kp00_'+teffstr[i]+'g40.txt',lambda,flambda
    fnu = lambda^2*flambda/3.0d18
    norm = lstar[i]/4./3.14159/dist[i]^2/int_tabulated(lambda,flambda)
    fnu = norm*fnu

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;estimate photospheric fluxes integrated over MIPS bands (mJy)
    readcol,'2massj.txt',wavej,resj
    readcol,'2massh.txt',waveh,resh
    readcol,'2massk.txt',wavek,resk
    readcol,'MIPSfiltsumm24.csv',wave24,res24
    readcol,'MIPSfiltsumm70.csv',wave70,res70
    readcol,'MIPSfiltsumm160.csv',wave160,res160
    intflux=10^(interpol(alog10(fnu),alog10(lambda*1.0e-4),alog10(wavej)))
    fluxj=int_tabulated(wavej,intflux*resj)/int_tabulated(wavej,resj)*1.0e26
    
    intflux=10^(interpol(alog10(fnu),alog10(lambda*1.0e-4),alog10(waveh)))
    fluxh=int_tabulated(waveh,intflux*resh)/int_tabulated(waveh,resh)*1.0e26
    
    intflux=10^(interpol(alog10(fnu),alog10(lambda*1.0e-4),alog10(wavek)))
    fluxk=int_tabulated(wavek,intflux*resk)/$
      int_tabulated(wavek,resk)*1.0e26
      
    intflux=10^(interpol(alog10(fnu),alog10(lambda*1.0e-4),alog10(wave24)))
    flux24=int_tabulated(wave24,intflux*res24)/int_tabulated(wave24,res24)*1.0e26
    
    intflux=10^(interpol(alog10(fnu),alog10(lambda*1.0e-4),alog10(wave70)))
    flux70=int_tabulated(wave70,intflux*res70)/int_tabulated(wave70,res70)*1.0e26
    
    intflux=10^(interpol(alog10(fnu),alog10(lambda*1.0e-4),alog10(wave160)))
    flux160=int_tabulated(wave160,intflux*res160)/$
      int_tabulated(wave160,res160)*1.0e26

    intflux=10^(interpol(alog10(fnu),alog10(lambda*1.0e-4),alog10(1200.)))
    flux1200=intflux*1.0e26

    printf,1,'HIP ',strtrim(string(long(hip[i])),2),$
      ' Predicted MIPS Fluxes (mJy):',fluxj,fluxh,fluxk,flux24,flux70,flux160,flux1200
endfor
close,1

stop
end
