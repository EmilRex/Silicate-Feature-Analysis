pro diskspectrum, rin, rout, amin, amax, Teff, diskmass, $ 
                     folivine, fcrystalline, fforst, $
                  rlaw=rlaw, alaw=alaw, verbose=verbose, $
                  lambda, flux, $
                  ring=ring, mie=mie
; rin, rout in units of stellar radius
; if ring keyword parameter set, then rin is ring position, rout is width
; amin, amax in units of microns
; Teff = effective temperature of star
; diskmass = normalized disk mass 
; folivine = olivine fraction
; olivine fraction = folivine*(1-fcrystalline)
; pyroxene fraction = (1-folivine)*(1-fcrystalline)
; forsterite fraction = fforst*fcrystalline
; enstatite fraction = (1-fforst)*fcrystalline
; lambda = input wavelengths at which to look up spectrum
; OUTPUTS:
; lambda, flux

;COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, crystallineabs

COMMON disk_benchmarking, run, times
print, run


t0 = systime(1)
;bulk density of grains in g/cm^3
rho_s = 3.3
AU_in_cm = 1.496e13

if not keyword_set(rlaw) then rlaw=1.5
if not keyword_set(alaw) then alaw=3.5

; number of grains per surface area per grain size bin
; = N(a) = N_0 (a/a0)^{-alaw} (r/r0)^{-rlaw}
; then total mass = 
; = int_rin^rout 2*pi*r dr int_amin^amax da 4/3 pi a^3 rho_s N(a)
; = 8*pi^2*rho_s/3 * int_rin^rout r dr int_amin^amax da a^3 
;                              * N_0 (a/a0)^{-alaw} (r/r0)^{-rlaw}
; = 8*pi^2*rho_s/3 *N_0*a0^{alaw}*r0^{rlaw} * 
;       int_rin^rout r^{-0.5} dr int_amin^amax da a^{-0.5} 
; = 8*pi^2*rho_s/3 *N_0*a0^{alaw}*r0^{rlaw} * 
;       2.0*[rout^0.5 - rin^0.5] * 2.0*[amax^0.5 - amin^0.5]
; = 32*pi^2*rho_s/3 *N_0*a0^4*r0^2 * 
;       [(rout/r0)^0.5 - (rin/r0)^0.5] * [(amax/a0)^0.5 - (amin/a0)^0.5]
; so N_0 = 3*diskmass/ { 32*pi^2*rho_s * a0^4*r0^2 * 
;                        [(rout/r0)^0.5 - (rin/r0)^0.5] * 
;                        [(amax/a0)^0.5 - (amin/a0)^0.5] }
; a0 = 1 micron, rho_s in g/cm
; N_0*a0*r0^2 = 3*diskmass/ { 32*pi^2*rho_s*(1e-4)^3 *
;                        [(rout/r0)^0.5 - (rin/r0)^0.5] * 
;                        [(amax/a0)^0.5 - (amin/a0)^0.5] }
; THIS EQUATION ONLY WORKS FOR ALAW=3.5
; norm = 3.0*diskmass/( 32.0*!dpi^2*rho_s*(1e-4)^3 * $
;                      (rout^0.5 - rin^0.5) * (amax^0.5 - amin^0.5) )

; sigma T^4 = !pi * int B_nu dnu
; luminosity = 4!pi R^2 sigma T^4

dloga = 0.02d
NA = ceil(alog10(amax/amin)/dloga)+1
dloga = alog10(amax/amin)/(NA-1.0d)
aall = amin*10.0D^(dindgen(NA)*dloga)
dla = replicate(dloga*alog(10.0D), NA)
dla[0] = 0.5*dla[0]
dla[NA-1] = dla[0]

if keyword_set(ring) then begin
   radius = rin
   width = rout
   rin = (radius - 4.0*width) > 5.0
   rout = radius + 4.0*width
endif 


if (abs(rout/rin-1.0) le 1e-3) then begin 
; single radius ring, no width'
   NR = 1
   dlr = 1.0
   rall = [rin]
   sigma = 1.0/(2.0*!dpi*rall^2)
endif else begin
   dlogr = 0.02d
   NR = ceil(alog10(rout/rin)/dlogr) + 1
   dlogr = alog10(rout/rin)/(NR-1.0d)
   rall = rin*10.0D^(dindgen(NR)*dlogr)
   dlr = replicate(dlogr*alog(10.0D), NR)
   dlr[0] = 0.5*dlr[0]
   dlr[NR-1] = dlr[0]
; set surface density profile

  if keyword_set(ring) then $
      sigma = exp(-0.5*(rall-radius)^2/width^2) $
   else $
      sigma = rall^(-rlaw)
endelse

t1 = systime(1)


; get temperatures - don't include crystalline silicates
equiltemplookup, Teff, aall, rall, tempr, $
        folivine, 0.0, 0.0

t2 = systime(1)


;if keyword_set(mie) then $
;;;; get Qabs, from Mie calculation
   qlookup, aall, lambda, folivine, fcrystalline, fforst, qabsall 
;else begin

t3 = systime(1)


;;;; use powder data for crystalline silicates
;   qlookup, aall, lambda, folivine, 0.0, 0.0, qabsall
;;;; crystalline silicates
;   abscoeff = crystallineopacities(lambda, reform(tempr), fforst, /grid)
;endelse


; calculate spectrum
flux = dblarr(n_elements(lambda))
NL = n_elements(lambda)
if nl le 1 then print, 'OH WOE!'

rintegrand = 2.0*!dpi * rall * sigma*dlr*rall      ; multiply by  R_sun^2
crosssection = !dpi*(aall*1e-4)^2 * aall^(-alaw)*dla*aall ; Multiply by 1e-4 for cm


; integrate to get normalization 
totsigma = total(rintegrand)
totmass = 4.0/3.0 * rho_s * total(aall*1e-4 * crosssection) * totsigma
amorphnorm = diskmass/totmass
crysnorm = diskmass/totsigma/total(crosssection)
;stop
;if keyword_set(verbose) then begin
;   print, 'amorphnorm = ', amorphnorm
;   if not keyword_set(ring) then $
;      print, 'rlaw=', rlaw
;endif


   c  = 2.998d10
   h = 6.626d-27
   k = 1.381d-16

t4 = systime(1)

ttrans_1_1 = systime(1)
   temp_t1 =reform(transpose(tempr),NA*NR) ; Not an issue
ttrans_1_2 = systime(1)

   wave_new1= (REBIN(lambda,NA*NR,NL))*1.d-4
   temp_t2= (REBIN(temp_t1,NA*NR,NL))

t5 = systime(1)

   hnukT = ((h*c/k)/wave_new1)/temp_t2
   tmp_p1 = (2.0*h*c)/(wave_new1)^3
   brightness = (tmp_p1*(1.d23/206265.0^2))/(exp(hnukT)-1.0)

    ; include dust sublimation
    subl = where(temp_t2 gt 1000.0)
    if subl[0] ne -1 then brightness[subl] = 0.0

;    tmp_p1 = REBIN(rintegrand,1,n_elements(rintegrand)*n_elements(tempr[*,0]))
;    rintegrand2 = REBIN(tmp_p1,n_elements(tmp_p1),n_elements(lambda))
;  brightness2 = reform(brightness,NR,NA,NL)

ttrans_2_1 = systime(1)
  brightness = transpose(brightness) ; Takes ~10-15% of time - alternatives?
ttrans_2_2 = systime(1)

;  brightness2 = (brightness)*abscoeff
  tmp_p1 = FLTARR(NA,NL) 

t6 = systime(1)

for i=0,NA-1 do begin

    intFdr = (matrix_multiply(brightness[*,0+NR*(i):NR-1+NR*(i)], rintegrand)*qabsall[i,*])*(crosssection[i]*amorphnorm)
;    intFdr2 = matrix_multiply(brightness2[*,0+NR*(i):NR-1+NR*(i)], rintegrand)*(crosssection[i]*crysnorm)
    tmp_p1[i,*] = intFdr;*(1.0-fcrystalline)+intFdr2*fcrystalline

;stop 
endfor

t7 = systime(1)

    flux = total(tmp_p1,1)
;PRINT, SYSTIME(1) - T4,   ' Seconds - A2c2 '
flux = flux/AU_in_cm^2

t8 = systime(1)

; *************************************************** ;
; Do benchmarking

total_time = (t8-t0)
print,'Total: ',total_time
print,''

p1 = round(100*(t1-t0)/(t8-t0))
print,'Setup: '+string(p1)+'%'

p2 = round(100*(t2-t1)/(t8-t0))
print,'equiltemplookup: '+string(p2)+'%'

p3 = round(100*(t3-t2)/(t8-t0))
print,'qlookup: '+string(p3)+'%'

p4 = round(100*(t4-t3)/(t8-t0))
print,'integrand: '+string(p4)+'%'

p5 = round(100*(t5-t4)/(t8-t0))
print,'Reform/Rebin: '+string(p5)+'%'

  ttr1=round(100*(ttrans_1_2-ttrans_1_1)/(t8-t0))
  print,'Transform 1: '+string(ttr1)+'%'

p6 = round(100*(t6-t5)/(t8-t0))
print,'Blackbody: '+string(p6)+'%'

  ttr2=round(100*(ttrans_2_2-ttrans_2_1)/(t8-t0))
  print,'Transform 2: '+string(ttr2)+'%'

p7 = round(100*(t7-t6)/(t8-t0))
print,'Matrix Mult For: '+string(p7)+'%'

p8 = round(100*(t8-t7)/(t8-t0))
print,'End: '+string(p8)+'%'

times[run,*] = [total_time, p1, p2, p3, p4, p5, p6, p7, p8, ttr1, ttr2]

run = run+1


;    if keyword_set(mie) then $
;       flux[i] = amorphflux $
;    else begin

;t1 = systime(1)
;if keyword_set(verbose) then begin
;   plot, lambda, flux, /xlog, /ylog
;   print, 'run time (new):', t1-t0
;endif

return
end










;pro compareringfits
;COMMON grainprops, Qastrosil, Qolivine, Qpyroxene, Qenstatite, Qforsterite, $
;   crystallineabs
;COMMON stellarprops, temptable, folivine, effectiveTemp, lambdastar, fluxstar

;;;restore, 'qtables_withcrys.sav' ; qastrosil, qolivine, qpyroxene
;;;restore, '../modelgrids/Teff6000grains.sav'

;lambda = [qastrosil.lambda[where(qastrosil.lambda lt qolivine.lambda[0])], $
;          qolivine.lambda]
;lambda = lambda[where(lambda ge 1.0)]
;lambda2 = lambda

;Teff = effectiveTemp
;diskspectrum, 100.0, 0.0, 5.0, 5.1, Teff, 0.1, 0.0, 0.0, $
;              alaw=3.5, lambda2, flux2 , /verb

;lambda1 = lambda2
;diskspectrum, 100.0, 1.0, 5.0, 5.1, Teff, 0.1, 0.0, 0.0, $
;              alaw=3.5, lambda1, flux1, $
;              /ring, /verb

;
;ii = interpol(findgen(n_elements(temptable[0].agrain)), temptable[0].agrain, 5.0)
;jj = interpol(findgen(n_elements(temptable[0].radii)), temptable[0].radii, 100.0)
;mytemp = interpolate(temptable[0].temperatures, ii, jj)
;onegr = modelonegrain(lambda1, [mytemp, 5.0, 0.0,0.0, 1.0])*0.002

;ymin = min([flux1,flux2], max=ymax)
;plot, lambda1, flux1, /xlog,/ylog, yrange=[ymin,ymax],lines=0
;oplot, lambda2, flux2, lines=1
;print, (flux1-flux2)/flux1
;oplot, lambda1, onegr, lines=2

;return
;end
