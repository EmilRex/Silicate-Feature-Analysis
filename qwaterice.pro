; ****************************************************************************************************** ;
pro qwaterice
; ****************************************************************************************************** ;

; from 0.1 microns to 1 mm
agrain = 0.1D * 10.0D^(dindgen(401)*0.01d)
NA = n_elements(agrain)
; maximum grain size parameter allowed to mie_single is 12000.0
maxx = 12000.0

; read in optical data
readcol, 'IOP_2008_ASCIItable.dat', lwater_ice, m_re, m_im, format='(F,F,F)'

; calculate Qwater_ice
nl = n_elements(lwater_ice)
Qwater_ice = {agrain:agrain, lambda:lwater_ice, $
  Qext:dblarr(NA,nl), Qscat:dblarr(NA,nl),$
  logQext:dblarr(NA,nl), logQscat:dblarr(NA,nl)}
  
print, 'Water_ice'
for i=0, nl-1 do begin
  mie_single, maxx, complex(m_re[i], -m_im[i]), qext_lg, qscat_lg
  Qwater_ice.Qext [*,i] = qext_lg
  Qwater_ice.Qscat[*,i] = qscat_lg
  Qwater_ice.logQext[*,i] = alog(qext_lg)
  Qwater_ice.logQscat[*,i] = alog(qscat_lg)
  
  x = 2.0*!dpi*agrain/lwater_ice[i]
  smallgrains = where(x le maxx)
  mie_single, x[smallgrains], complex(m_re[i], -m_im[i]), qext, qscat
  Qwater_ice.Qext [smallgrains,i] = qext
  Qwater_ice.Qscat[smallgrains,i] = qscat
  Qwater_ice.logQext[smallgrains,i] = alog(qext)
  Qwater_ice.logQscat[smallgrains,i] = alog(qscat)
  print, 'lambda = ', lwater_ice[i]
endfor

save, file='qwaterice.sav', Qwater_ice

end

; ****************************************************************************************************** ;
pro newqtables_waterice
; ****************************************************************************************************** ;

restore, 'qwaterice.sav'
restore, 'qtables.sav'
save, file='qtables_withwaterice.sav', Qastrosil, Qolivine, Qpyroxene,$
Qforsterite, Qenstatite, Qwater_ice
return

end