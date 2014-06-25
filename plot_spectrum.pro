pro plot_spectrum,name,wave,spec,specerr
;plot_spectrum,'test',NEW_WAVE,NEW_SPEC,NEW_SPECERR

set_plot,'PS'
    
device, filename ='plots/'+name+'_spectrum.ps',/COLOR,/HELVETICA,XSIZE=15,YSIZE=12.5 & !p.font =0
;loadct,39
;!p.background=16777215
    
;Make IRS Spectrum Plot w/ Errorbars
;TVLCT,[0,255,0,0],[0,0,255,0],[0,0,0,255]
    
plot,wave,spec,title=name,ystyle=1,psym=-3,xstyle=1,xtitle=('Wavelength (\mum)'),ytitle=('F_\nu (Jy)'),charthick=1, thick=1, xthick=2, ythick=2, charsize=1,color=0
;oplot,wave,spec,color=0
oploterr,wave,spec,specerr


;LEGEND,TEXTOIDL('\chi^{2}/d.o.f. : ')+ strtrim(string(tmp1,format='(f18.2)'),1), /top, /left,color=1,linestyle=lines[1]
    
device,/close
set_plot,'x'

end

