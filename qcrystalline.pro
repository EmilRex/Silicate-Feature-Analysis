pro qcrystalline, Qforsterite, Qenstatite
; read in files of optical constants
; generate Qext, Qabs values, as a function of grain size & wavelength
; three structures: Qenstatite and Qforsterite
; each structure contains:
; - agrain: 0.1 to 1mm, 100 bins per decade
; - lambda: read in from files of refraction indices
; - Qext: 2D array, value for each agrain & lambda, index in that order.
; - Qscat: " "

; from 0.1 microns to 1 mm
agrain = 0.1D * 10.0D^(dindgen(401)*0.01d) 
NA = n_elements(agrain)

; maximum grain size parameter allowed to mie_single is 12000.0
maxx = 12000.0

; FORSTERITE: Mg2 Si O4 -- magnesium end member of crystalline olivine
; tables given for Mg1.9 Fe0.1 Si O4 

; read in optical data
; 3 files: one for each orientation. 
; Assume equal probability of orientations, spheres
; Qext = (1/3)*( Qext_x + Qext_y + Qext_z) 
readcol, 'opacities/oliv_nk_x.txt', w_forst_x, n_forst_x, k_forst_x, $
         format='(F,F,F)', comment='#'
readcol, 'opacities/oliv_nk_y.txt', w_forst_y, n_forst_y, k_forst_y, $
         format='(F,F,F)', comment='#'
readcol, 'opacities/oliv_nk_z.txt', w_forst_z, n_forst_z, k_forst_z, $
         format='(F,F,F)', comment='#'

l_forst_x = reverse(1.0/w_forst_x)*1d4
l_forst_y = reverse(1.0/w_forst_y)*1d4
l_forst_z = reverse(1.0/w_forst_z)*1d4

n_forst_x = reverse(n_forst_x)
n_forst_y = reverse(n_forst_y)
n_forst_z = reverse(n_forst_z)

k_forst_x = reverse(k_forst_x)
k_forst_y = reverse(k_forst_y)
k_forst_z = reverse(k_forst_z)

l_forst = l_forst_x[ where(l_forst_x gt 0) ]
n_forst = [ [interpol(n_forst_x, l_forst_x, l_forst)] , $
          [interpol(n_forst_y, l_forst_y, l_forst)] , $
          [interpol(n_forst_z, l_forst_z, l_forst)] ]
k_forst = [ [interpol(k_forst_x, l_forst_x, l_forst)] , $
          [interpol(k_forst_y, l_forst_y, l_forst)] , $
          [interpol(k_forst_z, l_forst_z, l_forst)] ]

; hack to include small wavelength components:
l_forst = [0.01, 0.1, 0.2, l_forst]

n_sm = [ [replicate(n_forst[0,0],3)], $
         [replicate(n_forst[0,1],3)], $
         [replicate(n_forst[0,2],3)] ]
n_forst = [ n_sm, n_forst]   ; constant at small lambda

k_sm = [ [findgen(3)/3.0*k_forst[0,0]], $
         [findgen(3)/3.0*k_forst[0,1]], $
         [findgen(3)/3.0*k_forst[0,2]] ]
k_forst = [k_sm, k_forst]           ; linear at small lambda

; calculate Qforsterite
nl = n_elements(l_forst)
Qforsterite = {agrain:agrain, lambda:l_forst, $
               Qext:dblarr(NA,nl), Qscat:dblarr(NA,nl)}

print, 'Forsterite'
for i=0, nl-1 do begin
                                
; calculate large grain limit for each wavelength
   qext_lg = 0.0D
   qscat_lg = 0.0D
   for j = 0,2 do begin
      mie_single, maxx, complex(n_forst[i,j], -k_forst[i,j]), qext, qscat
      qext_lg = qext_lg + qext/3.0
      qscat_lg = qscat_lg + qscat/3.0
   endfor
   Qforsterite.Qext [*,i] = qext_lg
   Qforsterite.Qscat[*,i] = qscat_lg

; calculate over all grain sizes for each wavelength
   x = 2.0*!dpi*agrain/l_forst[i]
   smallgrains = where(x le maxx)
   qext = 0.0D
   qscat = 0.0D
   for j=0,2 do begin
      mie_single, x[smallgrains], complex(n_forst[i], -k_forst[i]), qe, qs
      qext = qext + qe/3.0
      qscat = qscat + qs/3.0
   endfor
   Qforsterite.Qext [smallgrains,i] = qext
   Qforsterite.Qscat[smallgrains,i] = qscat
   print, 'lambda = ', l_forst[i]
endfor


; ENSTATITE: crystalline form of pyroxene

; read in optical data
; 3 files: one for each orientation. 
; Assume equal probability of orientations, spheres
; Qext = (1/3)*( Qext_x + Qext_y + Qext_z) 
readcol, 'opacities/ens_p.lnk',  l_ens_x, n_ens_x, k_ens_x, $
         format='(F,F,F)'
readcol, 'opacities/ens_s1.lnk', l_ens_y, n_ens_y, k_ens_y, $
         format='(F,F,F)'
readcol, 'opacities/ens_s2.lnk', l_ens_z, n_ens_z, k_ens_z, $
         format='(F,F,F)'

l_ens_all = [l_ens_x, l_ens_y, l_ens_z]
ii = uniq(l_ens_all, sort(l_ens_all))
l_ens = l_ens_all[ii]
n_ens = [ [interpol(n_ens_x, l_ens_x, l_ens)] , $
          [interpol(n_ens_y, l_ens_y, l_ens)] , $
          [interpol(n_ens_z, l_ens_z, l_ens)] ]
k_ens = [ [interpol(k_ens_x, l_ens_x, l_ens)] , $
          [interpol(k_ens_y, l_ens_y, l_ens)] , $
          [interpol(k_ens_z, l_ens_z, l_ens)] ]

; hack to include small wavelength components:
l_ens = [0.01, 0.1, 0.2, l_ens]

n_sm = [ [replicate(n_ens[0,0],3)], $
         [replicate(n_ens[0,1],3)], $
         [replicate(n_ens[0,2],3)] ]
n_ens = [ n_sm, n_ens]   ; constant at small lambda

k_sm = [ [findgen(3)/3.0*k_ens[0,0]], $
         [findgen(3)/3.0*k_ens[0,1]], $
         [findgen(3)/3.0*k_ens[0,2]] ]
k_ens = [k_sm, k_ens]           ; linear at small lambda

; calculate Qenstatite
nl = n_elements(l_ens)
Qenstatite = {agrain:agrain, lambda:l_ens, $
              Qext:dblarr(NA,nl), Qscat:dblarr(NA,nl)}

print, 'Enstatite'
for i=0, nl-1 do begin
                                
; calculate large grain limit for each wavelength
   qext_lg = 0.0D
   qscat_lg = 0.0D
   for j = 0,2 do begin
      mie_single, maxx, complex(n_ens[i,j], -k_ens[i,j]), qext, qscat
      qext_lg = qext_lg + qext/3.0
      qscat_lg = qscat_lg + qscat/3.0
   endfor
   Qenstatite.Qext [*,i] = qext_lg
   Qenstatite.Qscat[*,i] = qscat_lg

; calculate over all grain sizes for each wavelength
   x = 2.0*!dpi*agrain/l_ens[i]
   smallgrains = where(x le maxx)
   qext = 0.0D
   qscat = 0.0D
   for j=0,2 do begin
      mie_single, x[smallgrains], complex(n_ens[i], -k_ens[i]), qe, qs
      qext = qext + qe/3.0
      qscat = qscat + qs/3.0
   endfor
   Qenstatite.Qext [smallgrains,i] = qext
   Qenstatite.Qscat[smallgrains,i] = qscat
   print, 'lambda = ', l_ens[i]
endfor

plot, [0],[0], xrange=[1,100], yrange=[1e-5,10], /xlog, /ylog, /nodata

for i=0,4 do begin
   myind = (NA-1)/4*i
   mycol = !p.color + i*255/5
   print, agrain[myind]
   oplot, Qforsterite.lambda, Qforsterite.Qext[myind,*]/agrain[myind], $
          lines=0, color=mycol
   oplot, Qforsterite.lambda, Qforsterite.Qscat[myind,*]/agrain[myind], $
          lines=1, color=mycol
   oplot, Qenstatite.lambda, Qenstatite.Qext[myind,*]/agrain[myind], $
          lines=2, color=mycol
   oplot, Qenstatite.lambda, Qenstatite.Qscat[myind,*]/agrain[myind], $
          lines=3, color=mycol
endfor

return
end


pro newqtables
qcrystalline, Qforsterite, Qenstatite
restore, 'qtables.sav'
save, file='qtables_withcrys.sav', Qastrosil, Qolivine, Qpyroxene,$
      Qforsterite, Qenstatite
return
end


pro qcrys_short, Qforsterite, Qenstatite
; calculate short wavelength components, from values of indices of 
; refractions tabulated by Huffman & Stapp
; INPUT: existing Qforsterite, Qenstatite tables
; OUTPUT: modified tables

; maximum grain size parameter allowed to mie_single is 12000.0
maxx = 12000.0

; FORSTERITE: Mg2 Si O4 -- magnesium end member of crystalline olivine
agrain = Qforsterite.agrain
NA = n_elements(agrain)

; read in optical data
readcol, 'opacities/HuffmanStapp1973/crysoliv.txt', $
         format='(F,F,F)', comment='#', lambda, nfor, kfor

; calculate Qforsterite
nl = n_elements(lambda)

Qextf  = dblarr(NA,nl)
Qscatf = dblarr(NA,nl)
qe = 0.0D
qs = 0.0D

print, 'Forsterite'
for i=0, nl-1 do begin
                                
; calculate large grain limit for each wavelength
   mie_single, maxx, complex(nfor[i], -kfor[i]), qe, qs
   Qextf[*,i]  = qe
   Qscatf[*,i] = qs

; calculate over all grain sizes for each wavelength
   x = 2.0*!dpi*agrain/lambda[i]
   smallgrains = where(x le maxx)
   mie_single, x[smallgrains], complex(nfor[i], -kfor[i]), qe, qs
   Qextf [smallgrains,i] = qe
   Qscatf[smallgrains,i] = qs
   print, 'lambda = ', lambda[i]

endfor

i1 = where(lambda le 1.99)
i2 = where(Qforsterite.lambda ge 1.99)
alll = [lambda[i1], Qforsterite.lambda[i2]]
Qext = [[Qextf[*,i1]], [Qforsterite.Qext[*,i2]]]
Qscat = [[Qscatf[*,i1]], [Qforsterite.Qscat[*,i2]]]
Qforsterite = {agrain:agrain, lambda:alll, $
               Qext:Qext, Qscat:Qscat}


; ENSTATITE: crystalline form of pyroxene
agrain = Qenstatite.agrain
NA = n_elements(agrain)

; read in optical data
readcol, 'opacities/HuffmanStapp1971_crysenst.txt', $
         format='(F,F,F)', comment='#', w, nens, kens

ii = sort(1.0/w)
lambda = (1.0/w)[ii]
nens = nens[ii]
kens = kens[ii]

; calculate Qenstatite
nl = n_elements(lambda)

Qexte  = dblarr(NA,nl)
Qscate = dblarr(NA,nl)
qe = 0.0D
qs = 0.0D

print, 'Enstatite'
for i=0, nl-1 do begin
                                
; calculate large grain limit for each wavelength
   mie_single, maxx, complex(nens[i], -kens[i]), qe, qs
   Qexte[*,i]  = qe
   Qscate[*,i] = qs

; calculate over all grain sizes for each wavelength
   x = 2.0*!dpi*agrain/lambda[i]
   smallgrains = where(x le maxx)
   mie_single, x[smallgrains], complex(nens[i], -kens[i]), qe, qs
   Qexte [smallgrains,i] = qe
   Qscate[smallgrains,i] = qs
   print, 'lambda = ', lambda[i]
endfor

i1 = where(lambda le 1.99)
i2 = where(Qenstatite.lambda ge 1.99)
alll = [lambda[i1], Qenstatite.lambda[i2]]
Qext = [[Qexte[*,i1]], [Qenstatite.Qext[*,i2]]]
Qscat = [[Qscate[*,i1]], [Qenstatite.Qscat[*,i2]]]
Qenstatite = {agrain:agrain, lambda:alll, $
               Qext:Qext, Qscat:Qscat}

plot, [0],[0], xrange=[1,100], yrange=[1e-5,10], /xlog, /ylog, /nodata

for i=0,4 do begin
   myind = (NA-1)/4*i
   mycol = !p.color + i*255/5
   print, agrain[myind]
   oplot, Qforsterite.lambda, Qforsterite.Qext[myind,*]/agrain[myind], $
          lines=0, color=mycol
   oplot, Qforsterite.lambda, Qforsterite.Qscat[myind,*]/agrain[myind], $
          lines=1, color=mycol
   oplot, Qenstatite.lambda, Qenstatite.Qext[myind,*]/agrain[myind], $
          lines=2, color=mycol
   oplot, Qenstatite.lambda, Qenstatite.Qscat[myind,*]/agrain[myind], $
          lines=3, color=mycol
endfor

return
end



pro qtables_shortcrys
restore, 'qtables_withcrys.sav'
QCRYS_SHORT, Qforsterite, Qenstatite
read_crys_opac, crystallineabs
save, file='qtables_withcrys.sav', Qastrosil, Qolivine, Qpyroxene,$
      Qforsterite, Qenstatite, crystallineabs
return
end
