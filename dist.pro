FUNCTION tempfunc, r, qabs
  ;input r (array) and return tempfunc (array)
  ;estimates the amorphous olivine temperature as a function of distance.
  COMMON FUNC_LAM, lambda, dlambda
  
  h=6.63d-27
  c=3.0d10
  k=1.38d-16
  
  blambda=2.0*h*(c)^2/lambda^5/ $
    (exp(h*c/lambda/k/10000.0)-1.0)

  lhs = (2.5*6.96e10/r)^2*total(qabs*blambda*dlambda) ; get rid of first stellar term. Use the kurzz stellar atmosphere model instead of blamda  

  t=3.0*(findgen(1000) + 1.0)
  n=size(qabs,/n_elements)
  blambdad=dblarr(1000,n)
  rhs=dblarr(1000)

  for j=0,999 do begin
    blambdad[j,*]=2.0*h*(c)^2/lambda^5/ $
      (exp(h*c/lambda/k/t[j])-1.0)
    rhs[j] = total(4.0*qabs*blambdad[j,*]*dlambda)
  endfor
  
  tdust=interpol(t,rhs,lhs,/spline) ; keep

  return,tdust
end

; *************************************************** ;
pro dist
; Written by CC and copied by EC on 2/3/14
; Use to adapt modelonegrain.pro calculations
; Do for each graintype


readcol,'qabs/qabs_pyrmg50_1p5.txt',wave,qabs
;readcol,'bhmie/qabs_olivine5_1p7.dat',wave,qabs

tlambda=wave*1.0e-4 ; convert from microns to cm
r=findgen(500)*1.5d13
tlambdap=[tlambda,0,0]
tlambdam=[0,0,tlambda]
tdlambda=(tlambdap-tlambdam)/2.0
n=size(wave,/n_elements)
tdlambda[n]=(tlambda[n-1]-tlambda[n-2])/2.0

COMMON FUNC_LAM, lambda, dlambda
lambda=wave*1.0e-4
dlambda=tdlambda[1:n]

tdust=(dindgen(500)+1)
lam1=24.0d-4
lam2=70.0d-4
h=6.63d-27
c=3.0d10
k=1.38d-16
pi=3.14159

estdust = 100
print,'Estimated Dust Temperature:',estdust

tgr=tempfunc(r,qabs)
estdist=interpol(r,tgr,estdust)
print,'Estimated Dust Distance',estdist/1.5d13

stop
end