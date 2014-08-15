pro crysoliv_short
; read in data thief output
; create new file of optical constants 
readcol, 'crysoliv_n.txt', ln, nopt, format='(F,F)', comment='#'
readcol, 'crysoliv_k.txt', lk, kopt, format='(F,F)', comment='#'
lk = [0.0, lk]
kopt = [0.0, kopt]
wmin = 0.25
wmax = 12.5
dw = 0.05
w = wmin + dw*findgen( (wmax-wmin)/dw + 1)
lambda = reverse(1.0/w)
n = interpol(nopt, ln, 1.0/lambda)
k = interpol(kopt, lk, 1.0/lambda)
openw, 1, 'crysoliv.txt'
printf, 1, "# lambda  n   k"
printf, 1, transpose([ [lambda], [n], [k]])
close, 1
return
end
