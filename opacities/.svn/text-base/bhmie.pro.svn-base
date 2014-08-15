pro bhmie, x,refrel,nang,s1,s2,qext,qsca,qback,gsca
; WARNING your nang should be *smaller* by 1 than mxnang 
; because of IDL array indexinf (from 0).
;***********************************************************************
; Subroutine BHMIE is the Bohren-Huffman Mie scattering subroutine
;    to calculate scattering and absorption by a homogenous isotropic
;    sphere.
; Given:
;    X = 2*pi*a/lambda
;    REFREL = (complex refr. index of sphere)/(real index of medium)
;    NANG = number of angles between 0 and 90 degrees
;           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
;           if called with NANG<2, will set NANG=2 and will compute
;           scattering for theta=0,90,180.
; Returns:
;    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
;                                scatt. E perp. to scatt. plane)
;    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
;                                scatt. E parr. to scatt. plane)
;    QEXT = C_ext/pi*a**2 = efficiency factor for extinction
;    QSCA = C_sca/pi*a**2 = efficiency factor for scattering
;    QBACK = (dC_sca/domega)/pi*a**2
;          = backscattering efficiency [NB: this is (1/4*pi) smaller
;            than the "radar backscattering efficiency"; see Bohren &
;            Huffman 1983 pp. 120-123]
;    GSCA = <cos(theta)> for scattering
;
; Original program taken from Bohren and Huffman (1983), Appendix A
; Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
; in order to compute <cos(theta)>
; Converted to IDL by P. J. Flatau, Scripps Inst. Oceanography. UCSD 96/11/14
; 91/05/07 (BTD): Modified to allow NANG=1
; 91/08/15 (BTD): Corrected error (failure to initialize P)
; 91/08/15 (BTD): Modified to enhance vectorizability.
; 91/08/15 (BTD): Modified to make NANG=2 if called with NANG=1
; 91/08/15 (BTD): Changed definition of QBACK.
; 92/01/08 (BTD): Converted to full double precision and double complex
;                 eliminated 2 unneed lines of code
;                 eliminated redundant variables (e.g. APSI,APSI0)
;                 renamed RN -> EN = double precision N
;                 Note that DOUBLE COMPLEX and DCMPLX are not part
;                 of f77 standard, so this version may not be fully
;                 portable.  In event that portable version is
;                 needed, use src/bhmie_f77.f
; 93/06/01 (BTD): Changed AMAX1 to generic function MAX
; 96/11/14 (PJF): Converted to double precision (again) 
;                 added  declaration standarization (strong typing)
;                 some code "polishing" using Toolpack
;                 Converted to IDL
;***********************************************************************
;IDL example test code for "bhmie"
;pro test
; do not declare s1, s2 in the main code
; read,x,x1,x2
; refrel=dcomplex(x1,x2)
; bhmie, x,refrel,nang,s1,s2,qext,qsca,qback,gsca
; print,x, refrel, qext,qsca,qback,gsca
;stop
;end
;

;     .. Parameters ..
;      INTEGER mxnang, nmxx
;      PARAMETER (mxnang=1000,nmxx=150000)
       mxnang=1000
       nmxx=150000

;     .. Scalar Arguments ..
;      DOUBLE COMPLEX refrel
;      DOUBLE PRECISION gsca, qback, qext, qsca, x
;      INTEGER nang
;     .. Array Arguments ..
;      DOUBLE COMPLEX s1(2*mxnang-1), s2(2*mxnang-1)
; s1 and s2 are defined here not in the subroutine calling  "bhmie"
       s1=dcomplexarr(2*mxnang-1)
       s2=dcomplexarr(2*mxnang-1)
;     .. Local Scalars ..
;      DOUBLE COMPLEX an, an1, bn, bn1, drefrl, xi, xi1, y
;      DOUBLE PRECISION chi, chi0, chi1, dang, dx, en, fn, p, pii, psi,$
;                      psi0, psi1, theta, xstop, ymod
;      INTEGER j, jj, n, nmx, nn, nstop

;     .. Local Arrays ..
;     DOUBLE COMPLEX d(nmxx)    
      d=dcomplexarr(nmxx)      

;      DOUBLE PRECISION amu(mxnang), pi(mxnang), pi0(mxnang),$
;                      pi1(mxnang), tau(mxnang)

      amu=dblarr(mxnang)
      pi=dblarr(mxnang)
      pi0=dblarr(mxnang)
      pi1=dblarr(mxnang)
      tau=dblarr(mxnang)

;     .. Intrinsic Functions ..
;      INTRINSIC abs, atan, cos, dble, dcmplx, max, sin
      IF (nang GT mxnang) then begin
         print, ' error: nang > mxnang in bhmie'
         stop
      endif

      IF (nang LT 2) then nang = 2
;*** Obtain pi:
      pii = 4.D0*atan(1.D0)
      dx = x
;flatau I have converted everything to double so drefrl is really not needed
      drefrl = refrel
      y = x*drefrl
      ymod = abs(y)
;
;*** Series expansion terminated after NSTOP terms
;    Logarithmic derivatives calculated from NMX on down
      xstop = x + 4.D0*x^0.3333D0 + 2.D0
      nmx = max(xstop,ymod) + 15
;flatau  IDL fix (convert to integer) 
      nmx=fix(nmx)      
;      print, ' nmx ', nmx
; BTD experiment 91/1/15: add one more term to series and compare results
;      NMX=AMAX1(XSTOP,YMOD)+16
; test: compute 7001 wavelengths between .0001 and 1000 micron
; for a=1.0micron SiC grain.  When NMX increased by 1, only a single
; computed number changed (out of 4*7001) and it only changed by 1/8387
; conclusion: we are indeed retaining enough terms in series!
      nstop = xstop
;
      IF (nmx GT nmxx) THEN begin
          print, 'error: nmx > nmxx=', nmxx, ' for |m|x=', ymod
          STOP
      ENDIF
;*** Require NANG.GE.1 in order to calculate scattering intensities
      dang = 0.D0
      IF (nang GT 1) then dang = .5D0*pii/double(nang-1)
      for j=1, nang do begin ; DO 10 j = 1, nang
          theta = double(j-1)*dang
          amu(j) = cos(theta)
      endfor ; 10 CONTINUE
      for j=1, nang do begin ;DO 20 j = 1, nang
          pi0(j) = 0.D0
          pi1(j) = 1.D0
      endfor  ;20 CONTINUE
      nn = 2*nang - 1
      for j=1, nn do begin ;DO 30 j = 1, nn
          s1(j) = dcomplex(0.D0,0.D0)
          s2(j) = dcomplex(0.D0,0.D0)
      endfor ;30 CONTINUE
;
;*** Logarithmic derivative D(J) calculated by downward recurrence
;    beginning with initial value (0.,0.) at J=NMX
;
      d(nmx) = dcomplex(0.D0,0.D0)
      nn = nmx - 1
      for n=1, nn do begin ;DO 40 n = 1, nn
          en = nmx - n + 1
          d(nmx-n) = (en/y) - (1.D0/ (d(nmx-n+1)+en/y))
      endfor ;40 CONTINUE
;
;*** Riccati-Bessel functions with real argument X
;    calculated by upward recurrence
;
      psi0 = cos(dx)
      psi1 = sin(dx)
      chi0 = -sin(dx)
      chi1 = cos(dx)
      xi1 = dcomplex(psi1,-chi1)
      qsca = 0.D0
      gsca = 0.D0
      p = -1.D0
      for n=1, nstop do begin; DO 80 n = 1, nstop
          en = n
          fn = (2.D0*en+1.D0)/ (en* (en+1.D0))
; for given N, PSI  = psi_n        CHI  = chi_n
;              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
;              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
; Calculate psi_n and chi_n
          psi = (2.D0*en-1.D0)*psi1/dx - psi0
          chi = (2.D0*en-1.D0)*chi1/dx - chi0
          xi = dcomplex(psi,-chi)
;
;*** Store previous values of AN and BN for use
;    in computation of g=<cos(theta)>
          IF (n GT 1) THEN begin
              an1 = an
              bn1 = bn
          ENDIF
;
;*** Compute AN and BN:
          an = (d(n)/drefrl+en/dx)*psi - psi1
          an = an/ ((d(n)/drefrl+en/dx)*xi-xi1)
          bn = (drefrl*d(n)+en/dx)*psi - psi1
          bn = bn/ ((drefrl*d(n)+en/dx)*xi-xi1)
;
;*** Augment sums for Qsca and g=<cos(theta)>
          animag=double(an*complex(0.d0,-1.d0))
          bnimag=double(bn*complex(0.d0,-1.d0))
          qsca = qsca + (2.D0*en+1.D0)* (abs(an)^2+abs(bn)^2)
          gsca = gsca + ((2.D0*en+1.D0)/ (en* (en+1.D0)))* $
                (double(an)*double(bn)+animag*bnimag)
          IF (n GT 1) THEN begin
              an1imag=double(an1*complex(0d0,-1.d0))
              bn1imag=double(bn1*complex(0d0,-1.d0))
              animag=double(an*complex(0d0,-1.d0))
              bnimag=double(bn*complex(0d0,-1.d0))
              gsca = gsca + ((en-1.D0)* (en+1.D0)/en)*$
                    (double(an1)*double(an)+an1imag*animag+$
                    double(bn1)*double(bn)+bn1imag*bnimag)
          ENDIF
;
;*** Now calculate scattering intensity pattern
;    First do angles from 0 to 90
          for j=1, nang do begin ;DO 50 j = 1, nang
              jj = 2*nang - j
              pi(j) = pi1(j)
              tau(j) = en*amu(j)*pi(j) - (en+1.D0)*pi0(j)
              s1(j) = s1(j) + fn* (an*pi(j)+bn*tau(j))
              s2(j) = s2(j) + fn* (an*tau(j)+bn*pi(j))
          endfor ; 50     CONTINUE
;
;*** Now do angles greater than 90 using PI and TAU from
;    angles less than 90.
;    P=1 for N=1,3,...; P=-1 for N=2,4,...
          p = -p
          for j=1, nang-1 do begin ; DO 60 j = 1, nang - 1
              jj = 2*nang - j
              s1(jj) = s1(jj) + fn*p* (an*pi(j)-bn*tau(j))
              s2(jj) = s2(jj) + fn*p* (bn*pi(j)-an*tau(j))
          endfor ; 60     CONTINUE
          psi0 = psi1
          psi1 = psi
          chi0 = chi1
          chi1 = chi
          xi1 = dcomplex(psi1,-chi1)
;
;*** Compute pi_n for next value of n
;    For each angle J, compute pi_n+1
;    from PI = pi_n , PI0 = pi_n-1
          for j=1, nang do begin ; DO 70 j = 1, nang
              pi1(j) = ((2.D0*en+1.D0)*amu(j)*pi(j)- (en+1.D0)*pi0(j))/$
                      en
              pi0(j) = pi(j)
           endfor ;70     CONTINUE
           endfor ;   80 CONTINUE
;
;*** Have summed sufficient terms.
;    Now compute QSCA,QEXT,QBACK,and GSCA
      gsca = 2.D0*gsca/qsca
      qsca = (2.D0/ (dx*dx))*qsca
      qext = (4.D0/ (dx*dx))*double(s1(1))
      qback = (abs(s1(2*nang-1))/dx)^2/pii
      RETURN
      END

