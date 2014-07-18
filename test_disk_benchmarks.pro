; +
; NAME:
;  test_disk_benchmarks
;
; PURPOSE:
;  Implements analysis of computation time in the code disk_spectrum. 
;
; INPUTS:
;  NONE
;
; KEYWORDS:
;
; OUTPUTS:
;   plots
;
; AUTHORS:
;  Emil Christensen - chris2er@dukes.jmu.edu
;
; DISCLAIMER
;  This software is provided as is without any warranty whatsoever.
;  Permission to use, copy, modify, and distribute modified or
;  unmodified copies is granted, provided this disclaimer
;  is included unchanged.
;
; MODIFICATION HISTORY:
;  Written by EC (7/17/2014)
; -
; *************************************************** ;

pro test_disk_benchmarks

; *************************************************** ;
; Compute analytics for all of code

restore,filename='disk_benchmark_data_HD146897.sav'
; Format:
; times[run,*] = [total_time, p1, p2, p3, p4, p5, p6, p7, p8, ttr1, ttr2]
titles = ['Total Time','Setup','equiltemplookup','qlookup','integrand','Reform/Rebin','Blackbody','Matrix Mult','End','Transpose 1','Transpose 2']

; Stats
print, 'Average Completion Time: '+strcompress(string(mean(times[*,0])))
print, 'Estimated Total Run Time: '+string(round(total(times[*,0])))+' seconds'

plothist,times[*,0],xhist,yhist, bin=0.1, title=titles[0],xtitle='Time (s)', ytitle='counts'
;stop

; Plot histogram of total runtimes:
for i=1, n_elements(times[0,*])-1 do begin
  av = mean(times[*,i])
  print,string(titles[i])+' Average = '+string(av)
  if (av ne 0) then begin
    plothist,times[*,i],xhist,yhist, bin=1, title=titles[i],xtitle='%', ytitle='counts'
  endif
  ;stop
endfor

; Run 1:

;Average Completion Time:   0.42461147
;Estimated Total Run Time:  3397 seconds (~56 minutes)
;Setup Average =            0.0000000
;equiltemplookup Average =  16.499125
;qlookup Average =          8.8080000
;integrand Average =        0.0000000
;Reform/Rebin Average =     19.500375  --> bad
;Blackbody Average =        51.729750  --> main problem
;Matrix Mult Average =      3.3630000
;End Average =              0.0000000
;Transpose 1 Average =      0.0000000
;Transpose 2 Average =      16.010500

; Run 2

;Average Completion Time:  0.25574188
;Estimated Total Run Time: 2046 seconds (~34 minutes)
;Setup Average =           0.0000000
;equiltemplookup Average = 28.392625 --> double of last run
;qlookup Average =         16.510000 --> double of last run
;integrand Average =       0.0000000
;Reform/Rebin Average =    15.110125 --> Still bad
;Blackbody Average =       37.177125 --> This and Reform/Rebin = ~52% 
;Matrix Mult Average =     2.7301250
;End Average =             0.0000000
;Transpose 1 Average =     0.0000000
;Transpose 2 Average =     10.967750



; *************************************************** ;
; Compute analytics for intensive area (Reform/Rebin + Blackbody)

restore,filename='disk_benchmark_lines_HD146897.sav'
line_titles = ['Line 1','Line 2','Line 3','Line 4','Line 5','Line 6','Line 7','Line 8']

for i=0, n_elements(lines[0,*])-1 do begin
  av = mean(times[*,i])
  print,string(line_titles[i])+' Average = '+string(av)
  if (av ne 0) then begin
    plothist,lines[*,i],xhist,yhist, bin=1, title=line_titles[i],xtitle='% Time', ytitle='counts'
  endif
  ;stop
endfor

; Run 1

;Line 1 Average = 0.2557419 - wave_new1= (REBIN(lambda,NA*NR,NL))*1.d-4
;Line 2 Average = 0.0000000 - temp_t2= (REBIN(temp_t1,NA*NR,NL))
;Line 3 Average = 28.392625 - hnukT = ((h*c/k)/wave_new1)/temp_t2
;Line 4 Average = 16.510000 - tmp_p1 = (2.0*h*c)/(wave_new1)^3
;Line 5 Average = 0.0000000 - brightness = (tmp_p1*(1.d23/206265.0^2))/(exp(hnukT)-1.0)
;Line 6 Average = 15.110125 - subl = where(temp_t2 gt 1000.0)
;Line 7 Average = 37.177125 - if subl[0] ne -1 then brightness[subl] = 0.0
;Line 8 Average = 2.7301250 - brightness = transpose(brightness)


end