pro test_disk_benchmarks

restore,filename='disk_benchmark_dataHD146897.sav'
; Format:
; times[run,*] = [total_time, p1, p2, p3, p4, p5, p6, p7, p8, ttr1, ttr2]


; Plot histogram of total runtimes:
plothist,times[*,0],xhist,yhist


end