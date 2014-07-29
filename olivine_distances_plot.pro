pro olivine_distances_plot

readcol, 'olivine_distances.csv',F='a,f,f,f',Name,Temp,Dist,Loc,Diff

plothist,Diff,xhist,yhist,bin=7,xtitle='Difference',ytitle='Count'
print,'Average Difference: '+string(mean(Diff))
print,'Median Difference: '+string(median(Diff))

end