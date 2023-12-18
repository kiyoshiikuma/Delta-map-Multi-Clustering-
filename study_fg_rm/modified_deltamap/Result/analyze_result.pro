infile=dialog_pickfile(filter='*.dat')
r_input=0.0
readcol,infile,like,r,rest

mean_r = MEAN(r)
mean_like = MEAN(like)

var_r = VARIANCE(r)
var_like = VARIANCE(like)

sortedindex_r=SORT(r)
upper95=r[sortedindex_r[50]]
lower95=r[sortedindex_r[950]]

upper68=r[sortedindex_r[320]]
lower68=r[sortedindex_r[680]]

;upper=1.0
;lower=1.0

print,'    mean_r,      sqrt(var_r),   upper95,   lower95,  upper68  (if # of sample is 1000)'
print,mean_r,sqrt(var_r),upper95,lower95,upper68,lower68

;outfile = infile+'_analyzed'
;openw, lun, outfile, /get_lun
;printf,lun,r_input,mean_r,sqrt(var_r),sqrt(var_r)/sqrt(100.0)
;free_lun,lun
end
