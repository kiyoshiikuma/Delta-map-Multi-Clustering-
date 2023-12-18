infile = 'LensedCMB_noiseless/twocomp/Nside4/fwhm2220/r0.001/nulltest.table2'
readcol,infile,r140,a,r100

r = r100-r140
var = variance(r)
mean = mean(r)
print,mean,sqrt(var)

end
