#preliminaries----------------------------------------------------#
set style line 1 lt 1 linecolor rgb "black" lw 3 pt 1 ps 1 
set style fill transparent solid 0.1
#-----------------------------------------------------------------#

num   = 60.0
xlow  = -0.25
xhigh = 2.5
width = (xhigh-xlow)/num
bin(x)=width*floor((x)/width)
set boxwidth width

set xrange [xlow:xhigh]
set xlabel 'r_{est} ({/Symbol \264} 10^{3})'

set yrange [:550]

plot "LensedCMB_noiseless/twocomp/Nside4/fwhm2220/r0.001/nulltest.table" us (bin($2*1.0e3)):(1.0) \
smooth frequency with boxes fs transparent solid 0.1 lw 3 lc rgb "blue" title '100GHz'

replot "LensedCMB_noiseless/twocomp/Nside4/fwhm2220/r0.001/nulltest.table" us (bin($5*1.0e3)):(1.0) \
smooth frequency with boxes fs transparent solid 0.1 lw 3 lc rgb "red" title '140GHz'

replot "LensedCMB_noiseless/twocomp/Nside4/fwhm2220/r0.001/nulltest.table" us (bin(($5-$2)*1.0e3)):(1.0) \
smooth frequency with boxes fs transparent solid 0.1 lw 3 lc rgb "black" title '100GHz-140GHz'



set parametric
set trange [0:600]
replot 1,t title '' with lines ls 1 dt (10,10)
replot 0,t title '' with lines ls 1 dt (10,10)

set terminal pdfcairo color enhanced size 6,5 font "Helvetica,22"
set output 'Figs/fig4b_new.pdf' 

replot
