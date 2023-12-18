#preliminaries----------------------------------------------------#
set style line 1 lt 1 linecolor rgb "black" lw 3 pt 1 ps 1 
set style fill transparent solid 0.1
#-----------------------------------------------------------------#

num   = 30.0
xlow  = 0.0
xhigh = 5.0
width = (xhigh-xlow)/num
bin(x)=width*floor((x)/width)
set boxwidth width

#set xrange [xlow-0.02:xhigh]
set xrange [xlow-0.09:xhigh]
set xlabel 'r_{est} ({/Symbol \264} 10^{3})'


set yrange [:170]

plot "LensedCMB_noiseless/twocomp/Nside4/fwhm2220/r0.001/modified_delta_map_CMB_140_FGD_40_60_230_280_340.dat" us (bin($2*1.0e3)):(1.0) \
smooth frequency with boxes fs transparent solid 0.1 lw 3 lc rgb "blue" title '(40,60,230,280,340)' ,\
"LensedCMB_noiseless/twocomp/Nside4/fwhm2220/r0.001/modified_delta_map_CMB_140_FGD_40_60_230_340_400.dat" us (bin($2*1.0e3)):(1.0) \
smooth frequency with boxes fs transparent solid 0.1 lw 3 lc rgb "red" title '(40,60,230,340,400)' 

set parametric
set trange [0:170]
replot 1,t title '' with lines ls 1 lw 3 dt (10,10)

set terminal pdfcairo color enhanced size 6,5 font "Helvetica,22"
set output "Figs/fig3b_new.pdf"
replot
