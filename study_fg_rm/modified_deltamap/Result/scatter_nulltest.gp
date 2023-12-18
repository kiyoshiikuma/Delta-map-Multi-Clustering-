#preliminaries----------------------------------------------------#
reset
set fontpath '/opt/teTeX/share/texmf-dist/fonts/type1/bluesky/cm'
set tics font "Times-Roman,28"
set key right top
set pointsize 3
set style line 1 lt 1 linecolor rgb "black" lw 3 pt 1 ps 1 
set style line 2 lt 2 linecolor rgb "red" lw 3 pt 2 ps 1 
set style line 3 lt 3 linecolor rgb "blue" lw 3 pt 1 ps 1 
set style line 4 lt 5 linecolor rgb "green" lw 2 pt 1 ps 1 
set style line 5 lt 4 linecolor rgb "cyan" lw 2 pt 1 ps 1 
set style line 6 lt 6 linecolor rgb "magenta" lw 2 pt 1 ps 1 
#-----------------------------------------------------------------#
#set size square

set xrange [0:0.003]
set yrange [0:0.003]

set ylabel 'r_{140}'
set xlabel 'r_{100}'


plot 'LensedCMB_noiseless/twocomp/Nside4/fwhm2220/r0.001/nulltest.table' us 2:5 title ''

replot x title '' with lines ls 1 dt (10,10)
#replot 1,t title '' with lines ls 1 dt (10,10)

set term post eps enhanced color "Times-Roman" 26 fontfile \
"cmsy10.pfb" fontfile "cmmi10.pfb" fontfile "cmr12.pfb" \
fontfile "cmu10.pfb" fontfile "cmr8.pfb" fontfile "cmti10.pfb"

set output 'Figs/scatter_nulltest.eps' 

replot

set out
set term X11
reset
