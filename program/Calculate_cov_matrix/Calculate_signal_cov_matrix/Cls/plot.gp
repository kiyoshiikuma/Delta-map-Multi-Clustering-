#preliminaries----------------------------------------------------#
reset
set fontpath '/opt/teTeX/share/texmf-dist/fonts/type1/bluesky/cm'
set tics font "Times-Roman,26"
set key right top
set pointsize 3
set style line 1 lt 1 linecolor rgb "black" lw 1 pt 1 ps 1 
set style line 2 lt 2 linecolor rgb "black" lw 3 pt 1 ps 1 
set style line 3 lt 3 linecolor rgb "blue" lw 3 pt 1 ps 1 
set style line 4 lt 5 linecolor rgb "red" lw 2 pt 1 ps 1 
set style line 5 lt 4 linecolor rgb "green" lw 2 pt 1 ps 1 
set style line 6 lt 6 linecolor rgb "magenta" lw 2 pt 1 ps 1 
#-----------------------------------------------------------------#
set ylabel '{/CMMI10 \140}({/CMMI10 \140}+{1}){/CMMI10 C}_{/CMMI10 \140}/ {2}{/Symbol p} (TT)' offset 1 font ",30" 
set xlabel '{/CMMI10 \140}' font ",36"

set xrange [:1000]
set log 

plot 'lblike_tensCls.dat' us 1:3 title 'kk2011' w l
replot 'PlanckTTlowP_tensCls.dat' us 1:3 title 'PlanckTTlowP' w l
replot 'planck15lens_tensCls.dat' us 1:3 title 'planck15lens' w l

set term post eps enhanced color "Times-Roman" 26 fontfile "cmsy10.pfb" fontfile "cmmi10.pfb" fontfile "cmr12.pfb" fontfile "cmu10.pfb" fontfile "cmr8.pfb" fontfile "cmti10.pfb"
set output 'test.eps'
replot

set term qt
