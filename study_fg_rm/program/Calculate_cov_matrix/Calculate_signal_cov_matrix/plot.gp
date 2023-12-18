#preliminaries----------------------------------------------------#
reset
set fontpath '/opt/teTeX/share/texmf-dist/fonts/type1/bluesky/cm'
set tics font "Times-Roman,26"
set key right top
set pointsize 3
set style line 1 lt 1 linecolor rgb "black" lw 3 pt 1 ps 1 
set style line 2 lt 2 linecolor rgb "red" lw 3 pt 2 ps 1 
set style line 3 lt 3 linecolor rgb "blue" lw 3 pt 1 ps 1 
set style line 4 lt 5 linecolor rgb "green" lw 2 pt 1 ps 1 
set style line 5 lt 4 linecolor rgb "cyan" lw 2 pt 1 ps 1 
set style line 6 lt 6 linecolor rgb "magenta" lw 2 pt 1 ps 1 
#-----------------------------------------------------------------#

r=1.e-3
set title 'r=10^{-3} (N_{side} = 8)'
sig_noise = 0.2
pi = 3.14159

set ylabel 'C_l'
set xlabel 'l'

set log y
set format y "10^{%L}"

plot 'signal_cl_scalar.txt' us 1:2 w l ls 1 title 'scalar EE'
replot 'signal_cl_tensor.txt' us 1:($2*r) w l ls 2 title 'tensor EE'
replot 'signal_cl_tensor.txt' us 1:($3*r) w l ls 3 title 'tensor BB'
replot 'signal_cl_scalar_lens.txt' us 1:($3) w l ls 4 title 'lens BB'
replot ((pi/10800.0)*0.2)**2.0 w l title '0.2 muK'
replot ((pi/10800.0)*0.4)**2.0 w l title '0.4 muK'
replot '/home/ichiki/Collab/LiteBird/SimulateDiagEEnoise/DiagEEnoise/Nside8/cltest/1.txt'\
us 1:3 w l title 'diag noise simulation' 

set term post eps enhanced color "Times-Roman" 26 fontfile \
"cmsy10.pfb" fontfile "cmmi10.pfb" fontfile "cmr12.pfb" \
fontfile "cmu10.pfb" fontfile "cmr8.pfb" fontfile "cmti10.pfb"

set output 'test.eps'

replot

set out
set term X11
reset
