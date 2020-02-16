#!/usr/local/bin/gunuplot
unset key
set datafile separator ","
set term postscript eps enhanced color
set size 0.75,0.75
set size ratio -1
set border linewidth 1.5
#http://www.ss.scphys.kyoto-u.ac.jp/person/yonezawa/contents/program/gnuplot/index.html
#http://www2.ccs.tsukuba.ac.jp/Astro/Members/kirihara/gnuplot.html
set xlabel "{/Italic u}_{1} [{/Symbol-Oblique m}m]"
set ylabel "{/Italic r}_{2} [{/Symbol-Oblique m}m]" #set xlabel "{/Times-Italic x}_{2} [{/Symbol m}m]"
set pm3d map #set pm3d interpolate 4,4 map
set palette rgbformulae 33,13,10 #set palette defined (0 "blue", 1 "white", 2 "red")
set tics in #set tics out
set xtics 1.
set mxtics 5
set ytics 1.
set mytics 5
set grid xtics ytics lw 3
set ticscale 1.5
set xrange [-2.0:2.0]
set yrange [-5.0:5.0]
#
set cbrange[1:13]
set output 'epsi_image/Slice0_0.eps'
splot 'Slice0_0.dat' using ($1*0.1):($3*0.1):4 with pm3d
set output 'epsi_image/Slice0b_0.eps'
splot 'Slice0b_0.dat' using ($1*0.1):($3*0.1):4 with pm3d
#
set output 'epsi_image/Slice0_1.eps'
splot 'Slice0_1.dat' using ($1*0.1):($3*0.1):4 with pm3d
set output 'epsi_image/Slice0b_1.eps'
splot 'Slice0b_1.dat' using ($1*0.1):($3*0.1):4 with pm3d
#
set output 'epsi_image/Slice0_2.eps'
splot 'Slice0_2.dat' using ($1*0.1):($3*0.1):4 with pm3d
set output 'epsi_image/Slice0b_2.eps'
splot 'Slice0b_2.dat' using ($1*0.1):($3*0.1):4 with pm3d

