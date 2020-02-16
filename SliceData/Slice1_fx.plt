#!/usr/local/bin/gunuplot
unset key
set datafile separator ","
set term postscript eps enhanced color
set size 0.75,0.75
set size ratio -1
set border linewidth 1.5
#http://www.ss.scphys.kyoto-u.ac.jp/person/yonezawa/contents/program/gnuplot/index.html
#http://www2.ccs.tsukuba.ac.jp/Astro/Members/kirihara/gnuplot.html
set xlabel "{/Italic x}_{2} [{/Symbol-Oblique m}m]" #set xlabel "{/Times-Italic x}_{2} [{/Symbol m}m]"
set ylabel "{/Italic x}_{0} [{/Symbol-Oblique m}m]"
set pm3d map #set pm3d interpolate 4,4 map
set palette rgbformulae 33,13,10 #set palette defined (0 "blue", 1 "white", 2 "red")
set tics in #set tics out
set xtics 1.
set mxtics 2
set ytics 1.
set mytics 2
set grid xtics ytics lw 3
set ticscale 1.5
set xrange [-5.0:5.0]
set yrange [-2.0:4.0]
#
set cbrange[0:2]
set output 'fx_image/Slice1_0fx.eps'
splot 'Slice1_0.dat' using ($2*0.1):($4*0.1):7 with pm3d
set output 'fx_image/Slice1b_0fx.eps'
splot 'Slice1b_0.dat' using ($2*0.1):($4*0.1):7 with pm3d
#
set cbrange[0:2]
set output 'fx_image/Slice1_1fx.eps'
splot 'Slice1_1.dat' using ($2*0.1):($4*0.1):7 with pm3d
set output 'fx_image/Slice1b_1fx.eps'
splot 'Slice1b_1.dat' using ($2*0.1):($4*0.1):7 with pm3d
#
set cbrange[0:40]
set output 'fx_image/Slice1_2fx.eps'
splot 'Slice1_2.dat' using ($2*0.1):($4*0.1):7 with pm3d
set output 'fx_image/Slice1b_2fx.eps'
splot 'Slice1b_2.dat' using ($2*0.1):($4*0.1):7 with pm3d

