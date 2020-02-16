#!/usr/local/bin/gunuplot
unset key
set datafile separator ","
set term postscript eps enhanced color
set size 0.75,0.75
set size ratio -1
set border linewidth 1.5
#http://www.ss.scphys.kyoto-u.ac.jp/person/yonezawa/contents/program/gnuplot/index.html
#http://www2.ccs.tsukuba.ac.jp/Astro/Members/kirihara/gnuplot.html
set xlabel "{/Italic r}_{0} [{/Symbol-Oblique m}m]" #set xlabel "{/Times-Italic x}_{2} [{/Symbol m}m]"
set ylabel "{/Italic u}_{1} [{/Symbol-Oblique m}m]"
set pm3d map #set pm3d interpolate 4,4 map
set palette rgbformulae 33,13,10 #set palette defined (0 "blue", 1 "white", 2 "red")
set tics in #set tics out
set xtics 1.
set mxtics 5
set ytics 1.
set mytics 5
set grid xtics ytics lw 3
set ticscale 1.5
set xrange [-1.2:1.2]
set yrange [-1.6:1.6]
#
set cbrange[0:10]
set output 'fx_image/Slice2_0fx.eps'
splot 'Slice2_0.dat' using ($2*0.1):($3*0.1):6 with pm3d
#
set cbrange[0:10]
set output 'fx_image/Slice2_1fx.eps'
splot 'Slice2_1.dat' using ($2*0.1):($3*0.1):6 with pm3d
#
set cbrange[0:14]
set output 'fx_image/Slice2_2fx.eps'
splot 'Slice2_2.dat' using ($2*0.1):($3*0.1):6 with pm3d
#
set xrange [-2.4:2.4]
set yrange [-1.6:1.6]
#
set cbrange[0:30]
set output 'fx_image/SliceB_0fx.eps'
splot 'SliceB_0.dat' using ($2*0.1):($3*0.1):6 with pm3d
set output 'fx_image/SliceT_0fx.eps'
splot 'SliceT_0.dat' using ($2*0.1):($3*0.1):6 with pm3d
#
set cbrange[0:10]
set output 'fx_image/SliceB_1fx.eps'
splot 'SliceB_1.dat' using ($2*0.1):($3*0.1):6 with pm3d
set output 'fx_image/SliceT_1fx.eps'
splot 'SliceT_1.dat' using ($2*0.1):($3*0.1):6 with pm3d
#
set cbrange[0:5]
set output 'fx_image/SliceB_2fx.eps'
splot 'SliceB_2.dat' using ($2*0.1):($3*0.1):6 with pm3d
set output 'fx_image/SliceT_2fx.eps'
splot 'SliceT_2.dat' using ($2*0.1):($3*0.1):6 with pm3d

