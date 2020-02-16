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
#set pm3d map #set pm3d interpolate 4,4 map
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
#set style line 1 pointtype 7 pointsize 0.075 linecolor 'red'
set style line 1 pointtype 7 pointsize 0.2 linecolor 'red'
set style line 2 pointtype 7 pointsize 0.2 linecolor 'green'
set style line 3 pointtype 7 pointsize 0.2 linecolor 'blue'
#
set output 'grid_image/Slice2_0g.eps'
plot 'Slice2_0.dat' using ($2*0.1):($3*0.1) linestyle 1
#
set output 'grid_image/Slice2_1g.eps'
plot 'Slice2_1.dat' using ($2*0.1):($3*0.1) linestyle 2
#
set output 'grid_image/Slice2_2g.eps'
plot 'Slice2_2.dat' using ($2*0.1):($3*0.1) linestyle 3
#
set xrange [-2.4:2.4]
#
set output 'grid_image/SliceB_0g.eps'
plot 'SliceB_0.dat' using ($2*0.1):($3*0.1) linestyle 1
set output 'grid_image/SliceT_0g.eps'
plot 'SliceT_0.dat' using ($2*0.1):($3*0.1) linestyle 1
#
set output 'grid_image/SliceB_1g.eps'
plot 'SliceB_1.dat' using ($2*0.1):($3*0.1) linestyle 2
set output 'grid_image/SliceT_1g.eps'
plot 'SliceT_1.dat' using ($2*0.1):($3*0.1) linestyle 2
#
set output 'grid_image/SliceB_2g.eps'
plot 'SliceB_2.dat' using ($2*0.1):($3*0.1) linestyle 3
set output 'grid_image/SliceT_2g.eps'
plot 'SliceT_2.dat' using ($2*0.1):($3*0.1) linestyle 3

