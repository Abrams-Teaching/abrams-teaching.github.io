#!/usr/bin/gnuplot
set term gif
set out "hhd.gif"
unset border
unset xtics
unset ytics
set size square
unset key
set style data line

p "0.hhd" u 1:($2+6), "1.hhd" u 1:($2+4), "2.hhd" u 1:($2+2),"3.hhd" u 1:($2),\
"4.hhd" u ($1+2):($2+6), "5.hhd" u ($1+2):($2+4), "6.hhd" u ($1+2):($2+2),"7.hhd" u ($1+2):($2),\
"8.hhd" u ($1+4):($2+6), "9.hhd" u ($1+4):($2+4), "10.hhd" u ($1+4):($2+2),"11.hhd" u ($1+4):($2)
