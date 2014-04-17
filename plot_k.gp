## the bash command is
# $ gnuplot -e "variable=value" -p plot_kp.gp

##### command line -e input ##################

#if (!exists("pngIO")) pngIO = 0 ; 
#print 'pngIO = '.pngIO

##### SETTING ################################
filename = "paraOut_16O_JA1.0_JB1.5_angk120_angNN70_Sp10.0.dat"

##### Plotting script ########################
reset #reset all graphic setting
print '============================================'
print 'file : ',filename
pngName = 'compare'.filename[9:strlen(filename)-4].'.png'

set term x11 size 900,600
set grid   #set graphic grid
set style data lines
set key autotitle columnhead
#set xrange [0:180]
#unset xrange
#unset xtics

#-------------- set mulit plot <row>, <col>
set multiplot layout 2, 1 title filename
set tmargin 0
#set bmargin 0
set lmargin 10
set rmargin 3
set grid

#--- 1st plot
set origin 0,0.47
set size 1, 0.47
unset xlabel
set ylabel 'Pn000'
#set xtics 0, 20, 180
set xrange[32:40]
set yrange[50:58]

# the sytex (A?B:C) is read as (if A then B else C)
plot  0 lc rgb "#000000", \
     filename using 7:(-1*$9) lc rgb "#FF0000" lw 2, \
           "" using 7:(-1*$9) lc rgb "#FF6600" lw 2, \
           "" using 7:(-1*$9) lc rgb "#999900" lw 2, \
           "" using 7:(-1*$9) lc rgb "#00FF00" lw 2, \
           "" using 7:(-1*$9) lc rgb "#00FFFF" lw 2, \
           "" using 7:(-1*$9) lc rgb "#0000FF" lw 2  

#--- 2nd plot
set origin 0,0
set size 1, 0.47
#set xlabel 'k [MeV/c]'
set ylabel 'P0n00'
#set xtics 0,20,300.
#set yrange [-1:1]
#set ytics -1,0.2,1

plot  0 lc rgb "#000000", \
     filename using ($7):($13) lc rgb "#FF0000" lw 2, \
           "" using ($7):($17) lc rgb "#FF6600" lw 2, \
           "" using ($7):($21) lc rgb "#999900" lw 2, \
           "" using ($7):($25) lc rgb "#00FF00" lw 2, \
           "" using ($7):($29) lc rgb "#00FFFF" lw 2, \
           "" using ($7):($33) lc rgb "#0000FF" lw 2

unset multiplot


#---------- pause and wait for eneter to stop

pause mouse keypress (" ==== press s to save, any for exit \n")
#print 'keypress ASCII = '.MOUSE_KEY

if (MOUSE_KEY == 115) \
set term png small size 900, 600;\
set output pngName;\
print 'png output -> '.pngName;\
else \
print 'no output';