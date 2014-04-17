## the bash command is
# $ gnuplot -e "variable=value" -p plot_kp.gp

##### command line -e input ##################

#if (!exists("pngIO")) pngIO = 0 ; 
#print 'pngIO = '.pngIO

##### SETTING ################################
filename = "paraOut_16O_JA1.0_JB1.5_angk60_angNN70_Sp10.0.dat"
sizeX = 600
sizeY = 600
##### Plotting script ########################
reset #reset all graphic setting
print '============================================'
print 'file : ',filename
pngName = filename[9:strlen(filename)-4].'.png'

set term x11 size sizeX,sizeY
set style data lines
#set key autotitle columnhead
set xrange [0:300]
set yrange [0:180]
set xlabel "k [MeV/c]"
set ylabel "angk [deg]"
set dgrid3d 300/5,180/5

splot filename using 1:2:11 title "DWIA 1s1/2" pal


#---------- pause and wait for eneter to stop

pause mouse keypress (" ==== press s to save, any for exit \n")
#print 'keypress ASCII = '.MOUSE_KEY

if (MOUSE_KEY == 115){
   set term png small size sizeX, sizeY
   set output pngName
   print 'png output -> '.pngName
   
}
if (MOUSE_KEY != 115){
   print 'no output'
}