## the bash command is
# $ gnuplot -e "variable=value" -p plot_kp.gp

##### command line -e input ##################

#if (!exists("pngIO")) pngIO = 0 ; 
#print 'pngIO = '.pngIO

##### SETTING ################################
filename = "paraOut_16O_JA1.0_JB1.5_angNN60_Sp10.0.dat"
sizeX = 1200
sizeY = 900
##### Plotting script ########################
reset #reset all graphic setting
print '============================================'
print 'file : ',filename
pngName = filename[9:strlen(filename)-4].'.png'

set term x11 size sizeX,sizeY
set style data lines

#---------- plot script
load "plot_k_angk_multi.gp"

#---------- pause and wait for eneter to stop

pause mouse keypress (" ==== press s to save, any for exit \n")
#print 'keypress ASCII = '.MOUSE_KEY

if (MOUSE_KEY == 115){
   set term png small size sizeX, sizeY
   set output pngName
   print 'png output -> '.pngName
   load "plot_k_angk_multi.gp"   
}
if (MOUSE_KEY != 115){
   print 'no output'
}