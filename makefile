CC=g++
#CFLAGS = -Wno-write-strings

all: 3DeeGen_Tc_angc_angd.o 3DeeGen_Tc_angc_angd_inv.o 3DeeGen_Tc_angc_angd_betad_inv.o

3DeeGen_Tc_angc_angd.o: 3DeeGen_Tc_angc_angd.cpp knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_Tc_angc_angd.cpp -o 3DeeGen_Tc_angc_angd.o

3DeeGen_Tc_angc_angd_inv.o: 3DeeGen_Tc_angc_angd_inv.cpp knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_Tc_angc_angd_inv.cpp -o 3DeeGen_Tc_angc_angd_inv.o

3DeeGen_Tc_angc_angd_betad_inv.o: 3DeeGen_Tc_angc_angd_betad_inv.cpp knockout3D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_Tc_angc_angd_betad_inv.cpp -o 3DeeGen_Tc_angc_angd_betad_inv.o

clean:
	rm -rfv *.o 
