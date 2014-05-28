CC=g++

all: make_infile.o read_outfile.o 3DeeGen_k.o 3DeeGen_angk_angNN.o 3DeeGen_k_angk.o 3DeeGen_k_angk_angNN.o

3DeeGen_k.o: 3DeeGen_k.cpp XsecTransform.h knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_k.cpp -o 3DeeGen_k.o

3DeeGen_angk_angNN.o: 3DeeGen_angk_angNN.cpp XsecTransform.h knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_angk_angNN.cpp -o 3DeeGen_angk_angNN.o

3DeeGen_k_angk.o: 3DeeGen_k_angk.cpp XsecTransform.h knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_k_angk.cpp -o 3DeeGen_k_angk.o

3DeeGen_k_angk_angNN.o: 3DeeGen_k_angk_angNN.cpp XsecTransform.h knockout2D.h 3DeeGenLibrary.h nuclei_mass.h
	$(CC) 3DeeGen_k_angk_angNN.cpp -o 3DeeGen_k_angk_angNN.o

make_infile.o: make_infile.cpp
	$(CC) make_infile.cpp -o make_infile.o

read_outfile.o: read_outfile.cpp
	$(CC) read_outfile.cpp -o read_outfile.o

clean:
	rm -rfv *.o 
