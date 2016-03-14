#all:
#	g++ -I ./ -std=c++0x -O4 library.h constant.h global.h structure.h functions.h  main.cpp adv_edit_distance_fast.cpp align_allreads.cpp alignment.cpp print_alignment.cpp sam_format.cpp utility.cpp -o advtestrun
CC=g++ 
LDFLAGS=-lm
CFLAGS=-I ./ -std=c++0x -c -O4

all: details

details: main.o align_allreads.o alignment.o seed_selection.o edit_distance.o optimize_kband.o  sam_format.o utility.o
	$(CC) $(LDFLAGS) main.o align_allreads.o alignment.o seed_selection.o edit_distance.o optimize_kband.o sam_format.o utility.o -o nanoblaster

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

align_allreads.o: align_allreads.cpp
	$(CC) $(CFLAGS) align_allreads.cpp

alignment.o: alignment.cpp
	$(CC) $(CFLAGS) alignment.cpp

seed_selection.o: seed_selection.cpp
	$(CC) $(CFLAGS) seed_selection.cpp

edit_distance.o: edit_distance.cpp
	$(CC) $(CFLAGS) edit_distance.cpp

optimize_kband.o: optimize_kband.cpp
	$(CC) $(CFLAGS) optimize_kband.cpp

sam_format.o: sam_format.cpp
	$(CC) $(CFLAGS) sam_format.cpp

utility.o: utility.cpp
	$(CC) $(CFLAGS) utility.cpp

clean:
	rm *.o nanoblaster

