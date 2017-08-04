stitch: stitch.c stitch.h
	gcc -g -Wall -std=gnu99 -fopenmp -O2 -o stitch stitch.c -lz
