stitch: stitch.c stitch.h
	gcc -g -Wall -std=gnu99 -fopenmp -o stitch stitch.c -lz
