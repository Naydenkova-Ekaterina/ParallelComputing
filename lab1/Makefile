all:
	gcc -O3 -Wall -Werror -o lab1-seq lab1.c -lm
	gcc -O3 -Wall -Werror -floop-parallelize-all -ftree-parallelize-loops=1 lab1.c -o lab1-par-1 -lm
	gcc -O3 -Wall -Werror -floop-parallelize-all -ftree-parallelize-loops=2 lab1.c -o lab1-par-2 -lm
	gcc -O3 -Wall -Werror -floop-parallelize-all -ftree-parallelize-loops=4 lab1.c -o lab1-par-4 -lm
	gcc -O3 -Wall -Werror -floop-parallelize-all -ftree-parallelize-loops=6 lab1.c -o lab1-par-6 -lm

cc:
	cc -O3 -Wall -Werror -o lab1-cc-seq.out lab1.c -lm
	cc -O3 -Wall -Werror -floop-parallelize-all -ftree-parallelize-loops=1 lab1.c -o lab1-cc-par-1.out -lm
	cc -O3 -Wall -Werror -floop-parallelize-all -ftree-parallelize-loops=2 lab1.c -o lab1-cc-par-2.out -lm
	cc -O3 -Wall -Werror -floop-parallelize-all -ftree-parallelize-loops=4 lab1.c -o lab1-cc-par-4.out -lm
	cc -O3 -Wall -Werror -floop-parallelize-all -ftree-parallelize-loops=6 lab1.c -o lab1-cc-par-6.out -lm

icc:
	icc -O3 -Wall -Werror -o lab1-icc-seq lab1.c -lm
	icc -parallel -par-threshold0 -par-num-threads=1 -o lab1-icc-par-1 lab1.c
	icc -parallel -par-threshold0 -par-num-threads=2 -o lab1-icc-par-2 lab1.c
	icc -parallel -par-threshold0 -par-num-threads=4 -o lab1-icc-par-4 lab1.c
	icc -parallel -par-threshold0 -par-num-threads=6 -o lab1-icc-par-6 lab1.c

clean:
	rm -f lab1-*
