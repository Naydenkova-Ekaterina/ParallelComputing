noomp:
	gcc -O3 -Wall -Werror -o lab3-noomp -Wno-unknown-pragmas lab3.c -lm  
  
omp: 
	gcc -O3 -Wall -Werror -o lab3-omp -fopenmp lab3.c -lm

sh-static-1: 
	gcc -O3 -Wall -Werror -o lab3-sh-static-1 -fopenmp -DSCHEDULE_TYPE=static -DCHUNCK_SIZE=1 lab3.c -lm

sh-static-2: 
	gcc -O3 -Wall -Werror -o lab3-sh-static-2 -fopenmp -DSCHEDULE_TYPE=static -DCHUNCK_SIZE=2 lab3.c -lm

sh-static-4: 
	gcc -O3 -Wall -Werror -o lab3-sh-static-4 -fopenmp -DSCHEDULE_TYPE=static -DCHUNCK_SIZE=4 lab3.c -lm

sh-static-10: 
	gcc -O3 -Wall -Werror -o lab3-sh-static-10 -fopenmp -DSCHEDULE_TYPE=static -DCHUNCK_SIZE=10 lab3.c -lm


sh-dynamic-1: 
	gcc -O3 -Wall -Werror -o lab3-sh-dynamic-1 -fopenmp -DSCHEDULE_TYPE=dynamic -DCHUNCK_SIZE=1 lab3.c -lm

sh-dynamic-2: 
	gcc -O3 -Wall -Werror -o lab3-sh-dynamic-2 -fopenmp -DSCHEDULE_TYPE=dynamic -DCHUNCK_SIZE=2 lab3.c -lm

sh-dynamic-4: 
	gcc -O3 -Wall -Werror -o lab3-sh-dynamic-4 -fopenmp -DSCHEDULE_TYPE=dynamic -DCHUNCK_SIZE=4 lab3.c -lm

sh-dynamic-10: 
	gcc -O3 -Wall -Werror -o lab3-sh-dynamic-10 -fopenmp -DSCHEDULE_TYPE=dynamic -DCHUNCK_SIZE=10 lab3.c -lm


sh-guided-1: 
	gcc -O3 -Wall -Werror -o lab3-sh-guided-1 -fopenmp -DSCHEDULE_TYPE=guided -DCHUNCK_SIZE=1 lab3.c -lm

sh-guided-2: 
	gcc -O3 -Wall -Werror -o lab3-sh-guided-2 -fopenmp -DSCHEDULE_TYPE=guided -DCHUNCK_SIZE=2 lab3.c -lm

sh-guided-4: 
	gcc -O3 -Wall -Werror -o lab3-sh-guided-4 -fopenmp -DSCHEDULE_TYPE=guided -DCHUNCK_SIZE=4 lab3.c -lm

sh-guided-10: 
	gcc -O3 -Wall -Werror -o lab3-sh-guided-10 -fopenmp -DSCHEDULE_TYPE=guided -DCHUNCK_SIZE=10 lab3.c -lm

clean:
	rm -f lab3-*
