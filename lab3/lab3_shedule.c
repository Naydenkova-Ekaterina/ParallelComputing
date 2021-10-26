#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

#ifdef _OPENMP
#include "omp.h"
#else
void omp_set_num_threads(int M) { }
#endif

#define A 504

double * generate_array(double * array, unsigned int * seed, int size, int min, int max) {
  for (int i = 0; i < size; i++) {
    array[i] = min + ((double)(rand_r(seed) % (max - min)));
  }
  return array;
}

void map_array1(double * array1, int size) {
  int i;
  #pragma omp parallel for default(none) private(i) shared(array1, size) schedule(SCHEDULE_TYPE, CHUNCK_SIZE)
  for (i = 0; i < size; i++) {
    array1[i] = exp(sqrt(array1[i]));
  }
}

void copy_array2(double * array2, double * array2_copy, int size) {
  int i;
  #pragma omp parallel for default(none) private(i) shared(array2, array2_copy, size) schedule(SCHEDULE_TYPE, CHUNCK_SIZE)
  for (i = 0; i < size; i++) {
    array2_copy[i+1] = array2[i];
  }
}

void map_array2(double * array2, double * array2_copy, int size) {
  int i;

  #pragma omp parallel for default(none) private(i) shared(array2, array2_copy, size) schedule(SCHEDULE_TYPE, CHUNCK_SIZE)
  for (i = 0; i < size; i++) {

    if (i > 0) {
      array2[i] = array2[i] + array2_copy[i];
    } 
    array2[i] = fabs(tan(array2[i]));
  }
}

void merge(double * array1, double * array2, int size) {
  int i;
  #pragma omp parallel for default(none) private(i) shared(array1, array2, size) schedule(SCHEDULE_TYPE, CHUNCK_SIZE)
  for (i = 0; i < size; i++) {
    if (array1[i] < array2[i]) {
      array2[i] = array1[i];
    }
  }
}

void swap(double * x, double * y) {
  int temp = * x;
  * x = * y;
  * y = temp;
}

void selection_sort(double * array, int size) {
  int i, j, min_idx;

  // One by one move boundary of unsorted subarray
  for (i = 0; i < size - 1; i++) {
    // Find the minimum element in unsorted array
    min_idx = i;
    for (j = i + 1; j < size; j++)
      if (array[j] < array[min_idx])
        min_idx = j;

    // Swap the found minimum element with the first element
    swap( & array[min_idx], & array[i]);
  }
}

double reduce(double * array, int size) {
  double res = 0;
  double min = array[0];
  for (int i = 1; i < size; i++) {
    if (array[i] < min && array[i] != 0) {
      min = array[i];
    }
  }
  int i;
  #pragma omp parallel for default(none) private(i) shared(size, min, array)  reduction(+:res) schedule(SCHEDULE_TYPE, CHUNCK_SIZE)
  for (i = 0; i < size; i++) {
    if ((int)(array[i] / min) % 2 == 0) {
      res += sin(array[i]);
    }
  }

  return res;
}

void printArray(double * array, int size) {
  for (int i = 0; i < size; i++) {
    printf("%f ", array[i]);
  }
  printf("\n");
}

int main(int argc, char * argv[]) {
  int i, N;
  unsigned int seed;
  struct timeval T1, T2;
  long delta_ms;

  if (argc == 3) {
    N = atoi(argv[1]); /* N равен первому параметру командной строки */
    omp_set_num_threads(atoi(argv[2]));
  } else {
    printf("Usage: ./lab1 N M");
    return 0;
  }

  gettimeofday( & T1, NULL); /* запомнить текущее время T1 */

  double * first_array = malloc(sizeof(double) * N);
  double * second_array = malloc(sizeof(double) * (N / 2));
  double * second_array_copy = malloc(sizeof(double) * (N / 2) + 1);
  second_array_copy[0] = 0;

  double X;

  for (i = 0; i < 100; i++) /* 100 экспериментов */ {

    seed = i;
    /* Заполнить массив исходных данных размером N */
    generate_array(first_array, & seed, N, 1, A);
    generate_array(second_array, & seed, N / 2, A, 10 * A);
    //printArray(first_array, N);
    //printArray(second_array, N/2);

    /* Решить поставленную задачу, заполнить массив с результатами */
    map_array1(first_array, N);
    copy_array2(second_array, second_array_copy, N / 2);
    map_array2(second_array, second_array_copy, N / 2);

    merge(first_array, second_array, N / 2);

    /* Отсортировать массив с результатами указанным методом */
    selection_sort(second_array, N / 2);

    //printArray(second_array, N/2);
    X = reduce(second_array, N / 2);
    printf("\n%d: X=%f \n", i, X);
  }

  free(first_array);
  free(second_array);
  free(second_array_copy);

  gettimeofday( & T2, NULL); /* запомнить текущее время T2 */
  delta_ms = 1000 * (T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;

  printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms); /* T2 - T1 */

  return 0;
}
