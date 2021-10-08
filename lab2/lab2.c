#include <stdio.h>

#include <stdlib.h>

#include <sys/time.h>

#include <math.h>

#include "/home/ekaterina/Загрузки/FW_1.3.1_Lin64/fwBase.h"

#include "/home/ekaterina/Загрузки/FW_1.3.1_Lin64/fwSignal.h"

#define A 504

void printArray(double * array, int size) {
  for (int i = 0; i < size; i++) {
    printf("%f ", array[i]);
  }
  printf("\n");
}

double * generate(double * array, unsigned int * seed, int size, int min, int max) {
  for (int i = 0; i < size; i++) {
    array[i] = min + ((double)(rand_r(seed) % (max - min)));
  }
  return array;
}

void map_array1(double * array1, int size) {
  fwsSqrt_64f(array1, array1, size);
  fwsExp_64f(array1, array1, size);
}

void map_array2(double * array2, int size) {
  double array2_copy[(int) size + 1];
  array2_copy[0] = 0;
  fwsCopy_64f(array2, & array2_copy[1], size);
  fwsAdd_64f(array2, array2_copy, array2, size);
  fwsTan_64f_A50(array2, array2, size);
  fwsAbs_64f(array2, array2, size);
}

void merge(double * array1, double * array2, int size) {
  fwsMinEvery_32f_I((float * ) array1, (float * ) array2, size);
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
  for (int i = 0; i < size; i++) {
    if (array[i] < min && array[i] != 0) {
      min = array[i];
    }
  }

  for (int i = 0; i < size; i++) {
    if ((int)(array[i] / min) % 2 == 0) {
      res += sin(array[i]);
    }
  }

  return res;
}

int main(int argc, char * argv[]) {
  int i, N;
  unsigned int seed;
  struct timeval T1, T2;
  long delta_ms;

  if (argc == 3) {
    N = atoi(argv[1]); /* N равен первому параметру командной строки */
    fwSetNumThreads(atoi(argv[2]));

  } else {
    printf("Usage: ./lab1 N M");
    return 0;
  }

  gettimeofday( & T1, NULL); /* запомнить текущее время T1 */

  double * first_array = malloc(sizeof(double) * N);
  double * second_array = malloc(sizeof(double) * (N / 2));

  double X;

  for (i = 0; i < 5; i++) /* 100 экспериментов */ {
    seed = i;
    /* Заполнить массив исходных данных размером N */
    generate(first_array, & seed, N, 1, A);
    generate(second_array, & seed, N / 2, A, 10 * A);
    //printArray(first_array, N);
    //printArray(second_array, N/2);

    /* Решить поставленную задачу, заполнить массив с результатами */
    map_array1(first_array, N);
    map_array2(second_array, N / 2);

    merge(first_array, second_array, N / 2);

    /* Отсортировать массив с результатами указанным методом */
    selection_sort(second_array, N / 2);

    //printArray(second_array, N/2);
    X = reduce(second_array, N / 2);
    printf("\nX=%f \n", X);
  }

  gettimeofday( & T2, NULL); /* запомнить текущее время T2 */
  delta_ms = 1000 * (T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;

  printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms); /* T2 - T1 */

  return 0;
}
