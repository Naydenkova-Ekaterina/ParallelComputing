#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>

#if defined(_OPENMP)
#include "omp.h"
#else
struct timeval timeval;
double omp_get_wtime()
{
  gettimeofday(&timeval, NULL);
  return (double)timeval.tv_sec + (double)timeval.tv_usec / 1000000;
}

void omp_set_nested(int val) { }

void omp_set_num_threads(int M) { }
#endif

#define A 504

double * generate_array(double * array, unsigned int * seed, int size, int min, int max) {
  #pragma omp parallel for default(none) shared(array, seed, size, min, max)
  for (int i = 0; i < size; i++) {
    unsigned int seed_i = i + *seed;
    array[i] = min + ((double)(rand_r(&seed_i) % (max - min)));
  }
  return array;
}

void map_array1(double * array1, int size) {
  #pragma omp parallel for default(none) shared(array1, size)
  for (int i = 0; i < size; i++) {
    array1[i] = exp(sqrt(array1[i]));
  }
}

void copy_array2(double * array2, double * array2_copy, int size) {
  #pragma omp parallel for default(none) shared(array2, array2_copy, size)
  for (int i = 0; i < size; i++) {
    array2_copy[i+1] = array2[i];
  }
}

void map_array2(double * array2, double * array2_copy, int size) {
  #pragma omp parallel for default(none) shared(array2, array2_copy, size)
  for (int i = 0; i < size; i++) {

    if (i > 0) {
      array2[i] = array2[i] + array2_copy[i];
    } 
    array2[i] = fabs(tan(array2[i]));
  }
}

void merge(double * array1, double * array2, int size) {
  #pragma omp parallel for default(none) shared(array1, array2, size)
  for (int i = 0; i < size; i++) {
    if (array1[i] < array2[i]) {
      array2[i] = array1[i];
    }
  }
}

void swap(double * x, double * y) {
  double temp = * x;
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

void mergeArrays(double *dest, double *arr1, double *arr2, int len1, int len2)
{
  int i, j, k;
  i = j = k = 0;
  for (i = 0; i < len1 && j < len2;) {
    if (arr1[i] < arr2[j]) {
      dest[k] = arr1[i];
      k++;
      i++;
    } else {
      dest[k] = arr2[j];
      k++;
      j++;
    }
  }
  while (i < len1) {
    dest[k] = arr1[i];
    k++;
    i++;
  }
  while (j < len2) {
    dest[k] = arr2[j];
    k++;
    j++;
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

  #pragma omp parallel for default(none) shared(size, min, array)  reduction(+:res)
  for (int i = 0; i < size; i++) {
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

int main_work(int argc, char* argv[], int *progress) {
  int i, N;
  unsigned int seed;
  double T1, T2;
  long delta_ms;

  if (argc == 3) {
    N = atoi(argv[1]); /* N равен первому параметру командной строки */
    omp_set_num_threads(atoi(argv[2]));
  } else {
    printf("Usage: ./lab1 N M");
    return 0;
  }

  T1 = omp_get_wtime(); /* запомнить текущее время T1 */

  double * first_array = malloc(sizeof(double) * N);
  double * second_array = malloc(sizeof(double) * (N / 2));
  #if defined(_OPENMP)
  double * merged_array = malloc(sizeof(double) * (N / 2));
  #endif
  double * second_array_copy = malloc(sizeof(double) * (N / 2) + 1);
  second_array_copy[0] = 0;

  double X;
  int num_of_iterations = 100;
  for (i = 0; i < num_of_iterations; i++) /* 100 экспериментов */ {
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

    #if defined(_OPENMP)
    #pragma omp parallel sections
    {
        #pragma omp section
        selection_sort(second_array, N / 4);
        #pragma omp section
        selection_sort(second_array + (N / 4), N / 2 - N / 4);
    }
    mergeArrays(merged_array, second_array, second_array + (N / 4), N / 4, N / 2 - N / 4);
    #else
    selection_sort(second_array, N / 2);
    #endif

    #if defined(_OPENMP)
    // printArray(merged_array, N/2);
     X = reduce(merged_array, N / 2);
    #else
    // printArray(second_array, N/2);
     X = reduce(second_array, N / 2);
    #endif
     printf("\n%d: X=%f \n", i, X);

    #if defined(_OPENMP)
    #pragma omp critical
    {
      *progress = 100 * (i + 1) / num_of_iterations;
    }
    #endif
  }

  free(first_array);
  free(second_array);
  free(second_array_copy);
  #if defined(_OPENMP)
  free(merged_array);
  #endif

  T2 = omp_get_wtime(); /* запомнить текущее время T2 */
  delta_ms = 1000 * (T2 - T1);

  printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms); /* T2 - T1 */
  return 0;
}

#if defined(_OPENMP)
void time_progress(int *progress)
{
  int progress_value;

  while(progress_value < 100) {
    #pragma omp critical
    {
    progress_value = *progress;
    }
    printf("Progress = %d%%\n", progress_value);
    sleep(1);
  }
}
#endif

int main(int argc, char * argv[]) {

  int progress = 0;

  omp_set_nested(1);

  #if defined(_OPENMP)
  #pragma omp parallel sections shared(progress)
  {
    #pragma omp section
    time_progress(&progress);
    #pragma omp section
    main_work(argc, argv, &progress);
  }
  #else
    main_work(argc, argv, &progress);
  #endif

  return 0;
}
