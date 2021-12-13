#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>

#if defined(_OPENMP)
#include "omp.h"
#else
struct timeval timeval;
double omp_get_wtime()
{
  gettimeofday(&timeval, NULL);
  return (double)timeval.tv_sec + (double)timeval.tv_usec / 1000000;
}
#endif

#define A 504

struct thread_params {
  double *first_array;
  double *second_array;
  int thread_id;
  int size;
  int min;
  int max;
  double x;
  int N;
  int num_of_threads;
  unsigned int *seed;
  int* progress;
  int chunk;
};

enum task_steps {
  GENERATE = 0,
  MAP1,
  COPY,
  MAP2,
  MERGE,
  SORT,
  REDUCE,
};

double steps_time[7] = { 0 };

void print_steps_time() {
  for (int i = 0; i < 7; i++) {
    printf(" %f", steps_time[i]);
  }
  printf("\n");
}

void execute_in_parallel(double * array1, double * array2, int size_param, void * func_param, int num_of_threads) { 	 
  void *retval; 									       	
  pthread_t threads[num_of_threads];							       
  struct thread_params *params = malloc(sizeof(struct thread_params) * num_of_threads);	      
											       
  for (int i = 0; i < num_of_threads; i++) {						       
    params[i].first_array = array1;							       
    params[i].second_array = array2;					       
    params[i].size = size_param;							       
    params[i].thread_id = i;								       
    params[i].num_of_threads = num_of_threads; 						       
    pthread_create(&threads[i], NULL, func_param, params + i);				       
  }											       
											       
  for (int i = 0; i < num_of_threads; i++) {						       
    pthread_join(threads[i], &retval);							       
  }											       
}

pthread_mutex_t progress_lock;
pthread_mutex_t shedule_dynamic_lock;

void *generate_array_function(void * params) {
  struct thread_params *thread_params = params;

  int size = thread_params->size;
  int min = thread_params->min;
  int max = thread_params->max;
  double *array = thread_params->first_array;
  unsigned int *seed = thread_params->seed;
  int chunk = thread_params->chunk;
  int i, j;
  int *shared_i = thread_params->progress;

  while(1) {
    pthread_mutex_lock(&shedule_dynamic_lock);
    i = *shared_i;
    (*shared_i) += chunk;
    pthread_mutex_unlock(&shedule_dynamic_lock);

    if (i >= size) break;

    for (j = i; j < (i + chunk) && j < size; j++, i++) {
      unsigned int seed_i = i + *seed;
      array[i] = min + ((double)(rand_r(&seed_i) % (max - min)));
    }
  }

  pthread_exit(NULL);
}

void generate_array(double * array, unsigned int * seed, int size, int min, int max, int num_of_threads) {
  double T1, T2;
  T1 = omp_get_wtime();

  int shared_i = 0;
  void *retval;
  pthread_t threads[num_of_threads];
  int chunk_size = size / num_of_threads;
  struct thread_params *params = malloc(sizeof(struct thread_params) * num_of_threads);

  for (int i = 0; i < num_of_threads; i++) {
    params[i].first_array = array;
    params[i].min = min;
    params[i].max = max;
    params[i].seed = seed;
    params[i].size = size;
    params[i].progress = &shared_i;
    params[i].chunk = chunk_size;
    pthread_create(&threads[i], NULL, generate_array_function, params + i);
  }

  for (int i = 0; i < num_of_threads; i++) {
    pthread_join(threads[i], &retval);
  }

  T2 = omp_get_wtime();
  steps_time[GENERATE] += (T2 - T1) * 1000;
}

void *map_array1_function(void *params) {
  struct thread_params *thread_params = params;
  int thread_id = thread_params->thread_id;
  int size = thread_params->size;
  double *array1 = thread_params->first_array;

  for (int i = thread_id; i < size; i += thread_params->num_of_threads) {
    array1[i] = exp(sqrt(array1[i]));
  }

  pthread_exit(NULL);
}

void map_array1(double * array1, int size, int num_of_threads) {
  double T1, T2;
  T1 = omp_get_wtime();
  execute_in_parallel(array1, NULL, size, map_array1_function, num_of_threads);
  T2 = omp_get_wtime();
  steps_time[MAP1] += (T2 - T1) * 1000;
}

void *copy_array2_function(void *params) {
  struct thread_params *thread_params = params;
	
  int thread_id = thread_params->thread_id;
  int size = thread_params->size;
  double *src = thread_params->first_array;
  double *dst = thread_params->second_array;

  for (int i = thread_id; i < size; i += thread_params->num_of_threads) {
    dst[i+1] = src[i];
  }

  pthread_exit(NULL);
}

void copy_array2(double * array2, double * array2_copy, int size, int num_of_threads) {
  double T1, T2;
  T1 = omp_get_wtime();
  execute_in_parallel(array2, array2_copy, size, copy_array2_function, num_of_threads);
  T2 = omp_get_wtime();
  steps_time[COPY] += (T2 - T1) * 1000;
}

void *map_array2_function(void *params) {
  struct thread_params *thread_params = params;
	
  int thread_id = thread_params->thread_id;
  int size = thread_params->size;
  double *array2 = thread_params->first_array;
  double *array2_copy = thread_params->second_array;

  for (int i = thread_id; i < size; i += thread_params->num_of_threads) {
    if (i > 0) {
      array2[i] = array2[i] + array2_copy[i];
    }
    array2[i] = fabs(tan(array2[i]));
  }

  pthread_exit(NULL);
}

void map_array2(double * array2, double * array2_copy, int size, int num_of_threads) {
  double T1, T2;
  T1 = omp_get_wtime();
  execute_in_parallel(array2, array2_copy, size, map_array2_function, num_of_threads);
  T2 = omp_get_wtime();
  steps_time[MAP2] += (T2 - T1) * 1000;
}

void *merge_function(void *params) {
  int i;
  struct thread_params *thread_params = params;
	
  int thread_id = thread_params->thread_id;
  int size = thread_params->size;
  double *array1 = thread_params->first_array;
  double *array2 = thread_params->second_array;

  for (i = thread_id; i < size; i += thread_params->num_of_threads) {
    if (array1[i] < array2[i]) {
      array2[i] = array1[i];
    }
  }

  pthread_exit(NULL);
}

void merge(double * array1, double * array2, int size, int num_of_threads) {
  double T1, T2;
  T1 = omp_get_wtime();
  execute_in_parallel(array1, array2, size, merge_function, num_of_threads);
  T2 = omp_get_wtime();
  steps_time[MERGE] += (T2 - T1) * 1000;
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

void *selection_sort_function(void *params)
{
  struct thread_params *thread_params = params;
  int size = thread_params->size;
  double *array = thread_params->first_array;

  selection_sort(array, size);
  pthread_exit(NULL);
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

void *reduce_function(void *params) {

  struct thread_params *thread_params = params;
  int thread_id = thread_params->thread_id;
  int size = thread_params->size;
  double *array = thread_params->first_array;
  double min = thread_params->x;
  double res = 0;

  for (int i = thread_id; i < size; i += thread_params->num_of_threads) {
    if ((int)(array[i] / min) % 2 == 0) {
      res += sin(array[i]);
    }
  }

  thread_params->x = res;
  pthread_exit(NULL);
}

double reduce(double * array, int size, int num_of_threads) {

  double res = 0;
  double min = array[0];

  double T1, T2;
  T1 = omp_get_wtime();

  pthread_t threads[num_of_threads];
  struct thread_params *thread_params = malloc(sizeof(struct thread_params) * num_of_threads);
  void *retval;

  for (int i = 1; i < size; i++) {
    if (array[i] < min && array[i] != 0) {
      min = array[i];
    }
  }

  for (int i = 0; i < num_of_threads; i++) {
    thread_params[i].first_array = array;
    thread_params[i].x = min;
    thread_params[i].size = size;
    thread_params[i].thread_id = i;
    thread_params[i].num_of_threads = num_of_threads;
    pthread_create(&threads[i], NULL, reduce_function, thread_params + i);
  }

  for (int i = 0; i < num_of_threads; i++) {
    pthread_join(threads[i], &retval);
    res += thread_params[i].x;
  }

  T2 = omp_get_wtime();
  steps_time[REDUCE] += (T2 - T1) * 1000;

  return res;
}

void printArray(double * array, int size) {
  for (int i = 0; i < size; i++) {
    printf("%f ", array[i]);
  }
  printf("\n");
}

void *main_work(void *params) {
  int i, j, N;
  unsigned int seed;
  double T1, T2;
  long delta_ms;
  struct thread_params *thr_params = params;
  struct thread_params *params_for_sort = malloc(sizeof(struct thread_params) * 2);
  pthread_t threads[2];
  void *retval;

  N = thr_params->N; /* N равен первому параметру командной строки */
  int *progress = thr_params->progress;

  T1 = omp_get_wtime(); /* запомнить текущее время T1 */

  double * first_array = malloc(sizeof(double) * N);
  double * second_array = malloc(sizeof(double) * (N / 2));
  double * merged_array = malloc(sizeof(double) * (N / 2));
  double * second_array_copy = malloc(sizeof(double) * (N / 2) + 1);
  second_array_copy[0] = 0;

  double X;
  int num_of_threads = thr_params->num_of_threads;
  int num_of_iterations = 100;
  for (i = 0; i < num_of_iterations; i++) /* 100 экспериментов */ {
    seed = i;

    /* Заполнить массив исходных данных размером N */
    generate_array(first_array, & seed, N, 1, A, num_of_threads);
    generate_array(second_array, & seed, N / 2, A, 10 * A, num_of_threads);
    //printArray(first_array, N);
    //printArray(second_array, N/2);

    /* Решить поставленную задачу, заполнить массив с результатами */
    map_array1(first_array, N, num_of_threads);
    copy_array2(second_array, second_array_copy, N / 2, num_of_threads);
    map_array2(second_array, second_array_copy, N / 2, num_of_threads);

    merge(first_array, second_array, N / 2, num_of_threads);

    /* Отсортировать массив с результатами указанным методом */

    double T1, T2;
    T1 = omp_get_wtime();
    params_for_sort[0].first_array = second_array;
    params_for_sort[0].size = N / 4;
    pthread_create(&threads[0], NULL, selection_sort_function, params_for_sort);

    params_for_sort[1].first_array = second_array + (N / 4);
    params_for_sort[1].size = N / 2 - N / 4;
    pthread_create(&threads[1], NULL, selection_sort_function, params_for_sort + 1);

    for (j = 0; j < 2; j++) {
      pthread_join(threads[j], &retval);
    }

    mergeArrays(merged_array, second_array, second_array + (N / 4), N / 4, N / 2 - N / 4);
    T2 = omp_get_wtime();
    steps_time[SORT] += (T2 - T1) * 1000;

    X = reduce(merged_array, N / 2, num_of_threads);
    printf("\n%d: X=%f \n", i, X);


    pthread_mutex_lock(&progress_lock);
    *progress = 100 * (i + 1) / num_of_iterations;
    pthread_mutex_unlock(&progress_lock);

  }

  free(first_array);
  free(second_array);
  free(second_array_copy);
  free(merged_array);

  T2 = omp_get_wtime(); /* запомнить текущее время T2 */
  delta_ms = 1000 * (T2 - T1);

  printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms); /* T2 - T1 */
  print_steps_time();
  pthread_exit(NULL);
}

void *time_progress(void *progress) {
  int progress_value;
	
  pthread_mutex_lock(&progress_lock);
  progress_value = *((int *) progress);
  pthread_mutex_unlock(&progress_lock);

  while(progress_value < 100) {
    pthread_mutex_lock(&progress_lock);
    progress_value = *((int *) progress);
    pthread_mutex_unlock(&progress_lock);
    printf("Progress = %d%%\n", progress_value);
    sleep(1);
  }

  pthread_exit(NULL);
}

int main(int argc, char * argv[]) {

  int progress = 0;
  void *retval;
  pthread_t t1, t2;

  struct thread_params *params = malloc(sizeof(struct thread_params));
  params->N = atoi(argv[1]);
  params->num_of_threads = atoi(argv[2]);
  params->progress = &progress;

  pthread_create(&t1, NULL, time_progress, &progress);
  pthread_create(&t2, NULL, main_work, params);

  pthread_join(t1, &retval);
  pthread_join(t2, &retval);

  return 0;
}
