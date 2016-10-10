/******************************************************************************
 * Parallell quick sort algorithm
 * By Alessandro Piccolo
 *****************************************************************************/
 
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void quick_sort_serial(double *vector, int a, int b);
void quick_sort_par(double *data, int length, MPI_Comm comm, int *last_length);
int search_split_point(double *data, int length, double pivot);
void merge(double *d1, int len1, double *d2, int len2, double *data_sort);
void copyArray(double *d1, double *d2, int start, int length);
void swap(double *vector, double temp, int i, int j);
double findMean(double *vector, int length);
int dbgSorted(double* data, int length);
void print_vector(double *vector, int length);
void print_vector_int(int* vector, int length);
 
int main(int argc, char *argv[]) {  
  int rank;                                   /* Id of processors            */
  int nprocs;                                 /* Number of processor         */
  int n = atoi(argv[1]);                      /* n vector long               */
  int root = 0, i;
  double pivot, time_init, time_end;
  double *data, *data_sub, *data_sorted;
  int *last_length, *recieve_counts, *receive_displacements;
  last_length = (int *)malloc(1*sizeof(int));
  
  MPI_Init(&argc, &argv);                    /* Initialize MPI               */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);    /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);      /* Get my number                */
  
  int elements_per_proc = n/nprocs;          /* Elements per processor       */

  /* Generate random numbers, pivot and start time */
  if (rank == 0) {
    if (n <= 1) {
      printf("Already sorted, only one element");
      return 1;    
    }
    data = (double *)malloc(n*sizeof(double));
    data_sorted = (double *)malloc(n*sizeof(double));
    recieve_counts = (int *)malloc(nprocs*sizeof(int));
    receive_displacements = (int *)malloc(nprocs*sizeof(int));

    srand(time(NULL));
    for (i = 0; i < n; i++) {
        data[i] = drand48();
    }
    time_init = MPI_Wtime();  
  }

  /* Scatter chunk of data into smaller data_sub vectors */  
  data_sub = (double *)malloc(n*sizeof(double)); 
  MPI_Scatter(data, elements_per_proc, MPI_DOUBLE, data_sub, 
               elements_per_proc, MPI_DOUBLE, root, MPI_COMM_WORLD);
  
  /* Quick sort in serial followed by parallel quick sort algorithm */
  quick_sort_serial(data_sub, 0, elements_per_proc-1);
  quick_sort_par(data_sub, elements_per_proc, MPI_COMM_WORLD, last_length);
  elements_per_proc = *last_length;

  /* Number of elements in each processor */  
  MPI_Gather(&elements_per_proc, 1, MPI_INT, recieve_counts, 1, 
             MPI_INT, root, MPI_COMM_WORLD);
  
  /* Fill displacement vector for gatherv function */
  if (rank == 0) {
    int index = 0; receive_displacements[0] = index;
    for (i = 1; i < nprocs; i++) {
      index = index + recieve_counts[i-1];
      receive_displacements[i] = index;
    }
  }
 
 /* Gather all sub vectors into one large sorted vector  */
  MPI_Gatherv(data_sub, elements_per_proc, MPI_DOUBLE, data_sorted,  
      recieve_counts, receive_displacements, MPI_DOUBLE, root, MPI_COMM_WORLD);
  
  /* End time and free datatypes */
  if (rank == 0) {
    time_end = MPI_Wtime() - time_init;
    printf("\nElapsed time: %.16f s\n", time_end);

    if (!dbgSorted(data_sorted, n)) {
       printf("Error not sorted \n");
    }
  }

  free(last_length);
  free(data);
  free(data_sub);
  free(data_sorted);
  free(receive_displacements);
  free(recieve_counts);
  MPI_Finalize();                            /* Shut up and clean down MPI   */
  return 0;
}

/*
 * Sorts a vector with recursive quick sort algorithm, a and b defines
 * the interval of what to sort. Variable b is also the index of our 
 * pivot element.
 */
void quick_sort_serial(double *vector, int a, int b) {
  double pivot;
  int i, j;

  if(a < b) {
    pivot = vector[b];
    i = a;
    j = b;

    while(i < j) {
      while(vector[i] < pivot) {
        i++;
      }
      while(vector[j] >= pivot) {
        j--;
      }
      if(i < j) {
        swap(vector, vector[i], i, j);
      }
    }
    swap(vector, pivot, b, i);
    
    quick_sort_serial(vector, a, i-1);
    quick_sort_serial(vector, i+1, b);
  }
}

/* 
 * Quick sort parallel, length is the number of elements in 
 * ata last_length is pointer to the number of elements in 
 * the last recursive functions data vector.
 */
void quick_sort_par(double *data, int length, MPI_Comm comm, int *last_length){
  int nprocs, rank, pivot_index, number_amount, length_keep, i;
  double pivot, mean_local[2] = {0.0,length}, mean_global[2]={0.0,0.0}; 
  double *data_recieve, *data_keep;
  MPI_Status status;
  MPI_Comm new_comm;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
    
  if (nprocs == 1) { /* We are done, do not have to proceed */
    *last_length = length;
    return;
  }

  /* Calculate pivot as mean value across all processors */
  for (i = 0; i < length; i++) { /* Sum the numbers locally */
    mean_local[0] = mean_local[0] + data[i];
  }
  MPI_Reduce(&mean_local, &mean_global, 2, MPI_DOUBLE, MPI_SUM, 0, comm);
  if (rank == 0) {
    pivot = mean_global[0]/mean_global[1] ; /* Pivot is meanvalue */
  }
  MPI_Bcast(&pivot, 1, MPI_DOUBLE, 0, comm);

  /* Find local index where to split the vector into lesser and greater half */
  pivot_index = search_split_point(data, length, pivot);

  /* Exchange halves */
  if (rank < nprocs/2) {
    /* Send elem higher than pivot from left hand side processor */
    MPI_Send(data + (pivot_index + 1), (length - 1) - pivot_index, MPI_DOUBLE,
             rank + nprocs/2, 00, comm);
    /* Recieve elem lower than pivot from right hand side processor */
    MPI_Probe(rank + nprocs/2, 11, comm, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &number_amount);
    data_recieve = (double *)malloc(number_amount*sizeof(double));
    MPI_Recv(data_recieve, number_amount, MPI_DOUBLE, rank + nprocs/2, 
             11, comm, MPI_STATUS_IGNORE);
  }
  else {
    /* Recieve elem higher than pivot from left hand side processor */
    MPI_Probe(rank - nprocs/2, 00, comm, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &number_amount);
    data_recieve = (double *)malloc(number_amount*sizeof(double));
    MPI_Recv(data_recieve, number_amount, MPI_DOUBLE,
             rank - nprocs/2, 00, comm, MPI_STATUS_IGNORE);
    /* Send elem lower than pivot from right hand side processor */
    MPI_Send(data, pivot_index + 1, MPI_DOUBLE, rank - nprocs/2, 11, comm); 
  }

  /* Create new array data to be kept, and to be merged with data recieved */
  if (rank < nprocs/2) {
    length_keep = pivot_index + 1;
    data_keep = (double *)malloc(length_keep*sizeof(double));
    copyArray(data, data_keep, 0, length_keep);
  }
  else {
    length_keep = (length - 1) - pivot_index;
    data_keep = (double *)malloc(length_keep*sizeof(double));
    copyArray(data, data_keep, pivot_index+1, length_keep);
  }

  /* Merge data into one */
  merge(data_recieve, number_amount, data_keep, length_keep, data);

  /* Split communicator into two parts, left and right */
  int color = rank/(nprocs/2);
  MPI_Comm_split(comm, color, rank, &new_comm);

  quick_sort_par(data, length_keep + number_amount, new_comm, last_length);
}

/* 
 * Search for split point index. Returns the pivot index, [0,pivot_index]
 * are considered as the lower half <= pivot value.
 */
int search_split_point(double *data, int length, double pivot) {
  int i = 0;
  for (i = 0; i < length; i++) {
    if (data[i] > pivot) {
      return i-1;
    }
  } 
  return i-1; /* Zero entry returns -1 */
}

/*
 * Merge two sorted arrays into one Sort d1 and d2 into data_sort
 * d1 and d2 have to be already sorted incrementally.
 */ 
void merge(double *d1, int len1, double *d2, int len2, double *data_sort) {
  int i = 0, j = 0, index = 0;
  
  while(i < len1 && j < len2) {
    if (d1[i] < d2[j]) {
      data_sort[index] = d1[i]; 
      i++;
    }
    else {
      data_sort[index] = d2[j];
      j++;
    }
    index++;
  }  
  if (i >= len1) {
    while (j < len2) {
      data_sort[index] = d2[j];
      j++;
      index++;
    }
  }
  if (j >= len2) {
    while (i < len1) {
      data_sort[index] = d1[i];
      i++;
      index++;
    }
  }
}

/* 
 * Copy d1 elements from index start with length steps into vector d2
 */ 
void copyArray(double *d1, double *d2, int start, int length) {
  int i, j = start;
  for (i = 0; i < length; i++) {
    d2[i] = d1[j];
    j++;
  }
}

/*
 * Swaps elem in place i with j in vector
 */
void swap(double *vector, double temp, int i, int j) {
  vector[i] = vector[j];
  vector[j] = temp;
}

/*
 * Checks if the data is sorted
 */
int dbgSorted(double* data, int length) {
  int i;
  for (i = 0; i < length-1; ++i) {
    if (data[i] > data[i+1]) {
      return 0;
    }
  }
  return 1;
}

/*
 * Prints vector double
 */
void print_vector(double* vector, int length) {
  int i;
  for (i = 0; i < length; i++) {
    printf("%f ", vector[i]);
  }
  printf("\n");
}

/*
 * Prints vector int
 */
void print_vector_int(int* vector, int length) {
  int i;
  for (i = 0; i < length; i++) {
    printf("%d ", vector[i]);
  }
  printf("\n");
}
