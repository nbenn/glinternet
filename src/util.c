#include "util.h"

#include <stdbool.h>
#include <stdlib.h>
#include <numa.h>
#include <numaif.h>
#include <sys/types.h>
#include <unistd.h>
#include <sched.h>

typedef struct allocator_data {
  size_t size;
  size_t offset;
  size_t alignment;
  int node;
} allocator_data;

/* call this function to start a nanosecond-resolution timer */
struct timespec timer_start(void){
    struct timespec start_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    return start_time;
}

/* call this function to end a timer, returning seconds elapsed as double */
double timer_end(struct timespec start_time){
    struct timespec end_time;
    double accum;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time);
    accum = ( end_time.tv_sec - start_time.tv_sec )
          + ( end_time.tv_nsec - start_time.tv_nsec )
            / (double)1000000000;
    return accum;
}

SEXP get_my_numa_node() {
  int cpu = sched_getcpu();
  int node = numa_node_of_cpu(cpu);

  SEXP result;
  PROTECT(result = allocVector(INTSXP, 1));
  INTEGER(result)[0] = node;

  UNPROTECT(1);
  return result;
}

SEXP get_cpu_node_usage() {

  SEXP max_num_cpus, max_num_nodes;

  PROTECT(max_num_cpus   = allocVector(INTSXP, 1));
  PROTECT(max_num_nodes  = allocVector(INTSXP, 1));

  INTEGER(max_num_cpus)[0]  = numa_num_configured_cpus();
  INTEGER(max_num_nodes)[0] = numa_max_node()+1;

  SEXP cpu_to_node, node_used;

  PROTECT(cpu_to_node = allocVector(INTSXP, INTEGER(max_num_cpus)[0]));
  PROTECT(node_used   = allocVector(INTSXP, INTEGER(max_num_nodes)[0]));

  int i,j;

  for (i = 0; i < INTEGER(max_num_cpus)[0]; ++i) {
    INTEGER(cpu_to_node)[i] = -1;
  }
  for (i = 0; i < INTEGER(max_num_nodes)[0]; ++i) {
    INTEGER(node_used)[i] = 0;
  }

  struct bitmask *cpubuf;

  cpubuf = numa_allocate_cpumask();
  if (numa_sched_getaffinity(getppid(), cpubuf) < 0) {
    numa_free_cpumask(cpubuf);
    Rf_error("could not run numa_sched_getaffinity in get_cpu_node_usage()");
  }
  for (i=0; i<cpubuf->size && i<INTEGER(max_num_cpus)[0]; i++) {
    INTEGER(cpu_to_node)[i] = numa_bitmask_isbitset(cpubuf, i) ? numa_node_of_cpu(i) : -1;
  }
  numa_free_cpumask(cpubuf);

  for (i = 0; i < INTEGER(max_num_nodes)[0]; ++i) {
    for (j = 0; j < INTEGER(max_num_cpus)[0]; ++j) {
      if (INTEGER(cpu_to_node)[j] == i) {
        INTEGER(node_used)[i]++;
      }
    }
  }

  SEXP num_used_cpus, num_used_nodes;

  PROTECT(num_used_cpus  = allocVector(INTSXP, 1));
  PROTECT(num_used_nodes = allocVector(INTSXP, 1));

  INTEGER(num_used_cpus)[0] = 0;
  INTEGER(num_used_nodes)[0] = 0;

  for (i = 0; i < INTEGER(max_num_cpus)[0]; ++i) {
    if(INTEGER(cpu_to_node)[i] >= 0) INTEGER(num_used_cpus)[0]++;
  }
  for (i = 0; i < INTEGER(max_num_nodes)[0]; ++i) {
    if(INTEGER(node_used)[i] != 0) INTEGER(num_used_nodes)[0]++;
  }

  SEXP result, result_names;
  char *names[6] = {
    "cpu_to_node",
    "node_used",
    "max_num_cpus",
    "max_num_nodes",
    "num_used_cpus",
    "num_used_nodes"
  };

  PROTECT(result = allocVector(VECSXP, 6));  
  PROTECT(result_names = allocVector(STRSXP, 6));

  SET_VECTOR_ELT(result, 0, cpu_to_node);
  SET_VECTOR_ELT(result, 1, node_used);
  SET_VECTOR_ELT(result, 2, max_num_cpus);
  SET_VECTOR_ELT(result, 3, max_num_nodes);
  SET_VECTOR_ELT(result, 4, num_used_cpus);
  SET_VECTOR_ELT(result, 5, num_used_nodes);

  for(i = 0; i < 6; i++){
    SET_STRING_ELT(result_names, i, mkChar(names[i])); 
  }
  setAttrib(result, R_NamesSymbol, result_names);

  UNPROTECT(8);
  return result;
}

/*void printmask(char *name, struct bitmask *mask) {
  int i;

  printf("%s: ", name);
  for (i = 0; i <= mask->size; i++)
    if (numa_bitmask_isbitset(mask, i))
      printf("%d ", i);
  putchar('\n');
}*/

void* custom_alloc(R_allocator_t *allocator, size_t size) {
  ((allocator_data*)allocator->data)->size = size;
  int node = ((allocator_data*)allocator->data)->node;
  size_t offset = ((allocator_data*)allocator->data)->offset;
  void* orig_addr = NULL;
  if (node == -1) {
    SEXP cpuInfo = get_cpu_node_usage();
    SEXP node_used, max_num_nodes;
    PROTECT(node_used      = VECTOR_ELT(cpuInfo, 1));
    PROTECT(max_num_nodes  = VECTOR_ELT(cpuInfo, 3));

    struct bitmask *nodemask = numa_allocate_nodemask();

    //Rprintf("get_cpu_node_usage: ");
    for (int i = 0; i < INTEGER(max_num_nodes)[0]; ++i) {
      //Rprintf("%d ", INTEGER(node_used)[i]);
      if (INTEGER(node_used)[i] > 0) {
        numa_bitmask_setbit(nodemask, i);
      }
    }
    //Rprintf("\n");
    //printmask("alloc interleaved among nodes", nodemask);

    orig_addr = numa_alloc_interleaved_subset(size+offset, nodemask);    

    numa_free_nodemask(nodemask);
    UNPROTECT(2);
  } 
  else if (node >= 0) {
    orig_addr = numa_alloc_onnode(size+offset, (size_t)node);
  } 
  else {
    Rf_error("impossible node specification in custom_alloc");
  }
  if (orig_addr == NULL) {
    Rf_error("could not allocate memory in custom_alloc");
  }
  void* shifted_addr = (void*)((uintptr_t)orig_addr + offset);
  //Rprintf("alloc on node %d: %p (%d)\n", node, orig_addr, size+offset);
  return (void*) shifted_addr;
}

void custom_free(R_allocator_t *allocator, void* shifted_addr) {
  size_t size = ((allocator_data*)allocator->data)->size;
  size_t offset = ((allocator_data*)allocator->data)->offset;
  void* orig_addr = (void*)((uintptr_t)shifted_addr-offset);
  //Rprintf("freeing: %p (%d)\n", orig_addr, size+offset);
  numa_free(orig_addr, size+offset);
}

SEXP alloc(int alignment, R_xlen_t length, int node) {
  allocator_data* custom_allocator_data = malloc(sizeof(allocator_data));
  custom_allocator_data->offset = 56;
  custom_allocator_data->alignment = alignment;
  custom_allocator_data->size = 0;
  custom_allocator_data->node = node;

  R_allocator_t* custom_allocator = malloc(sizeof(R_allocator_t));
  custom_allocator->mem_alloc = &custom_alloc;
  custom_allocator->mem_free = &custom_free;
  custom_allocator->res = NULL;
  custom_allocator->data = custom_allocator_data;

  SEXP result = PROTECT(allocVector3(REALSXP, length, custom_allocator));

  uintptr_t addr = (uintptr_t)REAL(result);
  int counter = 0;

  while (max_alignment(addr) < custom_allocator_data->alignment) {
    if (counter >= 10) {
      Rf_error("memory allocation failed: could not get suitably aligned memory");
    }
    uintptr_t new_addr = addr;
    Rprintf("old address: %p, max alig: %d\n--> repeating allocation.\n",
      addr, max_alignment(addr));
    while (max_alignment(new_addr) < custom_allocator_data->alignment) {
      new_addr += sizeof(double);
    }
    Rprintf("new address: %p, max alig: %d\n",
      new_addr, max_alignment(new_addr));
    custom_allocator_data->offset = custom_allocator_data->offset + (new_addr - addr);
    UNPROTECT(1);
    result = PROTECT(allocVector3(REALSXP, length, custom_allocator));
    addr = (uintptr_t)REAL(result);
    counter++;
  }

  //Rprintf("Successful alloc; address: %p, max alig: %d\n",
  //  addr, max_alignment(addr));

  UNPROTECT(1);
  return result;
}

static void finalize_singles_array_pointer(SEXP ext) {
  if (NULL == R_ExternalPtrAddr(ext))
    return;
  //Rprintf("finalizing singles array\n");
  float *ptr = (float *) R_ExternalPtrAddr(ext);
  double size = REAL(R_ExternalPtrTag(ext))[0];
  numa_free(ptr, (size_t)size);
  R_ClearExternalPtr(ext);
}

SEXP alloc_z_interleaved (SEXP a, SEXP b, SEXP x, SEXP info) {
  R_xlen_t n = xlength(x);
  R_xlen_t nrows = (R_xlen_t)asInteger(a);
  R_xlen_t ncols = (R_xlen_t)asInteger(b);
  int alignment = 64;
  SEXP result, result_names, node_used, max_num_nodes, num_used_nodes;

  if(nrows * ncols != n) {
    Rf_error("dimensions mismatch: %d rows times %d cols != %d elements", 
      nrows, ncols, n);
  }

  PROTECT(node_used      = VECTOR_ELT(info, 1));
  PROTECT(max_num_nodes  = VECTOR_ELT(info, 3));
  PROTECT(num_used_nodes = VECTOR_ELT(info, 5));

  /* the data will be allocated once per node in floats and once in an
     interleaved fashion among all nodes */
  PROTECT(result = allocVector(VECSXP, INTEGER(max_num_nodes)[0]+1));

  /* set up names for individual matrices: "node_numNode" and "interl" */
  char names[INTEGER(max_num_nodes)[0]+1][10];
  for (int i = 0; i < INTEGER(max_num_nodes)[0]; ++i) {
    int index_size = snprintf(NULL, 0, "%d", i) + 1;
    char index[index_size];
    char name[] = "node_";
    int name_size = sizeof(name)/sizeof(char);
    if ((index_size+name_size+1) > (sizeof(names[i])/sizeof(char))) {
      Rf_error("too many nodes in alloc_z: strbuf too small.");
    }
    snprintf(index, sizeof(index), "%d", i); 
    strncpy(names[i], name, name_size);
    strncat(names[i], index, index_size);
  }
  char name[] = "interl";
  int name_size = sizeof(name)/sizeof(char);
  strncpy(names[INTEGER(max_num_nodes)[0]], name, name_size);

  PROTECT(result_names = allocVector(STRSXP, INTEGER(max_num_nodes)[0]+1));
  for (int i = 0; i < INTEGER(max_num_nodes)[0]+1; ++i) {
    SET_STRING_ELT(result_names, i, mkChar(names[i])); 
  }
  setAttrib(result, R_NamesSymbol, result_names);

  SEXP interl_array;
  PROTECT(interl_array = alloc(alignment, n, -1));
  SET_VECTOR_ELT(result, INTEGER(max_num_nodes)[0], interl_array);
  UNPROTECT(1);

  /* make 2d array out of used matrices */
  SEXP dims2;
  PROTECT(dims2 = allocVector(INTSXP, 2));
  INTEGER(dims2)[0] = nrows;
  INTEGER(dims2)[1] = ncols;
  setAttrib(VECTOR_ELT(result, INTEGER(max_num_nodes)[0]), R_DimSymbol, dims2);

  UNPROTECT(6);
  return result;
}

SEXP alloc_z_float (SEXP z, SEXP info) {
  
  SEXP node_used, max_num_nodes;
  PROTECT(node_used      = VECTOR_ELT(info, 1));
  PROTECT(max_num_nodes  = VECTOR_ELT(info, 3));

  SEXP original;
  PROTECT(original = VECTOR_ELT(z, INTEGER(max_num_nodes)[0]));
  R_xlen_t n = xlength(original);

  /* allocate space for x matrix on each numa node which is used;
     allocate a dummy int vector of length 1 on unused nodes */
  for (int i=0; i<INTEGER(max_num_nodes)[0]; ++i) {
    SEXP tmp_ptr, tmp_tag;
    if (INTEGER(node_used)[i] > 0) {
      /* the size of the allocated memory is needed in numa_free
         -> use the ExternalPtrTag to store this info; could be larger
         than MAXINT therefore store as double */
      tmp_tag = PROTECT(allocVector(REALSXP, 1));
      REAL(tmp_tag)[0] = n*sizeof(float);
      Rprintf("allocating matrix (size %d) on node %d.\n",
        (size_t)REAL(tmp_tag)[0], i);
      float* array = (float*)numa_alloc_onnode(n*sizeof(float), (size_t)i);
      PROTECT(tmp_ptr = R_MakeExternalPtr(array, tmp_tag, R_NilValue));
      R_RegisterCFinalizerEx(tmp_ptr, finalize_singles_array_pointer, TRUE);
    }
    else {
      tmp_tag = PROTECT(allocVector(REALSXP, 1));
      REAL(tmp_tag)[0] = 1*sizeof(int);
      int* array = (int*)numa_alloc_onnode(sizeof(int), (size_t)i);
      PROTECT(tmp_ptr = R_MakeExternalPtr(array, tmp_tag, R_NilValue));
      R_RegisterCFinalizerEx(tmp_ptr, finalize_singles_array_pointer, TRUE);
    }
    if (tmp_ptr == NULL) Rf_error("could not allocate float array.");
    SET_VECTOR_ELT(z, i, tmp_ptr);
    UNPROTECT(2);
  }

  /* pupulate dummy vectors with -1, populate others with data from 
     interleaved double array */
  for (int i=0; i<INTEGER(max_num_nodes)[0]; ++i) {
    if (INTEGER(node_used)[i] <= 0) {
      *(int*)R_ExternalPtrAddr(VECTOR_ELT(z, i)) = -1;
    }
    else if (INTEGER(node_used)[i] > 0) {
      float *flt_ptr = R_ExternalPtrAddr(VECTOR_ELT(z, i));
      for (R_xlen_t j = 0; j < n; ++j) {
        flt_ptr[j]    = (float)REAL(original)[j];
      }
    }
    else Rf_error("cannot determine if node used or not.");
  }

  UNPROTECT(3);
  return z;
}

SEXP extract_col(SEXP src, SEXP i) {
  SEXP dims2;
  PROTECT(dims2 = getAttrib(src, R_DimSymbol));
  R_xlen_t nrows = (R_xlen_t)INTEGER(dims2)[0];
  R_xlen_t ncols = (R_xlen_t)INTEGER(dims2)[1];

  PROTECT(i = coerceVector(i, INTSXP));
  R_xlen_t col_idx = (R_xlen_t)INTEGER(i)[0];

  if(col_idx < 1 || col_idx > ncols) {
    Rf_error("in extract_col: 1 <= %d <= %d needed", 
      col_idx, ncols);
  }

  col_idx -= 1;

  SEXP res;
  PROTECT(res = allocVector(REALSXP, nrows));

  for (R_xlen_t i = 0; i < nrows; ++i) {
    REAL(res)[i] = REAL(src)[i+col_idx*nrows];
  }

  UNPROTECT(3);
  return res;
}

SEXP import_col(SEXP col, SEXP dest, SEXP i, SEXP info) {
  
  R_xlen_t nelem = xlength(col);

  PROTECT(i = coerceVector(i, INTSXP));
  R_xlen_t col_idx = (R_xlen_t)INTEGER(i)[0];
  col_idx -= 1;

  SEXP max_num_nodes, dims2;
  PROTECT(max_num_nodes  = VECTOR_ELT(info, 3));
  PROTECT(dims2 = getAttrib(VECTOR_ELT(dest, INTEGER(max_num_nodes)[0]), 
    R_DimSymbol));
  R_xlen_t nrows = (R_xlen_t)INTEGER(dims2)[0];
  R_xlen_t ncols = (R_xlen_t)INTEGER(dims2)[1];

  if(nrows != nelem) {
    Rf_error("in import_col: %d != %d", nrows, nelem);
  }

  if(col_idx < 0 || col_idx > (ncols-1)) {
    Rf_error("in import_col: 1 <= %d <= %d needed", 
      col_idx, ncols);
  }

  SEXP interl;
  PROTECT(interl = VECTOR_ELT(dest, INTEGER(max_num_nodes)[0]));
  for (R_xlen_t j = 0; j < nelem; ++j) {
    REAL(interl)[j+col_idx*nrows] = REAL(col)[j];
  }

  UNPROTECT(4);
  return dest;
}
