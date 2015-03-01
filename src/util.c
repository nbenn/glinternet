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
  size_t node;
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

void* custom_alloc(R_allocator_t *allocator, size_t size) {
  ((allocator_data*)allocator->data)->size = size;
  size_t node = ((allocator_data*)allocator->data)->node;
  size_t offset = ((allocator_data*)allocator->data)->offset;
  void* orig_addr = numa_alloc_onnode(size+offset, node);
  //Rprintf("original addr: %p\n", orig_addr);
  //Rprintf("max alig     : %d\n", max_alignment((uintptr_t)orig_addr));
  void* shifted_addr = (void*)((uintptr_t)orig_addr + offset);
  //Rprintf("shifted addr : %p\n", shifted_addr);
  //Rprintf("max alig     : %d\n", max_alignment((uintptr_t)shifted_addr));
  return (void*) shifted_addr;
}

void custom_free(R_allocator_t *allocator, void* shifted_addr) {
  size_t size = ((allocator_data*)allocator->data)->size;
  size_t offset = ((allocator_data*)allocator->data)->offset;
  //Rprintf("shifted addr : %p\n", shifted_addr);
  //Rprintf("max alig     : %d\n", max_alignment((uintptr_t)shifted_addr));
  void* orig_addr = (void*)((uintptr_t)shifted_addr-offset);
  //Rprintf("original addr: %p\n", orig_addr);
  //Rprintf("max alig     : %d\n", max_alignment((uintptr_t)orig_addr));
  numa_free(orig_addr, size+offset);
}

SEXP alloc(int alignment, R_xlen_t length, size_t node) {
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

SEXP alloc_z(SEXP a, SEXP b, SEXP x, SEXP info) {
  R_xlen_t n = xlength(x);
  R_xlen_t nrows = (R_xlen_t)asInteger(a);
  R_xlen_t ncols = (R_xlen_t)asInteger(b);
  int alignment = 64;
  SEXP result, result_names, node_used, max_num_nodes;

  if(nrows * ncols != n) {
    Rf_error("dimensions mismatch: %d rows times %d cols != %d elements", 
      nrows, ncols, n);
  }

  PROTECT(node_used      = VECTOR_ELT(info, 1));
  PROTECT(max_num_nodes  = VECTOR_ELT(info, 3));

  PROTECT(result = allocVector(VECSXP, INTEGER(max_num_nodes)[0]));

  int numa_node_x = -1;
  double size_mb = 1024*1024;
  get_mempolicy(&numa_node_x, NULL, 0, (void*)REAL(x), MPOL_F_NODE | MPOL_F_ADDR);
  long free_mem = 0;
  long size = numa_node_size(numa_node_x, &free_mem);

  if (free_mem/size_mb < (1.1*n*sizeof(double)/size_mb)) {
    Rprintf("node %d: %f [MB] free of %f [MB]\n", numa_node_x, free_mem/size_mb, size/size_mb);
    Rprintf("projected size of z is %f [MB]\n", n*sizeof(double)/size_mb);
    Rprintf("--> the input matrix has to be stored to disc");
    if (numa_node_x >= 0 && numa_node_x <= INTEGER(max_num_nodes)[0]) {
      INTEGER(node_used)[numa_node_x] = INTEGER(node_used)[numa_node_x]*(-1);
    }
    else {
      Rf_error("confusing node specification in alloc_z");
    }
  }

  int i, count = 0;
  for (i = 0; i < INTEGER(max_num_nodes)[0]; ++i) {
    if(INTEGER(node_used)[i] > 0) count++;
  }
  if(count == 0) {
    Rf_error("only node %d can be used, which does not have enough memory!\n", numa_node_x);
  }

  /* set up names for individual matrices: node_numNode */
  char names[INTEGER(max_num_nodes)[0]][10];
  for (i = 0; i < INTEGER(max_num_nodes)[0]; ++i) {
    int index_size = snprintf(NULL, 0, "%d", i) + 1;
    char index[index_size];
    char name[] = "node_";
    int name_size = sizeof(name)/sizeof(char);
    if ((index_size+name_size+1) > (sizeof(names[i])/sizeof(char))) {
      Rf_error("too many nodes in alloc_z");
    }
    snprintf(index, sizeof(index), "%d", i); 
    strncpy(names[i], name, name_size);
    strncat(names[i], index, index_size);
  }
  PROTECT(result_names = allocVector(STRSXP, INTEGER(max_num_nodes)[0]));
  for (i = 0; i < INTEGER(max_num_nodes)[0]; ++i) {
    SET_STRING_ELT(result_names, i, mkChar(names[i])); 
  }
  setAttrib(result, R_NamesSymbol, result_names);

  /* allocate space for x matrix on each numa node which is used;
     allocate a dummy int vector of length 1 on unused nodes */
  for (i=0; i<INTEGER(max_num_nodes)[0]; ++i) {
    SEXP tmp;
    if (INTEGER(node_used)[i] > 0) {
      //Rprintf("allocating matrix on node %d.\n", i);
      PROTECT(tmp = alloc(alignment, n, i));
    }
    else {
      PROTECT(tmp  = allocVector(INTSXP, 1));
    }
    SET_VECTOR_ELT(result, i, tmp);
    UNPROTECT(1);
  }

  /* pupulate dummy vectors with -1 */
  for (i=0; i<INTEGER(max_num_nodes)[0]; ++i) {
    if (INTEGER(node_used)[i] <= 0) {
      INTEGER(VECTOR_ELT(result, i))[0] = -1;
    }
  }

  /* make 2d array out of used matrices */
  for (i=0; i<INTEGER(max_num_nodes)[0]; ++i) {
    if (INTEGER(node_used)[i] > 0) {
      SEXP dims2;
      PROTECT(dims2 = allocVector(INTSXP, 2));
      INTEGER(dims2)[0] = nrows;
      INTEGER(dims2)[1] = ncols;
      setAttrib(VECTOR_ELT(result, i), R_DimSymbol, dims2);
      UNPROTECT(1);
    }
  }

  UNPROTECT(4);
  return result;
}

SEXP retry_alloc_z(SEXP z, SEXP info) {
  SEXP node_used, max_num_nodes;
  PROTECT(node_used      = VECTOR_ELT(info, 1));
  PROTECT(max_num_nodes  = VECTOR_ELT(info, 3));
  int affected_node = -1;
  int allocated_node = -1;

  for (int i=0; i<INTEGER(max_num_nodes)[0]; ++i) {
    if (INTEGER(node_used)[i] < 0) {
      affected_node = i;
    }
    if (INTEGER(node_used)[i] > 0) {
      allocated_node = i;
    }
  }

  if (affected_node == -1 || allocated_node == -1) {
    Rf_error("node to be allocated: %d\nnode to be used: %d\n",
       affected_node, allocated_node);
  }

  long free_mem = 0;
  long size = numa_node_size(affected_node, &free_mem);
  double size_mb = 1024*1024;

  SEXP tmp;
  SEXP original;
  PROTECT(original = VECTOR_ELT(z, allocated_node));
  
  R_xlen_t n = xlength(original);
  int alignment = 64;

  Rprintf("node %d: %f [MB] free of %f [MB]\n", affected_node, free_mem/size_mb, size/size_mb);
  Rprintf("projected size of z is %f [MB]\n", n*sizeof(double)/size_mb);
  Rprintf("allocating matrix of size %d on node %d.\n", n, affected_node);

  PROTECT(tmp = alloc(alignment, n, affected_node));

  for (size_t i = 0; i < n; ++i) {
    REAL(tmp)[i] = REAL(original)[i];
  }

  SEXP dims2;
  PROTECT(dims2 = getAttrib(original, R_DimSymbol));
  setAttrib(tmp, R_DimSymbol, dims2);

  SET_VECTOR_ELT(z, affected_node, tmp);

  INTEGER(node_used)[affected_node] = INTEGER(node_used)[affected_node]*(-1);

  UNPROTECT(5);
  return z;
}

SEXP extract_col(SEXP src, SEXP i) {
  SEXP dims2;
  PROTECT(dims2 = getAttrib(src, R_DimSymbol));
  R_xlen_t nrows = (R_xlen_t)INTEGER(dims2)[0];
  R_xlen_t ncols = (R_xlen_t)INTEGER(dims2)[1];

  PROTECT(i = coerceVector(i, INTSXP));
  int col_idx = *INTEGER(i);

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
  
  SEXP dims2;
  R_xlen_t nelem = xlength(col);

  PROTECT(i = coerceVector(i, INTSXP));
  int col_idx = *INTEGER(i);
  col_idx -= 1;

  SEXP node_used, max_num_nodes;
  PROTECT(node_used      = VECTOR_ELT(info, 1));
  PROTECT(max_num_nodes  = VECTOR_ELT(info, 3));

  for (int i=0; i<INTEGER(max_num_nodes)[0]; ++i) {
    if (INTEGER(node_used)[i] > 0) {
      for (R_xlen_t j = 0; j < nelem; ++j) {
        PROTECT(dims2 = getAttrib(VECTOR_ELT(dest, i), R_DimSymbol));
        R_xlen_t nrows = (R_xlen_t)INTEGER(dims2)[0];
        R_xlen_t ncols = (R_xlen_t)INTEGER(dims2)[1];

        if(nrows != nelem) {
          Rf_error("in import_col: %d != %d", nrows, nelem);
        }
        if(col_idx < 0 || col_idx > (ncols-1)) {
          Rf_error("in import_col: 1 <= %d <= %d needed", 
            col_idx, ncols);
        }

        REAL(VECTOR_ELT(dest, i))[j+col_idx*nrows] = REAL(col)[j];
        UNPROTECT(1);
      }
    }
  }

  UNPROTECT(3);
  return dest;
}
