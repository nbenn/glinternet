#include "util.h"

#include <stdbool.h>
#include <stdlib.h>

typedef struct allocator_data {
    size_t offset;
    size_t alignment;
} allocator_data;

// call this function to start a nanosecond-resolution timer
struct timespec timer_start(void){
    struct timespec start_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
    return start_time;
}

// call this function to end a timer, returning seconds elapsed as double
double timer_end(struct timespec start_time){
    struct timespec end_time;
    double accum;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time);
    accum = ( end_time.tv_sec - start_time.tv_sec )
          + ( end_time.tv_nsec - start_time.tv_nsec )
            / (double)1000000000;
    return accum;
}

void* custom_alloc(R_allocator_t *allocator, size_t size) {
  //size_t alignment = ((allocator_data*)allocator->data)->alignment;
  /*size_t offset = ((allocator_data*)allocator->data)->offset;
  void* orig_addr = _mm_malloc(size+offset, 4096);
  Rprintf("original addr: %p\n", orig_addr);
  Rprintf("max alig     : %d\n", max_alignment((uintptr_t)orig_addr));
  void* shifted_addr = (void*)((uintptr_t)orig_addr + offset);
  Rprintf("shifted addr : %p\n", shifted_addr);
  Rprintf("max alig     : %d\n", max_alignment((uintptr_t)shifted_addr));
  return (void*) shifted_addr;*/
  return (void*) _mm_malloc(size, 64);
}

void custom_free(R_allocator_t *allocator, void* shifted_addr) {
  /*size_t offset = ((allocator_data*)allocator->data)->offset;
  Rprintf("shifted addr : %p\n", shifted_addr);
  Rprintf("max alig     : %d\n", max_alignment((uintptr_t)shifted_addr));
  void* orig_addr = (void*)((uintptr_t)shifted_addr-offset);
  Rprintf("original addr: %p\n", orig_addr);
  Rprintf("max alig     : %d\n", max_alignment((uintptr_t)orig_addr));
  _mm_free(orig_addr);*/
  _mm_free(shifted_addr);
}

SEXP alloc(int alignment, R_xlen_t length) {
  allocator_data* custom_allocator_data = malloc(sizeof(allocator_data));
  custom_allocator_data->offset = 56;
  custom_allocator_data->alignment = alignment;

  R_allocator_t* custom_allocator = malloc(sizeof(R_allocator_t));
  custom_allocator->mem_alloc = &custom_alloc;
  custom_allocator->mem_free = &custom_free;
  custom_allocator->res = NULL;
  custom_allocator->data = custom_allocator_data;

  SEXP result = PROTECT(allocVector3(REALSXP, length, custom_allocator));

  /*uintptr_t addr = (uintptr_t)REAL(result);
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
  }*/

  UNPROTECT(1);
  return result;
}

SEXP alloc_z(SEXP a, SEXP b, SEXP x) {
  R_xlen_t n = xlength(x);
  int nrows = asInteger(a);
  int ncols = asInteger(b);
  int alignment = 64;
  SEXP result;

  int i, j;
  
  if((R_xlen_t)nrows * (R_xlen_t)ncols != n) {
    Rf_error("dimensions mismatch: %d rows times %d cols != %d elements", 
      nrows, ncols, n);
  }

  PROTECT(result = alloc(alignment, n));

  // populate used matrices with normalized x values
  for (i=0; i<ncols; ++i) {
    double first = REAL(x)[(size_t)i*nrows];      
    bool allEqual = false;
    for (j=1; j<nrows; ++j) {
      if (REAL(x)[j+(size_t)i*nrows] != first) {
        break;
      }
      else if (j == nrows-1) allEqual = true;
    }
    if (allEqual) {
      for (j=0; j<nrows; ++j) {
        REAL(result)[j+(size_t)i*nrows] = 0;
      }
      continue;
    }
    double mean = 0;
    double *restrict newvec = malloc(nrows * sizeof (double));
    for (j=0; j<nrows; ++j) {
      mean += REAL(x)[j+(size_t)i*nrows];
    }
    mean /= nrows;
    for (j=0; j<nrows; ++j) {
      newvec[j] = REAL(x)[j+(size_t)i*nrows] - mean;
    }
    double norm = 0;
    for (j=0; j<nrows; ++j) {
      norm += newvec[j]*newvec[j];
    }
    norm = 1/sqrt(norm);
    for (int j=0; j<nrows; ++j) {
      REAL(result)[j+(size_t)i*nrows] = newvec[j] * norm;
    }
    free(newvec);
  }

  // make 2d array out of used matrices
  SEXP dims2;
  PROTECT(dims2 = allocVector(INTSXP, 2));
  INTEGER(dims2)[0] = nrows;
  INTEGER(dims2)[1] = ncols;
  setAttrib(result, R_DimSymbol, dims2);

  UNPROTECT(2);
  return result;
}

SEXP alloc_res(SEXP y) {
  R_xlen_t nrows = length(y);
  int alignment = 64;
  SEXP result;

  PROTECT(result = alloc(alignment, nrows));

  double mean = 0;
    
  for (int i=0; i<nrows; ++i) {
    mean += REAL(y)[i];
  }
  mean /= nrows;
  for (int i=0; i<nrows; ++i) {
    REAL(result)[i] = REAL(y)[i]-mean;
  }

  UNPROTECT(1);
  return result;
}

SEXP copy_vec(SEXP src, SEXP dst) {
  R_xlen_t n_src = xlength(src);
  R_xlen_t n_dst = xlength(dst);

  if(n_dst != n_dst) {
    Rf_error("dimensions mismatch: %d (source) != %d (destination)", 
      n_src, n_dst);
  }

  PROTECT(src = coerceVector(src, REALSXP));
  PROTECT(dst = coerceVector(dst, REALSXP));

  double *restrict source = REAL(src);
  double *restrict dstntn = REAL(dst);

  for (size_t i = 0; i < n_src; ++i) {
    dstntn[i] = source[i];
  }

  UNPROTECT(2);
  return dst;
}

