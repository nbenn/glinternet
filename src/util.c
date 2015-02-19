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

// call this function to end a timer, returning elapsed seconds
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
  size_t offset = ((allocator_data*)allocator->data)->offset;
  void* orig_addr = _mm_malloc(size+offset, 4096);
  //Rprintf("original addr: %p\n", orig_addr);
  //Rprintf("max alig     : %d\n", max_alignment((uintptr_t)orig_addr));
  void* shifted_addr = (void*)((uintptr_t)orig_addr + offset);
  //Rprintf("shifted addr : %p\n", shifted_addr);
  //Rprintf("max alig     : %d\n", max_alignment((uintptr_t)shifted_addr));
  return (void*) shifted_addr;
  //return (void*) _mm_malloc(size, 64);
}

void custom_free(R_allocator_t *allocator, void* shifted_addr) {
  size_t offset = ((allocator_data*)allocator->data)->offset;
  //Rprintf("shifted addr : %p\n", shifted_addr);
  //Rprintf("max alig     : %d\n", max_alignment((uintptr_t)shifted_addr));
  void* orig_addr = (void*)((uintptr_t)shifted_addr-offset);
  //Rprintf("original addr: %p\n", orig_addr);
  //Rprintf("max alig     : %d\n", max_alignment((uintptr_t)orig_addr));
  _mm_free(orig_addr);
  //_mm_free(shifted_addr);
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

  UNPROTECT(1);
  return result;
}

SEXP alloc_z(SEXP a, SEXP b, SEXP x) {
  R_xlen_t n = xlength(x);
  R_xlen_t nrows = (R_xlen_t)asInteger(a);
  R_xlen_t ncols = (R_xlen_t)asInteger(b);
  int alignment = 64;
  SEXP result;

  if(nrows * ncols != n) {
    Rf_error("dimensions mismatch: %d rows times %d cols != %d elements", 
      nrows, ncols, n);
  }

  PROTECT(result = alloc(alignment, n));

  // make 2d array out of double array
  SEXP dims2;
  PROTECT(dims2 = allocVector(INTSXP, 2));
  INTEGER(dims2)[0] = nrows;
  INTEGER(dims2)[1] = ncols;
  setAttrib(result, R_DimSymbol, dims2);

  UNPROTECT(2);
  return result;
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

SEXP import_col(SEXP col, SEXP dest, SEXP i) {
  SEXP dims2;
  PROTECT(dims2 = getAttrib(dest, R_DimSymbol));
  R_xlen_t nrows = (R_xlen_t)INTEGER(dims2)[0];
  R_xlen_t ncols = (R_xlen_t)INTEGER(dims2)[1];

  R_xlen_t nelem = xlength(col);

  if(nrows != nelem) {
    Rf_error("in import_col: %d != %d", nrows, nelem);
  }

  PROTECT(i = coerceVector(i, INTSXP));
  int col_idx = *INTEGER(i);

  if(col_idx < 1 || col_idx > ncols) {
    Rf_error("in import_col: 1 <= %d <= %d needed", 
      col_idx, ncols);
  }

  col_idx -= 1;

  for (R_xlen_t i = 0; i < nrows; ++i) {
    REAL(dest)[i+col_idx*nrows] = REAL(col)[i];
  }

  UNPROTECT(2);
  return dest;
}
