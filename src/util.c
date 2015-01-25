#include "util.h"

#include <stdbool.h>
#include <stdlib.h>
#include <R_ext/Applic.h> /* for dgemm */

#ifdef HAVE_LONG_DOUBLE
# define LDOUBLE long double
#else
# define LDOUBLE double
#endif

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
  size_t offset = ((allocator_data*)allocator->data)->offset;
  void* orig_addr = _mm_malloc(size+offset, 4096);
  Rprintf("original addr: %p\n", orig_addr);
  Rprintf("max alig     : %d\n", max_alignment((uintptr_t)orig_addr));
  void* shifted_addr = (void*)((uintptr_t)orig_addr + offset);
  Rprintf("shifted addr : %p\n", shifted_addr);
  Rprintf("max alig     : %d\n", max_alignment((uintptr_t)shifted_addr));
  return (void*) shifted_addr;
  //return (void*) _mm_malloc(size, 64);
}

void custom_free(R_allocator_t *allocator, void* shifted_addr) {
  size_t offset = ((allocator_data*)allocator->data)->offset;
  Rprintf("shifted addr : %p\n", shifted_addr);
  Rprintf("max alig     : %d\n", max_alignment((uintptr_t)shifted_addr));
  void* orig_addr = (void*)((uintptr_t)shifted_addr-offset);
  Rprintf("original addr: %p\n", orig_addr);
  Rprintf("max alig     : %d\n", max_alignment((uintptr_t)orig_addr));
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

  /*LDOUBLE s, t;
  R_xlen_t i, j;
  double mean;

  // populate used matrices with normalized x values
  for (i=0; i<ncols; ++i) {
    double first = REAL(x)[i*nrows];      
    bool allEqual = false;
    for (j=1; j<nrows; ++j) {
      if (REAL(x)[j+i*nrows] != first) {
        break;
      }
      else if (j == nrows-1) allEqual = true;
    }
    if (allEqual) {
      for (j=0; j<nrows; ++j) {
        REAL(result)[j+i*nrows] = 0.;
      }
      continue;
    }
    s = 0.;
    t = 0.;
    double *restrict newvec = malloc(nrows * sizeof (double));
    // calculation of mean according to R/src/main/summary.c
    // see function do_summary(...)
    for (j = 0; j < nrows; ++j) s += REAL(x)[j+i*nrows];
    s /= nrows;
    if(R_FINITE((double)s)) {
      for (j = 0; j < nrows; ++j) t += (REAL(x)[j+i*nrows] - s);
      s += t/nrows;
    }
    mean = (double) s;

    for (j=0; j<nrows; ++j) {
      newvec[j] = REAL(x)[j+i*nrows] - mean;
    }

    // calculate inner prod of norm col vec (newvec)
    int nrx = 1, ncx = nrows, nry = nrows, ncy = 1;
    char *transa = "N", *transb = "N";
    double ans = 0.0, one = 1.0, zero = 0.0;

    F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
        newvec, &nrx, newvec, &nry, &zero, &ans, &nrx);

    double norm = 1/sqrt(ans);
    for (int j=0; j<nrows; ++j) {
      REAL(result)[j+i*nrows] = newvec[j] * norm;
    }
    free(newvec);
  }*/

  // make 2d array out of result
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

SEXP alloc_res(SEXP y) {

  //Rprintf("allocating vector\n");

  R_xlen_t n = XLENGTH(y);
  SEXP result;
  int alignment = 64;
  PROTECT(result = alloc(alignment, n));

  LDOUBLE s = 0., t = 0.;
  R_xlen_t i;

  // calculation of mean according to R/src/main/summary.c
  // see function do_summary(...)
  for (i = 0; i < n; ++i) s += REAL(y)[i];
  s /= n;
  if(R_FINITE((double)s)) {
    for (i = 0; i < n; ++i) t += (REAL(y)[i] - s);
    s += t/n;
  }

  double mean_val = (double) s;

  for (i = 0; i < n; ++i) {
    REAL(result)[i] = REAL(y)[i]-mean_val;
  }

  UNPROTECT(1);
  return result;
}

SEXP copy_vec(SEXP src, SEXP dst) {
  R_xlen_t n_src = xlength(src);
  R_xlen_t n_dst = xlength(dst);

  //Rprintf("copying vector\n");

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

