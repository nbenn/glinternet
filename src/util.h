#ifndef UTIL_HEADER
#define UTIL_HEADER

#define _GNU_SOURCE 1

#include <time.h>
#include <Rinternals.h>
#include <R_ext/Rallocators.h>
#include <stdint.h>

static inline int max_alignment(uintptr_t pointer) {
  int alignment = 2;
  while(pointer % alignment == 0) {
    alignment *= 2;
  }
  return alignment/2;
}

struct timespec timer_start();
double timer_end(struct timespec);

SEXP get_my_numa_node();
SEXP get_cpu_node_usage();

void* custom_alloc(R_allocator_t *allocator, size_t size);
void custom_free(R_allocator_t *allocator, void * addr);
SEXP alloc(int alignment, R_xlen_t length, size_t node);

SEXP alloc_z(SEXP a, SEXP b, SEXP x, SEXP info);
SEXP retry_alloc_z(SEXP z, SEXP info);

SEXP extract_col(SEXP src, SEXP i);
SEXP import_col(SEXP col, SEXP dest, SEXP i, SEXP info);

SEXP alloc_res(SEXP y);
SEXP copy_vec(SEXP src, SEXP dst) ;

#endif