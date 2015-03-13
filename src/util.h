#ifndef UTIL_HEADER
#define UTIL_HEADER

#define _GNU_SOURCE 1

#include <time.h>
#include <Rinternals.h>
#include <R_ext/Rallocators.h>
#include <stdint.h>
#include <x86intrin.h>

static inline double sum_to_double(const __m256d x, const __m256d y) {
  const __m256d sum0 = _mm256_add_pd(x, y);
  const __m256d sum1 = _mm256_add_pd(sum0, _mm256_permute2f128_pd(sum0, sum0, 0x1));
  const __m128d sum2 = _mm_hadd_pd(_mm256_castpd256_pd128(sum1), _mm256_castpd256_pd128(sum1));
  return _mm_cvtsd_f64(sum2);
}

static inline int max_alignment(uintptr_t pointer) {
  int alignment = 2;
  while(pointer % alignment == 0) {
    alignment *= 2;
  }
  return alignment/2;
}

static inline float sum_to_float(const __m256 x) {
  const __m128 hiQuad = _mm256_extractf128_ps(x, 1);
  const __m128 loQuad = _mm256_castps256_ps128(x);
  const __m128 sumQuad = _mm_add_ps(loQuad, hiQuad);
  const __m128 loDual = sumQuad;
  const __m128 hiDual = _mm_movehl_ps(sumQuad, sumQuad);
  const __m128 sumDual = _mm_add_ps(loDual, hiDual);
  const __m128 lo = sumDual;
  const __m128 hi = _mm_shuffle_ps(sumDual, sumDual, 0x1);
  const __m128 sum = _mm_add_ss(lo, hi);
  return _mm_cvtss_f32(sum);
}

struct timespec timer_start();
double timer_end(struct timespec);

SEXP get_my_numa_node();
SEXP get_cpu_node_usage();

void* custom_alloc(R_allocator_t *allocator, size_t size);
void custom_free(R_allocator_t *allocator, void * addr);
SEXP alloc(int alignment, R_xlen_t length, int node);

static void finalize_singles_array_pointer(SEXP ext);

SEXP extract_col(SEXP src, SEXP i);
SEXP import_col(SEXP col, SEXP dest, SEXP i, SEXP info);

SEXP alloc_z_interleaved(SEXP a, SEXP b, SEXP x, SEXP info);
SEXP alloc_z_float(SEXP z, SEXP info);

#endif
