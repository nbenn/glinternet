#ifndef UTIL_HEADER
#define UTIL_HEADER

#define _POSIX_C_SOURCE 199309L

#include <time.h>
#include <x86intrin.h>

static inline double sum_to_double(const __m256d x, const __m256d y) {
  const __m256d sum0 = _mm256_add_pd(x, y);
  const __m256d sum1 = _mm256_add_pd(sum0, _mm256_permute2f128_pd(sum0, sum0, 0x1));
  const __m128d sum2 = _mm_hadd_pd(_mm256_castpd256_pd128(sum1), _mm256_castpd256_pd128(sum1));
  return _mm_cvtsd_f64(sum2);
}

struct timespec timer_start();
double timer_end(struct timespec);

#endif