#include "util.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <R.h>
#include <omp.h>
#include <numa.h>
#include <sched.h>

static const double eps = 0.0;

void x_times_beta(int *restrict x, double *restrict zz[], double *restrict beta, int *nRows, int *nVars, int *restrict numLevels, int *restrict catIndices, int *restrict contIndices, int *restrict catcatIndices, int *restrict contcontIndices, int *restrict catcontIndices, double *restrict result){
  int n = *nRows;
  int pCat=nVars[0], pCont=nVars[1], pCatCat=2*nVars[2], pContCont=2*nVars[3], pCatCont=2*nVars[4];
  int i, p, nLevels, allzero, offset = 0;
  int *restrict xOffsetPtr;
  double *restrict zOffsetPtr;
  double factor;

  // Setting up the localResults arrays
#ifdef __AVX__

  double *restrict*restrict localResOpt;
#pragma omp parallel
  {
    const int threadCount = omp_get_num_threads();
    const int threadId = omp_get_thread_num();
#pragma omp single
    localResOpt = malloc(threadCount * sizeof(double*));
#pragma omp barrier
    localResOpt[threadId] = _mm_malloc(n * sizeof(double), 64);
    memset(localResOpt[threadId], 0, n * sizeof(double));
  }

#else

  double *restrict localResMnt = malloc(n * sizeof(double));
  memset(localResMnt, 0, n * sizeof(double));

#endif

#pragma pomp inst begin(fista_x_times_beta)
  
  /* categorical */
  if (pCat > 0){
    Rf_error("categorical variables currently not supported. (can be restored!)");
    factor = sqrt(n);
    for (p=0; p<pCat; p++){
      nLevels = numLevels[catIndices[p]-1];
      /* check if beta for this variable is all zero. If so, move on */
      allzero = 1;
      for (i=0; i<nLevels; i++){
  if (fabs(beta[offset + i]) > eps){
    allzero = 0;
    break;
  }
      }
      if (allzero){
  offset += nLevels;
  continue;
      }
      xOffsetPtr = x + (catIndices[p]-1)*n;
      for (i=0; i<n; i++){
  result[i] += beta[offset + xOffsetPtr[i]] / factor;
      }
      offset += nLevels;
    }
  }
  
  /* continuous */
  if (pCont > 0){

#ifdef __AVX__

    int nDiv8 = n/8;

#pragma omp parallel private(p, i, zOffsetPtr)
    {

      const int threadId = omp_get_thread_num();
      
      int cpu = sched_getcpu();
      int node = numa_node_of_cpu(cpu);
      double *restrict z = zz[node];

      if(max_alignment((uintptr_t)z) < 64) {
        Rf_error("alignment of z is %d; need at least 64 byte alignment",
          max_alignment((uintptr_t)z));
      }
      if(max_alignment((uintptr_t)localResOpt[threadId]) < 64) {
        Rf_error("alignment of localResOpt is %d; need at least 64 byte alignment",
          max_alignment((uintptr_t)localResOpt[threadId]));
      }

#pragma omp for
      for (p=0; p<pCont; p++){
        
        int localOffset = offset + p;
        /* check if beta is zero */
        if (fabs(beta[localOffset]) >= eps) {
          zOffsetPtr = z + ((size_t)contIndices[p]-1)*n;
          
          __m256d beta0 = _mm256_broadcast_sd(beta+localOffset);
          for (i=0; i<nDiv8; ++i){
            __m256d z0 = _mm256_load_pd(zOffsetPtr+i*8);
            __m256d z1 = _mm256_load_pd(zOffsetPtr+i*8+4);
            __m256d res0 = _mm256_load_pd(localResOpt[threadId]+i*8);
            __m256d res1 = _mm256_load_pd(localResOpt[threadId]+i*8+4);
          
            __m256d prd0 = _mm256_fmadd_pd(z0, beta0, res0);
            __m256d prd1 = _mm256_fmadd_pd(z1, beta0, res1);

            _mm256_store_pd(localResOpt[threadId]+i*8,   prd0);
            _mm256_store_pd(localResOpt[threadId]+i*8+4, prd1);
          }
          
          for (i=nDiv8*8; i<n; ++i){
            localResOpt[threadId][i] += zOffsetPtr[i] * beta[localOffset];
          }
        }
      }
    }

#else

    int cpu = sched_getcpu();
    int node = numa_node_of_cpu(cpu);
    double *restrict z = zz[node];

    for (p=0; p<pCont; p++){
      int localOffset = offset + p;
      /* check if beta is zero */
      if (fabs(beta[localOffset]) >= eps) {
        zOffsetPtr = z + ((size_t)contIndices[p]-1)*n;
        for (i=0; i<n; i++) {
          localResMnt[i] += zOffsetPtr[i] * beta[localOffset];
        }
      }
    }

#endif

    offset += pCont;
  }
  
  /* categorical-categorical */
  if (pCatCat > 0){
    Rf_error("categorical variables currently not supported. (can be restored!)");
    int len;
    int *restrict yOffsetPtr;
    factor = sqrt(n);
    for (p=0; p<pCatCat; p+=2){
      nLevels = numLevels[catcatIndices[p]-1];
      len = nLevels * numLevels[catcatIndices[p+1]-1];
      /* check if beta is zero */
      allzero = 1;
      for (i=0; i<len; i++){
  if (fabs(beta[offset + i]) > eps){
    allzero = 0;
    break;
  }
      }
      if (allzero){
  offset += len;
  continue;
      }
      xOffsetPtr = x + (catcatIndices[p]-1)*n;
      yOffsetPtr = x + (catcatIndices[p+1]-1)*n;
      for (i=0; i<n; i++){
  result[i] += beta[offset + xOffsetPtr[i] + nLevels*yOffsetPtr[i]] / factor;
      }
      offset += len;
    }
  }
  
  /* continuous-continuous */
  if (pContCont > 0){

    factor = sqrt(3);

#ifdef __AVX__

    int nDiv8 = n/8;

#pragma omp parallel private(p, i, allzero, zOffsetPtr)
    {
      const int threadId = omp_get_thread_num();

      int cpu = sched_getcpu();
      int node = numa_node_of_cpu(cpu);
      double *restrict z = zz[node];

      double *restrict wOffsetPtr;
      double mean, norm;
      double *restrict prdOpt = _mm_malloc(n * sizeof *prdOpt, 64);

      if(max_alignment((uintptr_t)localResOpt[threadId]) < 64) {
        Rf_error("alignment of localResOpt is %d; need at least 64 byte alignment",
          max_alignment((uintptr_t)localResOpt[threadId]));
      }
      if(max_alignment((uintptr_t)z) < 64) {
        Rf_error("alignment of z is %d; need at least 64 byte alignment",
          max_alignment((uintptr_t)z));
      }
      if(max_alignment((uintptr_t)prdOpt) < 64) {
        Rf_error("alignment of prdOpt is %d; need at least 64 byte alignment",
          max_alignment((uintptr_t)prdOpt));
      }

#pragma omp for
      for (p=0; p<pContCont; p+=2){
        int localOffset = offset + 3 * (p / 2);
        /* check if beta is zero */
        allzero = 1;
        for (i=0; i<3; i++){
          if (fabs(beta[localOffset + i]) > eps){
            allzero = 0;
            break;
          }
        }
        if (!allzero){
          wOffsetPtr = z + ((size_t)contcontIndices[p]-1)*n;
          zOffsetPtr = z + ((size_t)contcontIndices[p+1]-1)*n;
          __m256d mean1 = _mm256_setzero_pd();
          __m256d norm1 = _mm256_setzero_pd();
          __m256d mean2 = _mm256_setzero_pd();
          __m256d norm2 = _mm256_setzero_pd();
          __m256d beta0 = _mm256_broadcast_sd(beta+localOffset);
          __m256d beta1 = _mm256_broadcast_sd(beta+localOffset+1);
          __m256d factr = _mm256_set1_pd(1.0/factor);
          for (i=0; i<nDiv8; ++i){
            __m256d z0 = _mm256_load_pd(zOffsetPtr+i*8);
            __m256d z1 = _mm256_load_pd(zOffsetPtr+i*8+4);
            __m256d w0 = _mm256_load_pd(wOffsetPtr+i*8);
            __m256d w1 = _mm256_load_pd(wOffsetPtr+i*8+4);

            __m256d wb0 = _mm256_mul_pd(w0, beta0);
            __m256d wb1 = _mm256_mul_pd(w1, beta0);
            __m256d wzb0 = _mm256_fmadd_pd(z0, beta1, wb0);
            __m256d wzb1 = _mm256_fmadd_pd(z1, beta1, wb1);

            __m256d res0 = _mm256_load_pd(localResOpt[threadId]+i*8);
            __m256d res1 = _mm256_load_pd(localResOpt[threadId]+i*8+4);

            __m256d wzbfr0 = _mm256_fmadd_pd(wzb0, factr, res0);
            __m256d wzbfr1 = _mm256_fmadd_pd(wzb1, factr, res1);

            _mm256_store_pd(localResOpt[threadId]+i*8,   wzbfr0);
            _mm256_store_pd(localResOpt[threadId]+i*8+4, wzbfr1);

            __m256d prod1 = _mm256_mul_pd(w0, z0);
            __m256d prod2 = _mm256_mul_pd(w1, z1);

            mean1 = _mm256_add_pd(mean1, prod1);
            mean2 = _mm256_add_pd(mean2, prod2);
            norm1 = _mm256_fmadd_pd(prod1, prod1, norm1);
            norm2 = _mm256_fmadd_pd(prod2, prod2, norm2);

            _mm256_store_pd(prdOpt+i*8,   prod1);
            _mm256_store_pd(prdOpt+i*8+4, prod2);
          }

          mean = sum_to_double(mean1, mean2);
          norm = sum_to_double(norm1, norm2);

          for (i=nDiv8*8; i<n; ++i){
            localResOpt[threadId][i] += (wOffsetPtr[i]*beta[localOffset] + zOffsetPtr[i]*beta[localOffset+1]) / factor;
            prdOpt[i] = wOffsetPtr[i] * zOffsetPtr[i];
            mean += prdOpt[i];
            norm += prdOpt[i]*prdOpt[i];
          }

          if (norm > 0){
            mean /= n;
            norm = sqrt(3 * (norm-n*mean*mean));
            double betanorm = beta[localOffset+2] / norm;
            __m256d bnr = _mm256_set1_pd(betanorm);          
            __m256d men = _mm256_set1_pd(mean);

            for (i=0; i<nDiv8; ++i){
              __m256d prod0 = _mm256_load_pd(prdOpt+i*8);
              __m256d prod1 = _mm256_load_pd(prdOpt+i*8+4);
              __m256d res0 = _mm256_load_pd(localResOpt[threadId]+i*8);
              __m256d res1 = _mm256_load_pd(localResOpt[threadId]+i*8+4);

              __m256d prmn0 = _mm256_sub_pd(prod0, men);
              __m256d prmn1 = _mm256_sub_pd(prod1, men);

              __m256d pmbnr0 = _mm256_fmadd_pd(prmn0, bnr, res0);
              __m256d pmbnr1 = _mm256_fmadd_pd(prmn1, bnr, res1);

              _mm256_store_pd(localResOpt[threadId]+i*8,   pmbnr0);
              _mm256_store_pd(localResOpt[threadId]+i*8+4, pmbnr1);
            }

            for (i=nDiv8*8; i<n; ++i){
              localResOpt[threadId][i] += (prdOpt[i]-mean) * betanorm;
            }
          }
        }
      }

      _mm_free(prdOpt);
    }

#else

    int cpu = sched_getcpu();
    int node = numa_node_of_cpu(cpu);
    double *restrict z = zz[node];

    double mean, norm;
    double *restrict wOffsetPtr;
    double *restrict prdMnt = malloc(n * sizeof *prdMnt);

    for (p=0; p<pContCont; p+=2){
      int localOffset = offset + 3 * (p / 2);
      /* check if beta is zero */
      allzero = 1;
      for (i=0; i<3; i++){
        if (fabs(beta[localOffset + i]) > eps){
          allzero = 0;
          break;
        }
      }
      if (!allzero){
        wOffsetPtr = z + ((size_t)contcontIndices[p]-1)*n;
        zOffsetPtr = z + ((size_t)contcontIndices[p+1]-1)*n;
        mean = norm = 0.0;
        for (i=0; i<n; i++){
          localResMnt[i] += (wOffsetPtr[i]*beta[localOffset] + zOffsetPtr[i]*beta[localOffset+1]) / factor;
          prdMnt[i] = wOffsetPtr[i] * zOffsetPtr[i];
          mean += prdMnt[i];
          norm += prdMnt[i]*prdMnt[i];
        }
        if (norm > 0){
          mean /= n;
          norm = sqrt(3 * (norm-n*pow(mean, 2)));
          for (i=0; i<n; i++){
            localResMnt[i] += (prdMnt[i]-mean) * beta[localOffset+2] / norm;
          }
        }
      }
    }

    free(prdMnt);

#endif

    offset += 3 * (pContCont / 2);
  }
  
  /* categorical-continuous */
  if (pCatCont > 0){
    Rf_error("categorical variables currently not supported. (can be restored!)");
    
    int cpu = sched_getcpu();
    int node = numa_node_of_cpu(cpu);
    double *restrict z = zz[node];

    factor = sqrt(2*n);
    double factorZ = sqrt(2);
    for (p=0; p<pCatCont; p+=2){
      nLevels = numLevels[catcontIndices[p]-1];
      /* check if beta is zero */
      allzero = 1;
      for (i=0; i<2*nLevels; i++){
  if (fabs(beta[offset + i]) > eps){
    allzero = 0;
    break;
  }
      }
      if (allzero){
  offset += 2*nLevels;
  continue;
      }
      xOffsetPtr = x + (catcontIndices[p]-1)*n;
      zOffsetPtr = z + (catcontIndices[p+1]-1)*n;
      for (i=0; i<n; i++){
  result[i] += beta[offset + xOffsetPtr[i]] / factor;
  result[i] += zOffsetPtr[i] * beta[offset + nLevels + xOffsetPtr[i]] / factorZ;
      }
      offset += 2*nLevels;
    }
  }

  /*double *restrict sum = malloc(n * sizeof *sum);
  memset(sum, 0, n * sizeof(double));
#pragma omp parallel private(i) 
  {
    const int threadCount = omp_get_num_threads();
#pragma omp for
    for (i=0; i<n; ++i) {
      for (int j=0; j<threadCount; j++)
        sum[i] += localResOpt[j][i];
    }
#pragma omp barrier
#pragma omp single
    for (i=0; i<n; ++i) {
      if (fabs(localResMnt[i]-sum[i]) > 1.0e-12) {
        Rprintf("%.20f\n%.20f\n%.20f\n%.20f\n\n", 
          localResMnt[i], sum[i],
          localResMnt[i]-sum[i], 1.0e-12);
        Rf_error("i've seen enough from x_times_beta");
      }
    }
  }
  free(sum);*/
  
  /* aggregate local results */
#ifdef __AVX__

#pragma omp parallel private(i) 
  {
    const int threadCount = omp_get_num_threads();
    const int threadId = omp_get_thread_num();
#pragma omp for
      for (i=0; i<n; i++) {
        for (int j=0; j<threadCount; j++)
            result[i] += localResOpt[j][i];
      }
      _mm_free(localResOpt[threadId]);
#pragma omp barrier
#pragma omp single 
      free((void*) localResOpt);
  }

#else

  for (i=0; i<n; i++) {
    result[i] += localResMnt[i];
  }
  free(localResMnt);

#endif
  
#pragma pomp inst end(fista_x_times_beta)

}

double compute_loglik(const double *restrict y, const double *restrict linear, const double *restrict intercept, const int *restrict nRows, const int *restrict family){
  double result = 0.0, mu = *intercept;
  int i, n = *nRows;
#pragma pomp inst begin(fista_compute_loglik)
  if (*family == 0){
    for (i=0; i<n; i++){
      result += pow(y[i]-mu-linear[i], 2);
    }
    result /= (2*n);
  }
  else {
    for (i=0; i<n; i++){
      result += -y[i]*(mu+linear[i]) + log(1+exp(mu+linear[i]));
    }
    result /= n;
  }
#pragma pomp inst end(fista_compute_loglik)
  return result;
}

void compute_objective(const double *restrict y, const double *restrict res, const double *restrict linear, const double *restrict intercept, const double *restrict beta, const int *restrict nRows, const int *restrict numGroups, const int *restrict groupSizes, const double *restrict lambda, double *restrict objValue, const int *restrict family){
  int i, j, size, n = *nRows, numgroups = *numGroups, offset = 0;
  double loglik = 0.0, penalty = 0.0, mu = *intercept, temp;
#pragma pomp inst begin(fista_compute_objective)
  if (*family == 0){
    for (i=0; i<n; i++) loglik += res[i]*res[i];
    loglik /= (2*n);
  }
  else {
    for (i=0; i<n; i++) loglik += -y[i]*(mu+linear[i]) + log(1+exp(mu+linear[i]));
    loglik /= n;
  }
  for (i=0; i<numgroups; i++){
    size = groupSizes[i];
    temp = 0.0;
    for (j=0; j<size; j++){
      temp += pow(beta[offset+j], 2);
    }
    penalty += sqrt(temp);
    offset += size;
  }
  *objValue = loglik + penalty*(*lambda);
#pragma pomp inst end(fista_compute_objective)
}

void compute_gradient(int *restrict x, double *restrict zz[], double *restrict r, int *restrict nRows, int *restrict nVars, int *restrict numLevels, int *restrict catIndices, int *restrict contIndices, int *restrict catcatIndices, int *restrict contcontIndices, int *restrict catcontIndices, double *restrict gradient){
  int n = *nRows;
  int pCat=nVars[0], pCont=nVars[1], pCatCat=2*nVars[2], pContCont=2*nVars[3], pCatCont=2*nVars[4];
  int i, p, offset = 0;
  int *restrict xOffsetPtr;
  double *restrict zOffsetPtr;
  double factor;

  size_t gradLength = pCont + 3 * pContCont/2;
#ifdef __AVX__
  double *restrict gradOpt = malloc(gradLength * sizeof *gradOpt);
#else
  double *restrict gradMnt = malloc(gradLength * sizeof *gradMnt);
#endif

#pragma pomp inst begin(fista_compute_gradient)
  /* categorical */
  if (pCat > 0){
    Rf_error("categorical variables currently not supported. (can be restored!)");
    factor = sqrt(n);
    for (p=0; p<pCat; p++){
      xOffsetPtr = x + (catIndices[p]-1)*n;
      for (i=0; i<n; i++){
  gradient[offset + xOffsetPtr[i]] += r[i];
      }
      offset += numLevels[catIndices[p]-1];
    }
    for (i=0; i<offset; i++){
      gradient[i] /= factor;
    }
  }

  /* continuous */
  if (pCont > 0){

#ifdef __AVX__

#pragma omp parallel private(p, i, zOffsetPtr)
    {

      int cpu = sched_getcpu();
      int node = numa_node_of_cpu(cpu);
      double *restrict z = zz[node];

      int nDiv8 = n/8;

      if(max_alignment((uintptr_t)z) < 64) {
        Rf_error("alignment of z is %d; need at least 64 byte alignment",
          max_alignment((uintptr_t)z));
      }
      if(max_alignment((uintptr_t)r) < 64) {
        Rf_error("alignment of r is %d; need at least 64 byte alignment",
          max_alignment((uintptr_t)r));
      }

#pragma omp for
      for (p=0; p<pCont; p++){

        int localOffset = offset + p;
        zOffsetPtr = z + ((size_t)contIndices[p]-1)*n;
        __m256d grd0 = _mm256_setzero_pd();
        __m256d grd1 = _mm256_setzero_pd();

        for (i=0; i<nDiv8; ++i){
          __m256d z0 = _mm256_load_pd(zOffsetPtr+i*8);
          __m256d z1 = _mm256_load_pd(zOffsetPtr+i*8+4);
          __m256d r0 = _mm256_load_pd(r+i*8);
          __m256d r1 = _mm256_load_pd(r+i*8+4);

          grd0 = _mm256_fmadd_pd(z0, r0, grd0);
          grd1 = _mm256_fmadd_pd(z1, r1, grd1);
        }

        double gradient0 = sum_to_double(grd0, grd1);

        for (i=nDiv8*8; i<n; ++i){
          gradient0 += zOffsetPtr[i] * r[i];
        }

        gradOpt[localOffset] = gradient[localOffset] + gradient0;
      }
    }

#else

    int cpu = sched_getcpu();
    int node = numa_node_of_cpu(cpu);
    double *restrict z = zz[node];

    for (p=0; p<pCont; p++){
      int localOffset = offset + p;
      zOffsetPtr = z + ((size_t)contIndices[p]-1)*n;
      double gradient0 = gradient[localOffset];
      for (i=0; i<n; i++){
        gradient0 += zOffsetPtr[i] * r[i];
      }
      gradMnt[localOffset] = gradient0;
    }

#endif

    offset += pCont;
  }

  /* categorical-categorical */
  if (pCatCat > 0){
    Rf_error("categorical variables currently not supported. (can be restored!)");
    factor = sqrt(n);
    int nLevels, start = offset;
    int *restrict yOffsetPtr;
    for (p=0; p<pCatCat; p+=2){
      xOffsetPtr = x + (catcatIndices[p]-1)*n;
      yOffsetPtr = x + (catcatIndices[p+1]-1)*n;
      nLevels = numLevels[catcatIndices[p]-1];
      for (i=0; i<n; i++){
  gradient[offset + xOffsetPtr[i] + nLevels*yOffsetPtr[i]] += r[i];
      }
      offset += nLevels * numLevels[catcatIndices[p+1]-1];
    }
    for (i=start; i<offset; i++){
      gradient[i] /= factor;
    }
  }

  /* continuous-continuous */
  if (pContCont > 0){
    factor = sqrt(3);

#ifdef __AVX__

    int nDiv8 = n/8;

#pragma omp parallel private(p, i, zOffsetPtr)
    {

      int cpu = sched_getcpu();
      int node = numa_node_of_cpu(cpu);
      double *restrict z = zz[node];

      double *restrict wOffsetPtr;
      double mean, norm;

      double prod;
      double rprd, rsum;

      if(max_alignment((uintptr_t)z) < 64) {
        Rf_error("alignment of z is %d; need at least 64 byte alignment",
          max_alignment((uintptr_t)z));
      }
      if(max_alignment((uintptr_t)r) < 64) {
        Rf_error("alignment of r is %d; need at least 64 byte alignment",
          max_alignment((uintptr_t)r));
      }

#pragma omp for
      for (p=0; p<pContCont; p+=2){
        int localOffset = offset + 3 * (p / 2);
        wOffsetPtr = z + ((size_t)contcontIndices[p]-1)*n;
        zOffsetPtr = z + ((size_t)contcontIndices[p+1]-1)*n;

        __m256d grad201 = _mm256_setzero_pd();
        __m256d grad202 = _mm256_setzero_pd();
        __m256d grad211 = _mm256_setzero_pd();
        __m256d grad212 = _mm256_setzero_pd();

        __m256d mean1 = _mm256_setzero_pd();
        __m256d mean2 = _mm256_setzero_pd();
        __m256d norm1 = _mm256_setzero_pd();
        __m256d norm2 = _mm256_setzero_pd();
        __m256d rprd1 = _mm256_setzero_pd();
        __m256d rprd2 = _mm256_setzero_pd();
        __m256d rsum1 = _mm256_setzero_pd();
        __m256d rsum2 = _mm256_setzero_pd();

        for (i=0; i<nDiv8; ++i){
          __m256d w1 = _mm256_load_pd(wOffsetPtr+i*8);
          __m256d w2 = _mm256_load_pd(wOffsetPtr+i*8+4);
          __m256d z1 = _mm256_load_pd(zOffsetPtr+i*8);
          __m256d z2 = _mm256_load_pd(zOffsetPtr+i*8+4);
          __m256d r1 = _mm256_load_pd(r+i*8);
          __m256d r2 = _mm256_load_pd(r+i*8+4);

          grad201 = _mm256_fmadd_pd(w1, r1, grad201);
          grad202 = _mm256_fmadd_pd(w2, r2, grad202);
          grad211 = _mm256_fmadd_pd(z1, r1, grad211);
          grad212 = _mm256_fmadd_pd(z2, r2, grad212);

          __m256d prod1 = _mm256_mul_pd(w1, z1);
          __m256d prod2 = _mm256_mul_pd(w2, z2);

          mean1 = _mm256_add_pd(mean1, prod1);
          mean2 = _mm256_add_pd(mean2, prod2);
          norm1 = _mm256_fmadd_pd(prod1, prod1, norm1);
          norm2 = _mm256_fmadd_pd(prod2, prod2, norm2);
          rprd1 = _mm256_fmadd_pd(prod1, r1, rprd1);
          rprd2 = _mm256_fmadd_pd(prod2, r2, rprd2);
          rsum1 = _mm256_add_pd(rsum1, r1);
          rsum2 = _mm256_add_pd(rsum2, r2);
        }
        
        double grad20 = sum_to_double(grad201, grad202);
        double grad21 = sum_to_double(grad211, grad212);
        mean = sum_to_double(mean1, mean2);
        norm = sum_to_double(norm1, norm2);
        rprd = sum_to_double(rprd1, rprd2);
        rsum = sum_to_double(rsum1, rsum2);

        for (i=nDiv8*8; i<n; ++i){
          grad20 += wOffsetPtr[i] * r[i];
          grad21 += zOffsetPtr[i] * r[i];
          prod = wOffsetPtr[i] * zOffsetPtr[i];
          mean += prod;
          norm += prod * prod;
          rprd += prod * r[i];
          rsum += r[i];
        }

        gradOpt[localOffset]   = (grad20+gradient[localOffset])   / factor;
        gradOpt[localOffset+1] = (grad21+gradient[localOffset+1]) / factor;

        if (norm > 0){
          mean = mean / n;
          norm = sqrt(3 * (norm-n*mean*mean));
          gradOpt[localOffset+2] = (gradient[localOffset+2] + (rprd - mean * rsum)) / norm;
        }
      }
    }

#else

    int cpu = sched_getcpu();
    int node = numa_node_of_cpu(cpu);
    double *restrict z = zz[node];

    double *restrict wOffsetPtr;
    double mean, norm;

    double *restrict product = malloc(n * sizeof *product);
    double gradient0, gradient1, gradient2;

    for (p=0; p<pContCont; p+=2){
      int localOffset = offset + 3 * (p / 2);
      wOffsetPtr = z + ((size_t)contcontIndices[p]-1)*n;
      zOffsetPtr = z + ((size_t)contcontIndices[p+1]-1)*n;
      mean = 0.0;
      norm = 0.0;
      gradient0 = gradient[localOffset];
      gradient1 = gradient[localOffset+1];
      gradient2 = gradient[localOffset+2];
      for (i=0; i<n; i++){
        gradient0 += wOffsetPtr[i] * r[i];
        gradient1 += zOffsetPtr[i] * r[i];
        product[i] = wOffsetPtr[i] * zOffsetPtr[i];
        mean += product[i];
        norm += product[i]*product[i];
      }
      gradient0 = gradient0 / factor;
      gradient1 = gradient1 / factor;
      gradMnt[localOffset]   = gradient0;
      gradMnt[localOffset+1] = gradient1;

      if (norm > 0){
        mean /= n;
        norm = sqrt(3 * (norm-n*pow(mean, 2)));
        for (i=0; i<n; i++){
          gradient2 += r[i]*(product[i]-mean);
        }
        gradient2 = gradient2 / norm;
        gradMnt[localOffset+2] = gradient2;
      }
    }

    free(product);

#endif

    offset += 3 * (pContCont / 2);

  }
  /* categorical-continuous */
  if (pCatCont > 0){
    Rf_error("categorical variables currently not supported. (can be restored!)");

    int cpu = sched_getcpu();
    int node = numa_node_of_cpu(cpu);
    double *restrict z = zz[node];

    factor = sqrt(2*n);
    double factorZ = sqrt(2);
    int nLevels;
    for (p=0; p<pCatCont; p+=2){
      nLevels = numLevels[catcontIndices[p]-1];
      xOffsetPtr = x + (catcontIndices[p]-1)*n;
      zOffsetPtr = z + (catcontIndices[p+1]-1)*n;
      for (i=0; i<n; i++){
  gradient[offset + xOffsetPtr[i]] += r[i];
  gradient[offset + nLevels + xOffsetPtr[i]] += zOffsetPtr[i] *r[i];
      }
      for (i=offset; i<offset+nLevels; i++){
  gradient[i] /= factor;
      }
      for (i=offset+nLevels; i<offset+2*nLevels; i++){
  gradient[i] /= factorZ;
      }
      offset += 2*nLevels;
    }
  }

  /*for (i=0; i<offset; ++i) {
    if (fabs(gradMnt[i]-gradOpt[i]) > 1.0e-12) {
      Rprintf("%.20f\n%.20f\n%.20f\n%.20f\n\n", 
        gradMnt[i], gradOpt[i],
        gradMnt[i]-gradOpt[i], 1.0e-12);
      Rf_error("i've seen enough from compute_gradient");
    }
  }*/

  /* normalize by n */
  for (i=0; i<offset; i++){
#ifdef __AVX__
    gradient[i] = gradOpt[i];
#else
    gradient[i] = gradMnt[i];
#endif
    gradient[i] /= -n;
  }

#pragma pomp inst end(fista_compute_gradient)

#ifdef __AVX__
  free(gradOpt);
#else
  free(gradMnt);
#endif

}

void compute_group_info(const int *nVars, const int *restrict numLevels, const int *restrict catIndices, const int *restrict contIndices, const int *restrict catcatIndices, const int *restrict contcontIndices, const int *restrict catcontIndices, int *length, int *groupSizes){
  int pCat=nVars[0], pCont=nVars[1], pCatCat=2*nVars[2], pContCont=2*nVars[3], pCatCont=2*nVars[4];
  int p, counter = 0, len = 0;
#pragma pomp inst begin(fista_compute_group_info)
  if (pCat > 0){
    for (p=0; p<pCat; p++){
      groupSizes[counter] = numLevels[catIndices[p]-1];
      len += groupSizes[counter++];
    }
  }
  if (pCont > 0){
    for (p=0; p<pCont; p++){
      groupSizes[counter++] = 1;
      ++len;
    }
  }
  if (pCatCat > 0){
    for (p=0; p<pCatCat; p+=2){
      groupSizes[counter] = numLevels[catcatIndices[p]-1] * numLevels[catcatIndices[p+1]-1];
      len += groupSizes[counter++];
    }
  }
  if (pContCont > 0){
    for (p=0; p<pContCont; p+=2){
      groupSizes[counter++] = 3;
      len += 3;
    }
  }
  if (pCatCont > 0){
    for (p=0; p<pCatCont; p+=2){
      groupSizes[counter] = 2 * numLevels[catcontIndices[p]-1];
      len += groupSizes[counter++];
    }
  }
  *length = len;
#pragma pomp inst end(fista_compute_group_info)
}

void compute_update(const double *restrict beta, double *restrict betaUpdated, const double *restrict gradient, const int *restrict groupSizes, const int *restrict numGroups, const double *restrict stepsize, const double *restrict lambda){
  int i, j, size, offset = 0, numgroups = *numGroups;
  double step = *stepsize, factor = step * (*lambda);
  double norm;
#pragma pomp inst begin(fista_compute_update)
  for (i=0; i<numgroups; i++){
    size = groupSizes[i];
    norm = 0.0;
    for (j=0; j<size; j++){
      betaUpdated[offset+j] = beta[offset+j] - step*gradient[offset+j];
      norm += betaUpdated[offset+j]*betaUpdated[offset+j];
    }
    norm = sqrt(norm);
    for (j=0; j<size; j++){
      betaUpdated[offset+j] *= fmax(0.0, 1.0 - factor/norm);
    }
    offset += size;
  }
#pragma pomp inst end(fista_compute_update)
}

void optimize_step(int *restrict x, double *restrict zz[], const double *restrict y, const double *restrict residual, double *restrict linear, int *restrict nRows, int *restrict numGroups, int *restrict groupSizes, int *restrict gradientLength, const double *restrict intercept, double *restrict beta, double *restrict betaUpdated, const double *restrict gradient, double *restrict stepsize, const double *restrict lambda, const double *restrict alpha, int *restrict nVars, int *restrict numLevels, int *restrict catIndices, int *restrict contIndices, int *restrict catcatIndices, int *restrict contcontIndices, int *restrict catcontIndices, const int *restrict family){
  int i, n = *nRows, len = *gradientLength;
  double step = *stepsize;
  double loglik, loglikUpdated;
  loglik = compute_loglik(y, linear, intercept, nRows, family);
  double *restrict delta = malloc(len * sizeof *delta);
  double gradientTimesDelta, deltaTimesDelta;
  double factor = *alpha;
  double mu = 0.0;
#pragma pomp inst begin(fista_optimize_step)
  while (1){
    gradientTimesDelta = 0.0;
    deltaTimesDelta = 0.0;
    compute_update(beta, betaUpdated, gradient, groupSizes, numGroups, &step, lambda);
    for (i=0; i<len; i++){
      delta[i] = betaUpdated[i] - beta[i];
      gradientTimesDelta += gradient[i] * delta[i];
      deltaTimesDelta += delta[i]*delta[i];
    }
    memset(linear, 0, n * sizeof *linear);
    if (*family == 0){/* gaussian case */
      x_times_beta(x, zz, delta, nRows, nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, linear);
      loglikUpdated = compute_loglik(residual, linear, &mu, nRows, family); /* residual already had intercept subtracted from it */
    }
    else {/* binomial case */
      x_times_beta(x, zz, betaUpdated, nRows, nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, linear);
      loglikUpdated = compute_loglik(y, linear, intercept, nRows, family);
    }
    if (loglikUpdated <= loglik + gradientTimesDelta + deltaTimesDelta/(2*step)) break;
    step *= factor;
  }
  *stepsize = step;
  free(delta);
#pragma pomp inst end(fista_optimize_step)
}

int check_convergence(const double *restrict beta, const double *restrict gradient, const int *restrict groupSizes, const int *restrict numGroups, const double *restrict lambda, const double *restrict tol){
  int i, j, size, allzero, offset = 0, numgroups = *numGroups;
  double norm, factor = *lambda, tolerance = *tol;
  for (i=0; i<numgroups; i++){
    size = groupSizes[i];
    /* check if beta is zero and compute norm of the gradient */
    allzero = 1;
    for (j=0; j<size; j++){
      if (fabs(beta[offset+j]) > eps){
  allzero = 0;
  break;
      }
    }
    norm = 0.0;
    for (j=0; j<size; j++){
      norm += gradient[offset+j]*gradient[offset+j];
    }
    norm = sqrt(norm);
    /* if nonzero, check for kkt condition */
    if (!allzero && fabs(norm-factor)/factor>tolerance) return 0;
    /* if zero, make sure gradient norm is < lambda */
    if (allzero && norm>factor) return 0;
    offset += size;
  }
  return 1;
} 

double update_theta(const double *restrict beta, const double *restrict intermediate, const double *restrict intermediateOld, const int gradientLength, const double theta){
  int i;
  double value = 0.0;
#pragma pomp inst begin(fista_update_theta)
  for (i=0; i<gradientLength; i++){
    value += (beta[i]-intermediate[i]) * (intermediate[i]-intermediateOld[i]);
  }
#pragma pomp inst end(fista_update_theta)
  return value > 0.0 ? 1.0 : theta;
}

void update_intercept(const double *restrict y, const int *restrict nRows, const double *restrict linear, double *restrict intercept, double *restrict residual, const int *restrict family){
  int i, n = *nRows;
  double residualMean = 0.0, mu = *intercept;
#pragma pomp inst begin(fista_update_intercept)
  if (*family == 0){
    for (i=0; i<n; i++){
      residual[i] = y[i] - mu - linear[i];
      residualMean += residual[i];
    }
    residualMean /= n;
    *intercept += residualMean;
    for (i=0; i<n; i++){
      residual[i] -= residualMean;
    }
  }
  else {
    const double xmax = -log(DBL_EPSILON);
    const double xmin = log(DBL_MIN);
    double *restrict temp = malloc(n * sizeof *temp);
    double *restrict exponent = malloc(n * sizeof *exponent);
    double f = 0.0, fPrime, sumY = 0.0, expMu;
    double sum;
    expMu = exp(-mu);
    for (i=0; i<n; i++){
      exponent[i] = exp(-linear[i]);
      temp[i] = expMu * exponent[i];
      sumY += y[i];
      sum = mu + linear[i];
      f += y[i] - (sum > xmax ? 1.0 : (sum < xmin ? 0.0 : 1/(1+temp[i])));
    }
    int iter = 0;
    while (iter<1000 && fabs(f)>1e-2){
      fPrime = 0.0;
      for (i=0; i<n; i++){
  sum = mu + linear[i];
  fPrime -= (sum > xmax || sum < xmin) ? 0.0 : temp[i]/pow(1+temp[i], 2);
      }
      mu -= f/fPrime;
      expMu = exp(-mu);
      f = sumY;
      for (i=0; i<n; i++){
  temp[i] = expMu * exponent[i];
  sum = mu + linear[i];
  f -= sum > xmax ? 1.0 : (sum < xmin ? 0.0 : 1/(1+temp[i]));
      }
      ++iter;
    }
    *intercept = mu;
    for (i=0; i<n; i++){
      sum = mu + linear[i];
      residual[i] = y[i] - (sum > xmax ? 1.0 : (sum < xmin ? 0.0 : 1/(1+temp[i])));
    }
    free(temp);
    free(exponent);
  }
#pragma pomp inst end(fista_update_intercept)
}      

double compute_stepsize(const double *restrict gradient, const double *restrict gradientOld, const double *restrict beta, const double *restrict betaOld, const int gradientLength){
  int i;
  double normBeta = 0.0, normGradient = 0.0;
#pragma pomp inst begin(fista_compute_stepsize)
  for (i=0; i<gradientLength; i++){
    normBeta += (beta[i]-betaOld[i])*(beta[i]-betaOld[i]);
    normGradient += (gradient[i]-gradientOld[i])*(gradient[i]-gradientOld[i]);
  }
#pragma pomp inst end(fista_compute_stepsize)
  return sqrt(normBeta/normGradient);
}

void gl_solver(int *restrict x, double *restrict zz[], double *restrict y, int *restrict nRows, double *restrict intercept, double *restrict beta, double *restrict residual, double *restrict linear, int *restrict numLevels, int *restrict nVars, int *restrict catIndices, int *restrict contIndices, int *restrict catcatIndices, int *restrict contcontIndices, int *restrict catcontIndices, double *restrict lambda, double *restrict tol, double *restrict alpha, int *restrict maxIter, int *restrict convergedFlag, double *restrict objValue, double *restrict steps, int *restrict family, Rboolean verbose){
#pragma pomp inst begin(fista_gl_solver)
  /* initialize required variables */
  struct timespec timer_x_times_beta,
                  timer_compute_gradient,
                  timer_optimize_step;
  double accum_timer_x_times_beta = 0;
  double accum_timer_compute_gradient = 0;
  double accum_timer_optimize_step = 0;
  int i, gradientLength, iter, converged, n = *nRows;
  int numGroups = nVars[0] + nVars[1] + nVars[2] + nVars[3] + nVars[4];
  int *restrict groupSizes = malloc(numGroups * sizeof *groupSizes);
  double theta, thetaOld, momentum, stepsize;
  /* get grouping information and initialize gradient */
  compute_group_info(nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, &gradientLength, groupSizes);
  double *restrict gradient = malloc(gradientLength * sizeof *gradient);
  double *restrict intermediate = malloc(gradientLength * sizeof *intermediate);
  double *restrict intermediateOld = malloc(gradientLength * sizeof *intermediateOld);
  memcpy(intermediate, beta, gradientLength * sizeof *beta);
  double *restrict betaOld = malloc(gradientLength * sizeof *betaOld);
  double *restrict gradientOld = malloc(gradientLength * sizeof *gradientOld);
  /* compute residual from initialized beta */
  if (verbose) {
    timer_x_times_beta = timer_start();
  }
  x_times_beta(x, zz, beta, nRows, nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, linear);
  if (verbose) {
    Rprintf("---> timer (x_times_beta): %lf [s]\n", timer_end(timer_x_times_beta));
  }
  update_intercept(y, nRows, linear, intercept, residual, family);
  /* start accelerated FISTA */
  iter = 0;
  theta = 1.0;
  *convergedFlag = 0;
  while (iter < *maxIter) {
    /* compute gradient */
    memcpy(gradientOld, gradient, gradientLength * sizeof *gradient);
    memset(gradient, 0, gradientLength * sizeof *gradient);
    if (verbose) {
      timer_compute_gradient = timer_start();
    }
    compute_gradient(x, zz, residual, nRows, nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, gradient);
    if (verbose) {
      accum_timer_compute_gradient += timer_end(timer_compute_gradient);
    }
    /* check convergence */
    converged = check_convergence(beta, gradient, groupSizes, &numGroups, lambda, tol);
    if (converged){
      *convergedFlag = 1;
      break;
    }
    /* compute intermediate update and stepsize */
    memcpy(intermediateOld, intermediate, gradientLength * sizeof *intermediate);
    stepsize = iter > 0 ? compute_stepsize(gradient, gradientOld, beta, betaOld, gradientLength) : 1.0;
    if (verbose) {
      timer_optimize_step = timer_start();
    }
    optimize_step(x, zz, y, residual, linear, nRows, &numGroups, groupSizes, &gradientLength, intercept, beta, intermediate, gradient, &stepsize, lambda, alpha, nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, family);
    if (verbose) {
      accum_timer_optimize_step += timer_end(timer_optimize_step);
    }
    /* check if restart required  and update theta*/
    thetaOld = update_theta(beta, intermediate, intermediateOld, gradientLength, theta);
    /* update momentum */
    theta = (1+sqrt(1+4*pow(thetaOld, 2))) / 2;
    momentum = (thetaOld-1) / theta;
    /* update beta */
    memcpy(betaOld, beta, gradientLength * sizeof *beta);
    for (i=0; i<gradientLength; i++){
      beta[i] = intermediate[i] + momentum*(intermediate[i]-intermediateOld[i]);
    }
    /* update residual and mu */
    memset(linear, 0, n * sizeof *linear);
    if (verbose) {
      timer_x_times_beta = timer_start();
    }
    x_times_beta(x, zz, beta, nRows, nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, linear);
    if (verbose) {
      accum_timer_x_times_beta += timer_end(timer_x_times_beta);
    }
    update_intercept(y, nRows, linear, intercept, residual, family);
    /* update iteration count */
    //compute_objective(y, residual, linear, intercept, beta, nRows, &numGroups, groupSizes, lambda, objValue+iter, family);
    steps[iter] = stepsize;
    ++iter;
  }
  if (verbose) {
    Rprintf("---> %d iterations:\n     compute_gradient: %lf\n     optimize_step:    %lf\n     x_times_beta:     %lf\n", 
      iter,
      accum_timer_compute_gradient,
      accum_timer_optimize_step,
      accum_timer_x_times_beta
    );
  }
  compute_objective(y, residual, linear, intercept, beta, nRows, &numGroups, groupSizes, lambda, objValue, family);
  free(groupSizes);
  free(gradient);
  free(intermediate);
  free(intermediateOld);
  free(betaOld);
  free(gradientOld);
#pragma pomp inst end(fista_gl_solver)
}

SEXP R_gl_solver(SEXP R_x, SEXP R_z, SEXP R_y, SEXP R_nRows, SEXP R_intercept, SEXP R_beta, SEXP R_linear, SEXP R_numLevels, SEXP R_nVars, SEXP R_catIndices, SEXP R_contIndices, SEXP R_catcatIndices, SEXP R_contcontIndices, SEXP R_catcontIndices, SEXP R_lambda, SEXP R_tol, SEXP R_alpha, SEXP R_maxIter, SEXP R_convergedFlag, SEXP R_objValue, SEXP R_steps, SEXP R_family, SEXP R_cpuInfo, SEXP R_verbose){
  PROTECT(R_x = coerceVector(R_x, INTSXP));
  //PROTECT(R_z = coerceVector(R_z, REALSXP));
  PROTECT(R_y = coerceVector(R_y, REALSXP));
  PROTECT(R_nRows = coerceVector(R_nRows, INTSXP));
  PROTECT(R_intercept = coerceVector(R_intercept, REALSXP));
  PROTECT(R_beta = coerceVector(R_beta, REALSXP));
  //PROTECT(R_residual = coerceVector(R_residual, REALSXP));
  PROTECT(R_linear = coerceVector(R_linear, REALSXP));
  PROTECT(R_numLevels = coerceVector(R_numLevels, INTSXP));
  PROTECT(R_nVars = coerceVector(R_nVars, INTSXP));
  PROTECT(R_catIndices = coerceVector(R_catIndices, INTSXP));
  PROTECT(R_contIndices = coerceVector(R_contIndices, INTSXP));
  PROTECT(R_catcatIndices = coerceVector(R_catcatIndices, INTSXP));
  PROTECT(R_contcontIndices = coerceVector(R_contcontIndices, INTSXP));
  PROTECT(R_catcontIndices = coerceVector(R_catcontIndices, INTSXP));
  PROTECT(R_lambda = coerceVector(R_lambda, REALSXP));
  PROTECT(R_tol = coerceVector(R_tol, REALSXP));
  PROTECT(R_alpha = coerceVector(R_alpha, REALSXP));
  PROTECT(R_maxIter = coerceVector(R_maxIter, INTSXP));
  PROTECT(R_convergedFlag = coerceVector(R_convergedFlag, INTSXP));
  PROTECT(R_objValue = coerceVector(R_objValue, REALSXP));
  PROTECT(R_steps = coerceVector(R_steps, REALSXP));
  PROTECT(R_family = coerceVector(R_family, INTSXP));
  PROTECT(R_verbose = coerceVector(R_verbose, LGLSXP));
  int *restrict x = INTEGER(R_x);
  //double *restrict z = REAL(R_z);
  double *restrict y = REAL(R_y);
  int *restrict nRows = INTEGER(R_nRows);
  double *restrict intercept = REAL(R_intercept);
  double *restrict beta = REAL(R_beta);
  //double *restrict residual = REAL(R_residual);
  double *restrict linear = REAL(R_linear);
  int *restrict numLevels = INTEGER(R_numLevels);
  int *restrict nVars = INTEGER(R_nVars);
  int *restrict catIndices = INTEGER(R_catIndices);
  int *restrict contIndices = INTEGER(R_contIndices);
  int *restrict catcatIndices = INTEGER(R_catcatIndices);
  int *restrict contcontIndices = INTEGER(R_contcontIndices);
  int *restrict catcontIndices = INTEGER(R_catcontIndices);
  double *restrict lambda = REAL(R_lambda);
  double *restrict tol = REAL(R_tol);
  double *restrict alpha = REAL(R_alpha);
  int *restrict maxIter = INTEGER(R_maxIter);
  int *restrict convergedFlag = INTEGER(R_convergedFlag);
  double *restrict objValue = REAL(R_objValue);
  double *restrict steps = REAL(R_steps);
  int *restrict family = INTEGER(R_family);

  SEXP R_node_used, R_max_num_nodes, R_num_used_cpus;
  PROTECT(R_node_used     = VECTOR_ELT(R_cpuInfo, 1));
  PROTECT(R_max_num_nodes = VECTOR_ELT(R_cpuInfo, 3));
  PROTECT(R_num_used_cpus = VECTOR_ELT(R_cpuInfo, 4));
  int *restrict node_used     = INTEGER(R_node_used);
  int *restrict max_num_nodes = INTEGER(R_max_num_nodes);
  int *restrict num_used_cpus = INTEGER(R_num_used_cpus);

  SEXP R_zz[max_num_nodes[0]];
  double *restrict zz[max_num_nodes[0]];
  for(int i=0; i<max_num_nodes[0]; ++i) {
    PROTECT(R_zz[i] = VECTOR_ELT(R_z, i*2+1));
    if (node_used[i] > 0) {
      zz[i] = REAL(R_zz[i]);
    }
    else {
      zz[i] = (double*)INTEGER(R_zz[i]);
    }
    UNPROTECT(1);
  }

  int cpu = sched_getcpu();
  int node = numa_node_of_cpu(cpu);

  SEXP R_residual;
  PROTECT(R_residual = alloc(64, *nRows, node));
  double *restrict residual = REAL(R_residual);

  Rboolean verbose = FALSE;
  if (LOGICAL(R_verbose)[0] == TRUE) verbose = TRUE;

  gl_solver(x, zz, y, nRows, intercept, beta, residual, linear, numLevels, nVars, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, lambda, tol, alpha, maxIter, convergedFlag, objValue, steps, family, verbose);
  
  UNPROTECT(26);
  SEXP result = PROTECT(allocVector(VECSXP, 4));
  SET_VECTOR_ELT(result, 0, R_intercept);
  SET_VECTOR_ELT(result, 1, R_beta);
  SET_VECTOR_ELT(result, 2, R_residual);
  SET_VECTOR_ELT(result, 3, R_objValue);
  const char *names[4] = {"mu", "coefficients", "res", "objValue"};
  SEXP sNames = PROTECT(allocVector(STRSXP, 4));
  int i;
  for (i=0; i<4; i++) SET_STRING_ELT(sNames, i, mkChar(names[i]));
  setAttrib(result, R_NamesSymbol, sNames);
  UNPROTECT(2);
  return result;
}
