#include "util.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#ifdef _OPENMP
# include <omp.h>
#endif

static const double eps = 0.0;

void retrieve_beta(const double *restrict beta, const int *restrict groupSizes, const int *restrict numGroups, int *restrict idx, int *restrict betaIdx){
  int i, p, offset = 0, numgroups = *numGroups, size;
  for (p=0; p<numgroups; p++){
    size = groupSizes[p];
    for (i=0; i<size; i++){
      if (fabs(beta[offset+i]) > eps){
	memset(betaIdx+offset, 1, size * sizeof *betaIdx);
	idx[p] = 1;
	break;
      }
    }
    offset += size;
  }
}

SEXP R_retrieve_beta(SEXP R_beta, SEXP R_groupSizes, SEXP R_numGroups, SEXP R_idx, SEXP R_betaIdx){
  PROTECT(R_beta = coerceVector(R_beta, REALSXP));
  PROTECT(R_groupSizes = coerceVector(R_groupSizes, INTSXP));
  PROTECT(R_numGroups = coerceVector(R_numGroups, INTSXP));
  PROTECT(R_idx = coerceVector(R_idx, INTSXP));
  PROTECT(R_betaIdx = coerceVector(R_betaIdx, INTSXP));
  double *restrict beta = REAL(R_beta);
  int *restrict groupSizes = INTEGER(R_groupSizes);
  int *restrict numGroups = INTEGER(R_numGroups);
  int *restrict idx = INTEGER(R_idx);
  int *restrict betaIdx = INTEGER(R_betaIdx);
  retrieve_beta(beta, groupSizes, numGroups, idx, betaIdx);
  UNPROTECT(5);
  SEXP result = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(result, 0, R_idx);
  SET_VECTOR_ELT(result, 1, R_betaIdx);
  const char *names[2] = {"idx", "betaIdx"};
  SEXP sNames = PROTECT(allocVector(STRSXP, 2));
  int i;
  for (i=0; i<2; i++) SET_STRING_ELT(sNames, i, mkChar(names[i]));
  setAttrib(result, R_NamesSymbol, sNames);
  UNPROTECT(2);
  return result;
}

int compare(const void *x, const void *y){
  return *(int *)x - *(int *)y;
}

int match_pair(const int *restrict x, const int *restrict y){
  if ((x[0]==y[0] && x[1]==y[1]) || (x[0]==y[1] && x[1]==y[0])) return 1;
  return 0;
}

int match_pair_catcont(const int *restrict x, const int *restrict y){
  return (x[0]==y[0] && x[1]==y[1]) ? 1 : 0;
}

void initialize_beta(double *restrict beta, const double *restrict betaOld, const int *restrict nVars, const int *restrict nVarsOld, const int *restrict cumGroupSizes, const int *restrict numLevels, const int *restrict catIndices, const int *restrict catIndicesOld, const int *restrict contIndices, const int *restrict contIndicesOld, const int *restrict catcatIndices, const int *restrict catcatIndicesOld, const int *restrict contcontIndices, const int *restrict contcontIndicesOld, const int *restrict catcontIndices, const int *restrict catcontIndicesOld){
  int i, p, size, offset = 1;
  int pCat = nVars[0], pCont = nVars[1], pCatCat = 2*nVars[2], pContCont = 2*nVars[3], pCatCont = 2*nVars[4];
  int pCatOld = nVarsOld[0], pContOld = nVarsOld[1], pCatCatOld = 2*nVarsOld[2], pContContOld = 2*nVarsOld[3], pCatContOld = 2*nVarsOld[4];
  int *restrict match;
  /* intercept */
  beta[0] = betaOld[0];
  if (pCat>0){
    for (p=0; p<pCat; p++){
      size = numLevels[catIndices[p]-1];
      if (pCatOld > 0){
	match = bsearch(catIndices+p, catIndicesOld, pCatOld, sizeof *catIndicesOld, compare);
	if (match != NULL){
	  memcpy(beta+offset, betaOld+cumGroupSizes[match-catIndicesOld], size * sizeof *betaOld);
	}
      }
      offset += size;
    }
  }
  if (pCont > 0){
    for (p=0; p<pCont; p++){
      if (pContOld > 0){
	match = bsearch(contIndices+p, contIndicesOld, pContOld, sizeof *contIndicesOld, compare);
	if (match != NULL){
	  beta[offset] = betaOld[cumGroupSizes[pCatOld + (match-contIndicesOld)]];
	}
      }
      ++offset;
    }
  }
  if (pCatCat > 0){
    for (p=0; p<pCatCat; p+=2){
      size = numLevels[catcatIndices[p]-1] * numLevels[catcatIndices[p+1]-1];
      if (pCatCatOld > 0){
	for (i=0; i<pCatCatOld; i+=2){
	  if (match_pair(catcatIndices+p, catcatIndicesOld+i)){
	    memcpy(beta+offset, betaOld+cumGroupSizes[pCatOld+pContOld+i/2], size * sizeof *betaOld);
	    break;
	  }
	}
      }
      offset += size;
    }
  }
  if (pContCont > 0){
    for (p=0; p<pContCont; p+=2){
      if (pContContOld > 0){
	for (i=0; i<pContContOld; i+=2){
	  if (match_pair(contcontIndices+p, contcontIndicesOld+i)){
	    memcpy(beta+offset, betaOld+cumGroupSizes[pCatOld+pContOld+(pCatCatOld+i)/2], 3 * sizeof *betaOld);
	    break;
	  }
	}
      }
      offset += 3;
    }
  }
  if (pCatCont > 0){
    for (p=0; p<pCatCont; p+=2){
      size = 2 * numLevels[catcontIndices[p]-1];
      if (pCatContOld > 0){
	for (i=0; i<pCatContOld; i+=2){
	  if (match_pair_catcont(catcontIndices+p, catcontIndicesOld+i)){
	    memcpy(beta+offset, betaOld+cumGroupSizes[pCatOld+pContOld+(pCatCatOld+pContContOld+i)/2], size * sizeof *betaOld);
	    break;
	  }
	}
      }
      offset += size;
    }
  }
}

SEXP R_initialize_beta(SEXP R_beta, SEXP R_betaOld, SEXP R_nVars, SEXP R_nVarsOld, SEXP R_cumGroupSizes, SEXP R_numLevels, SEXP R_catIndices, SEXP R_catIndicesOld, SEXP R_contIndices, SEXP R_contIndicesOld, SEXP R_catcatIndices, SEXP R_catcatIndicesOld, SEXP R_contcontIndices, SEXP R_contcontIndicesOld, SEXP R_catcontIndices, SEXP R_catcontIndicesOld){
  PROTECT(R_beta = coerceVector(R_beta, REALSXP));
  PROTECT(R_betaOld = coerceVector(R_betaOld, REALSXP));
  PROTECT(R_nVars = coerceVector(R_nVars, INTSXP));
  PROTECT(R_nVarsOld = coerceVector(R_nVarsOld, INTSXP));
  PROTECT(R_cumGroupSizes = coerceVector(R_cumGroupSizes, INTSXP));
  PROTECT(R_numLevels = coerceVector(R_numLevels, INTSXP));
  PROTECT(R_catIndices = coerceVector(R_catIndices, INTSXP));
  PROTECT(R_catIndicesOld = coerceVector(R_catIndicesOld, INTSXP));
  PROTECT(R_contIndices = coerceVector(R_contIndices, INTSXP));
  PROTECT(R_contIndicesOld = coerceVector(R_contIndicesOld, INTSXP));
  PROTECT(R_catcatIndices = coerceVector(R_catcatIndices, INTSXP));
  PROTECT(R_catcatIndicesOld = coerceVector(R_catcatIndicesOld, INTSXP));
  PROTECT(R_contcontIndices = coerceVector(R_contcontIndices, INTSXP));
  PROTECT(R_contcontIndicesOld = coerceVector(R_contcontIndicesOld, INTSXP));
  PROTECT(R_catcontIndices = coerceVector(R_catcontIndices, INTSXP));
  PROTECT(R_catcontIndicesOld = coerceVector(R_catcontIndicesOld, INTSXP));
  double *restrict beta = REAL(R_beta);
  double *restrict betaOld = REAL(R_betaOld);
  int *restrict nVars = INTEGER(R_nVars);
  int *restrict nVarsOld = INTEGER(R_nVarsOld);
  int *restrict cumGroupSizes = INTEGER(R_cumGroupSizes);
  int *restrict numLevels = INTEGER(R_numLevels);
  int *restrict catIndices = INTEGER(R_catIndices);
  int *restrict catIndicesOld = INTEGER(R_catIndicesOld);
  int *restrict contIndices = INTEGER(R_contIndices);
  int *restrict contIndicesOld = INTEGER(R_contIndicesOld);
  int *restrict catcatIndices = INTEGER(R_catcatIndices);
  int *restrict catcatIndicesOld = INTEGER(R_catcatIndicesOld);
  int *restrict contcontIndices = INTEGER(R_contcontIndices);
  int *restrict contcontIndicesOld = INTEGER(R_contcontIndicesOld);
  int *restrict catcontIndices = INTEGER(R_catcontIndices);
  int *restrict catcontIndicesOld = INTEGER(R_catcontIndicesOld);
  initialize_beta(beta, betaOld, nVars, nVarsOld, cumGroupSizes, numLevels, catIndices, catIndicesOld, contIndices, contIndicesOld, catcatIndices, catcatIndicesOld, contcontIndices, contcontIndicesOld, catcontIndices, catcontIndicesOld);
  UNPROTECT(16);
  return R_beta;
}

void rescale_beta(int *restrict x, double *restrict z, const int *restrict nRows, const double *restrict beta, const int *restrict betaLen, const int *restrict nVars, const int *restrict numLevels, const int *restrict catIndices, const int *restrict contIndices, const int *restrict catcatIndices, const int *restrict contcontIndices, const int *restrict catcontIndices, double *restrict result){
  int i, p, size, n = *nRows, offset = 1;
  double factor, mean, norm;
  double *restrict zOffsetPtr;
  int pCat = nVars[0], pCont = nVars[1], pCatCat = 2*nVars[2], pContCont = 2*nVars[3], pCatCont = 2*nVars[4];
  memcpy(result, beta, *betaLen * sizeof *beta);
  if (pCat + pCont + pCatCat + pContCont + pCatCont == 0) return;
  if (pCat > 0){
    factor = sqrt(n);
    for (p=0; p<pCat; p++){
      size = numLevels[catIndices[p]-1];
      for (i=0; i<size; i++){
	result[offset + i] /= factor;
      }
      offset += size;
    }
  }
  if (pCont > 0){
    for (p=0; p<pCont; p++){
      zOffsetPtr = z + (contIndices[p]-1)*n;
      mean = 0.0;
      norm = 0.0;
      for (i=0; i<n; i++){
	mean += zOffsetPtr[i];
	norm += zOffsetPtr[i]*zOffsetPtr[i];
      }
      mean /= n;
      norm = (fabs(norm) > 1e-30) ? sqrt(norm-n*pow(mean, 2)) : 1.0;
      result[offset] /= norm;
      result[0] -= mean * result[offset];
      ++offset;
    }
  }
  if (pCatCat > 0){
    factor = sqrt(n);
    for (p=0; p<pCatCat; p+=2){
      size = numLevels[catcatIndices[p]-1] * numLevels[catcatIndices[p+1]-1];
      for (i=0; i<size; i++){
	result[offset + i] /= factor;
      }
      offset += size;
    }
  }
  if (pContCont > 0){
    factor = sqrt(3);
    double meanZ, normZ, meanProduct, normProduct;
    double *restrict wOffsetPtr;
    double *restrict product = malloc(n * sizeof *product);
    for (p=0; p<pContCont; p+=2){
      wOffsetPtr = z + (contcontIndices[p]-1)*n;
      zOffsetPtr = z + (contcontIndices[p+1]-1)*n;
      mean = norm = meanZ = normZ = 0.0;
      for (i=0; i<n; i++){
	mean += wOffsetPtr[i];
	norm += wOffsetPtr[i]*wOffsetPtr[i];
	meanZ += zOffsetPtr[i];
	normZ += zOffsetPtr[i]*zOffsetPtr[i];
      }
      mean /= n;
      meanZ /= n;
      norm = (fabs(norm) > 1e-30) ? sqrt(norm - n*pow(mean, 2)) : 1.0;
      normZ = (fabs(normZ) > 1e-30) ? sqrt(normZ - n*pow(meanZ, 2)) : 1.0;
      result[offset] /= (factor * norm);
      result[offset+1] /= (factor * normZ);
      result[0] -= mean*result[offset] + meanZ*result[offset+1];
      meanProduct = normProduct = 0.0;
      for (i=0; i<n; i++){
	product[i] = (wOffsetPtr[i]-mean) * (zOffsetPtr[i]-meanZ) / (norm*normZ);
	meanProduct += product[i];
	normProduct += product[i]*product[i];
      }
      meanProduct /= n;
      normProduct = (fabs(normProduct) > 1e-30) ? sqrt(normProduct - n*pow(meanProduct, 2)) : 1.0;
      result[offset+2] /= (factor * normProduct);
      result[0] -= meanProduct * result[offset+2];
      result[offset+2] /= norm * normZ;
      result[offset] -= meanZ * result[offset+2];
      result[offset+1] -= mean * result[offset+2];
      result[0] += mean * meanZ * result[offset+2];
      offset += 3;
    }
    free(product);
  }
  if (pCatCont > 0){
    factor = sqrt(2);
    double factor1 = sqrt(2*n);
    for (p=0; p<pCatCont; p+=2){
      zOffsetPtr = z + (catcontIndices[p+1]-1)*n;
      size = numLevels[catcontIndices[p]-1];
      mean = norm = 0.0;
      for (i=0; i<n; i++){
	mean += zOffsetPtr[i];
	norm += zOffsetPtr[i]*zOffsetPtr[i];
      }
      mean /= n;
      norm = (fabs(norm) > 1e-30) ? sqrt(norm - n*pow(mean, 2)) : 1.0;
      for (i=0; i<size; i++){
	result[offset+size+i] /= (factor*norm);
	result[offset+i] = result[offset+i]/factor1 - mean*result[offset+size+i];
      }
      offset += 2*size;
    }
  }
}

SEXP R_rescale_beta(SEXP R_x, SEXP R_z, SEXP R_nRows, SEXP R_beta, SEXP R_betaLen, SEXP R_nVars, SEXP R_numLevels, SEXP R_catIndices, SEXP R_contIndices, SEXP R_catcatIndices, SEXP R_contcontIndices, SEXP R_catcontIndices, SEXP R_result){
  PROTECT(R_x = coerceVector(R_x, INTSXP));
  PROTECT(R_z = coerceVector(R_z, REALSXP));
  PROTECT(R_nRows = coerceVector(R_nRows, INTSXP));
  PROTECT(R_beta = coerceVector(R_beta, REALSXP));
  PROTECT(R_betaLen = coerceVector(R_betaLen, INTSXP));
  PROTECT(R_nVars = coerceVector(R_nVars, INTSXP));
  PROTECT(R_numLevels = coerceVector(R_numLevels, INTSXP));
  PROTECT(R_catIndices = coerceVector(R_catIndices, INTSXP));
  PROTECT(R_contIndices = coerceVector(R_contIndices, INTSXP));
  PROTECT(R_catcatIndices = coerceVector(R_catcatIndices, INTSXP));
  PROTECT(R_contcontIndices = coerceVector(R_contcontIndices, INTSXP));
  PROTECT(R_catcontIndices = coerceVector(R_catcontIndices, INTSXP));
  PROTECT(R_result = coerceVector(R_result, REALSXP));
  int *restrict x = INTEGER(R_x);
  double *restrict z = REAL(R_z);
  int *restrict nRows = INTEGER(R_nRows);
  double *restrict beta = REAL(R_beta);
  int *restrict betaLen = INTEGER(R_betaLen);
  int *restrict nVars = INTEGER(R_nVars);
  int *restrict numLevels = INTEGER(R_numLevels);
  int *restrict catIndices = INTEGER(R_catIndices);
  int *restrict contIndices = INTEGER(R_contIndices);
  int *restrict catcatIndices = INTEGER(R_catcatIndices);
  int *restrict contcontIndices = INTEGER(R_contcontIndices);
  int *restrict catcontIndices = INTEGER(R_catcontIndices);
  double *restrict result = REAL(R_result);
  rescale_beta(x, z, nRows, beta, betaLen, nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, result);
  UNPROTECT(13);
  return R_result;
}

void x_times_rescaled_beta(int *restrict x, double *restrict z, const double *restrict beta, const int *restrict nRows, const int *restrict nVars, const int *restrict numLevels, const int *restrict catIndices, const int *restrict contIndices, const int *restrict catcatIndices, const int *restrict contcontIndices, const int *restrict catcontIndices, double *restrict result){
  int i, p, size, n = *nRows, offset = 1;
  int pCat = nVars[0], pCont = nVars[1], pCatCat = 2*nVars[2], pContCont = 2*nVars[3], pCatCont = 2*nVars[4];
  int *restrict xOffsetPtr;
  double *restrict zOffsetPtr;
  for (i=0; i<n; i++){
    result[i] = beta[0];
  }
  if (pCat > 0){
    for (p=0; p<pCat; p++){
      size = numLevels[catIndices[p]-1];
      xOffsetPtr = x + (catIndices[p]-1)*n;
      for (i=0; i<n; i++){
	result[i] += beta[offset + xOffsetPtr[i]];
      }
      offset += size;
    }
  }
  if (pCont > 0){
    for (p=0; p<pCont; p++){
      zOffsetPtr = z + (contIndices[p]-1)*n;
      for (i=0; i<n; i++){
	result[i] += zOffsetPtr[i] * beta[offset];
      }
      ++offset;
    }
  }
  if (pCatCat > 0){
    int *restrict yOffsetPtr;
    for (p=0; p<pCatCat; p+=2){
      size = numLevels[catcatIndices[p]-1];
      xOffsetPtr = x + (catcatIndices[p]-1)*n;
      yOffsetPtr = x + (catcatIndices[p+1]-1)*n;
      for (i=0; i<n; i++){
	result[i] += beta[offset + xOffsetPtr[i] + size*yOffsetPtr[i]];
      }
      offset += size * numLevels[catcatIndices[p+1]-1];
    }
  }
  if (pContCont > 0){
    double *restrict wOffsetPtr;
    for (p=0; p<pContCont; p+=2){
      wOffsetPtr = z + (contcontIndices[p]-1)*n;
      zOffsetPtr = z + (contcontIndices[p+1]-1)*n;
      for (i=0; i<n; i++){
	result[i] += wOffsetPtr[i]*beta[offset] + zOffsetPtr[i]*beta[offset+1] + wOffsetPtr[i]*zOffsetPtr[i]*beta[offset+2];
      }
      offset += 3;
    }
  }
  if (pCatCont > 0){
    for (p=0; p<pCatCont; p+=2){
      size = numLevels[catcontIndices[p]-1];
      xOffsetPtr = x + (catcontIndices[p]-1)*n;
      zOffsetPtr = z + (catcontIndices[p+1]-1)*n;
      for (i=0; i<n; i++){
	result[i] += beta[offset + xOffsetPtr[i]] + zOffsetPtr[i]*beta[offset + size + xOffsetPtr[i]];
      }
      offset += 2*size;
    }
  }
}	

SEXP R_x_times_rescaled_beta(SEXP R_x, SEXP R_z, SEXP R_beta, SEXP R_nRows, SEXP R_nVars, SEXP R_numLevels, SEXP R_catIndices, SEXP R_contIndices, SEXP R_catcatIndices, SEXP R_contcontIndices, SEXP R_catcontIndices, SEXP R_result){
  PROTECT(R_x = coerceVector(R_x, INTSXP));
  PROTECT(R_z = coerceVector(R_z, REALSXP));
  PROTECT(R_nRows = coerceVector(R_nRows, INTSXP));
  PROTECT(R_beta = coerceVector(R_beta, REALSXP));
  PROTECT(R_nVars = coerceVector(R_nVars, INTSXP));
  PROTECT(R_numLevels = coerceVector(R_numLevels, INTSXP));
  PROTECT(R_catIndices = coerceVector(R_catIndices, INTSXP));
  PROTECT(R_contIndices = coerceVector(R_contIndices, INTSXP));
  PROTECT(R_catcatIndices = coerceVector(R_catcatIndices, INTSXP));
  PROTECT(R_contcontIndices = coerceVector(R_contcontIndices, INTSXP));
  PROTECT(R_catcontIndices = coerceVector(R_catcontIndices, INTSXP));
  PROTECT(R_result = coerceVector(R_result, REALSXP));
  int *restrict x = INTEGER(R_x);
  double *restrict z = REAL(R_z);
  int *restrict nRows = INTEGER(R_nRows);
  double *restrict beta = REAL(R_beta);
  int *restrict nVars = INTEGER(R_nVars);
  int *restrict numLevels = INTEGER(R_numLevels);
  int *restrict catIndices = INTEGER(R_catIndices);
  int *restrict contIndices = INTEGER(R_contIndices);
  int *restrict catcatIndices = INTEGER(R_catcatIndices);
  int *restrict contcontIndices = INTEGER(R_contcontIndices);
  int *restrict catcontIndices = INTEGER(R_catcontIndices);
  double *restrict result = REAL(R_result);
  x_times_rescaled_beta(x, z, beta, nRows, nVars, numLevels, catIndices, contIndices, catcatIndices, contcontIndices, catcontIndices, result);
  UNPROTECT(12);
  return R_result;
}

void compute_norms_cat(int *restrict x, double *restrict r, int *restrict nRows, int *restrict nVars, int *restrict numLevels, int *restrict numCores, double *restrict result){
  int n = *nRows;
  int p = *nVars;
  int i, j, offset, len;
  double *restrict temp;
#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(*numCores);
#endif
#pragma pomp inst begin(crout_compute_norms_cat)
#pragma omp parallel for shared(x, r, n, p, result) private(i, j, offset, len, temp)
  for (j=0; j<p; j++){
    offset = j*n;
    len = numLevels[j];
    temp = calloc(len, sizeof *temp);
    for (i=0; i<n; i++){
      temp[x[offset+i]] += r[i];
    }
    for (i=0; i<len; i++){
      result[j] += temp[i]*temp[i];
    }
    result[j] = sqrt(result[j]/n)/n;
    free(temp);
  }
#pragma pomp inst end(crout_compute_norms_cat)
}

SEXP R_compute_norms_cat(SEXP R_x, SEXP R_r, SEXP R_nRows, SEXP R_nVars, SEXP R_numLevels, SEXP R_numCores, SEXP R_result){
  PROTECT(R_x = coerceVector(R_x, INTSXP));
  PROTECT(R_r = coerceVector(R_r, REALSXP));
  PROTECT(R_nRows = coerceVector(R_nRows, INTSXP));
  PROTECT(R_nVars = coerceVector(R_nVars, INTSXP));
  PROTECT(R_numLevels = coerceVector(R_numLevels, INTSXP));
  PROTECT(R_numCores = coerceVector(R_numCores, INTSXP));
  PROTECT(R_result = coerceVector(R_result, REALSXP));
  int *restrict x = INTEGER(R_x);
  double *restrict r = REAL(R_r);
  int *restrict nRows = INTEGER(R_nRows);
  int *restrict nVars = INTEGER(R_nVars);
  int *restrict numLevels = INTEGER(R_numLevels);
  int *restrict numCores = INTEGER(R_numCores);
  double *restrict result = REAL(R_result);
  compute_norms_cat(x, r, nRows, nVars, numLevels, numCores, result);
  UNPROTECT(7);
  return R_result;
}

void compute_norms_cat_cat(int *restrict x, double *restrict r, int *restrict nRows, int *restrict nVars, int *restrict numLevels, int *restrict xIndices, int *restrict yIndices, int *restrict numCores, double *restrict result){
  int n = *nRows;
  int p = *nVars;
  int i, j, xOffset, yOffset, len, xlevels;
  double *restrict temp;
#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(*numCores);
#endif
#pragma pomp inst begin(crout_compute_norms_cat_cat)
#pragma omp parallel for shared(x, r, n, p, numLevels, xIndices, yIndices, result) private(i, j, xOffset, yOffset, len, xlevels, temp)
  for (j=0; j<p; j++){
    xOffset = (xIndices[j] - 1)*n;  //R uses 1-based indexing
    yOffset = (yIndices[j] - 1)*n;
    xlevels = numLevels[xIndices[j]-1];
    len = xlevels * numLevels[yIndices[j]-1];
    temp = calloc(len, sizeof *temp);
    for (i=0; i<n; i++){
      temp[x[xOffset+i] + xlevels*x[yOffset+i]] += r[i];
    }
    for (i=0; i<len; i++){
      result[j] += temp[i]*temp[i];
    }
    result[j] = sqrt(result[j]/n)/n;
    free(temp);
  }
#pragma pomp inst end(crout_compute_norms_cat_cat)
}

SEXP R_compute_norms_cat_cat(SEXP R_x, SEXP R_r, SEXP R_nRows, SEXP R_nVars, SEXP R_numLevels, SEXP R_xIndices, SEXP R_yIndices, SEXP R_numCores, SEXP R_result){
  PROTECT(R_x = coerceVector(R_x, INTSXP));
  PROTECT(R_r = coerceVector(R_r, REALSXP));
  PROTECT(R_nRows = coerceVector(R_nRows, INTSXP));
  PROTECT(R_nVars = coerceVector(R_nVars, INTSXP));
  PROTECT(R_numLevels = coerceVector(R_numLevels, INTSXP));
  PROTECT(R_xIndices = coerceVector(R_xIndices, INTSXP));
  PROTECT(R_yIndices = coerceVector(R_yIndices, INTSXP));
  PROTECT(R_numCores = coerceVector(R_numCores, INTSXP));
  PROTECT(R_result = coerceVector(R_result, REALSXP));
  int *restrict x = INTEGER(R_x);
  double *restrict r = REAL(R_r);
  int *restrict nRows = INTEGER(R_nRows);
  int *restrict nVars = INTEGER(R_nVars);
  int *restrict numLevels = INTEGER(R_numLevels);
  int *restrict xIndices = INTEGER(R_xIndices);
  int *restrict yIndices = INTEGER(R_yIndices);
  int *restrict numCores = INTEGER(R_numCores);
  double *restrict result = REAL(R_result);
  compute_norms_cat_cat(x, r, nRows, nVars, numLevels, xIndices, yIndices, numCores, result);
  UNPROTECT(9);
  return R_result;
}

void compute_norms_cat_cont(int *restrict x, double *restrict z, double *restrict catNorms, double *restrict r, int *restrict nRows, int *restrict nVars, int *restrict numLevels, int *restrict xIndices, int *restrict zIndices, int *restrict numCores, double *restrict result){
  int n = *nRows;
  int p = *nVars;
  int i, j, xOffset, zOffset, levels;
  double *restrict temp;
#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(*numCores);
#endif
#pragma pomp inst begin(crout_compute_norms_cat_cont)
#pragma omp parallel for shared(x, z, catNorms, r, n, p, numLevels, xIndices, zIndices, result) private(i, j, xOffset, zOffset, levels, temp)
  for (j=0; j<p; j++){
    xOffset = (xIndices[j] - 1)*n;
    zOffset = (zIndices[j] - 1)*n;
    levels = numLevels[xIndices[j]-1];
    temp = calloc(levels, sizeof *temp);
    for (i=0; i<n; i++){
      temp[x[xOffset+i]] += z[zOffset+i]*r[i];
    }
    result[j] = pow(n*catNorms[xIndices[j]-1], 2);
    for (i=0; i<levels; i++){
      result[j] += temp[i]*temp[i];
    }
    result[j] = sqrt(result[j]/2)/n;
    free(temp);
  }
#pragma pomp inst end(crout_compute_norms_cat_cont)
}

SEXP R_compute_norms_cat_cont(SEXP R_x, SEXP R_z, SEXP R_catNorms, SEXP R_r, SEXP R_nRows, SEXP R_nVars, SEXP R_numLevels, SEXP R_xIndices, SEXP R_zIndices, SEXP R_numCores, SEXP R_result){
  PROTECT(R_x = coerceVector(R_x, INTSXP));
  PROTECT(R_z = coerceVector(R_z, REALSXP));
  PROTECT(R_catNorms = coerceVector(R_catNorms, REALSXP));
  PROTECT(R_r = coerceVector(R_r, REALSXP));
  PROTECT(R_nRows = coerceVector(R_nRows, INTSXP));
  PROTECT(R_nVars = coerceVector(R_nVars, INTSXP));
  PROTECT(R_numLevels = coerceVector(R_numLevels, INTSXP));
  PROTECT(R_xIndices = coerceVector(R_xIndices, INTSXP));
  PROTECT(R_zIndices = coerceVector(R_zIndices, INTSXP));
  PROTECT(R_numCores = coerceVector(R_numCores, INTSXP));
  PROTECT(R_result = coerceVector(R_result, REALSXP));
  int *restrict x = INTEGER(R_x);
  double *restrict z = REAL(R_z);
  double *restrict catNorms = REAL(R_catNorms);
  double *restrict r = REAL(R_r);
  int *restrict nRows = INTEGER(R_nRows);
  int *restrict nVars = INTEGER(R_nVars);
  int *restrict numLevels = INTEGER(R_numLevels);
  int *restrict xIndices = INTEGER(R_xIndices);
  int *restrict zIndices = INTEGER(R_zIndices);
  int *restrict numCores = INTEGER(R_numCores);
  double *restrict result = REAL(R_result);
  compute_norms_cat_cont(x, z, catNorms, r, nRows, nVars, numLevels, xIndices, zIndices, numCores, result);
  UNPROTECT(11);
  return R_result;
}

void compute_norms_cont_cont(double *restrict x, double *restrict contNorms, double *restrict r, int *restrict nRows, int *restrict nVars, int *restrict xIndices, int *restrict yIndices, int *restrict numCores, double *restrict result){
  int n = *nRows;
  int p = *nVars;
  int i, j;
  size_t xOffset, yOffset;
  double mean, norm, temp, rprd, rsum;
  double product;
  long n2 = (long)n * n;
#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(*numCores);
#endif
#pragma pomp inst begin(crout_compute_norms_cont_cont)
#ifdef __AVX__
  if(max_alignment(x) < 64) {
    Rf_error("alignment of x is %d; need at least 64bit alignment", max_alignment(x));
  }
  if(max_alignment(r) < 64) {
    Rf_error("alignment of r is %d; need at least 64bit alignment", max_alignment(r));
  }
  int nRowsDiv8 = n/8;
#pragma omp parallel for shared(x, contNorms, r, n, p, xIndices, yIndices, result) private(i, j, xOffset, yOffset, mean, norm, temp, product)
  for (i=0; i<p; ++i) {
    xOffset = ((size_t)xIndices[i] - 1)*n;
    yOffset = ((size_t)yIndices[i] - 1)*n;
    __m256d mean1_pd = _mm256_set_pd(0.0,0.0,0.0,0.0);
    __m256d mean2_pd = _mm256_set_pd(0.0,0.0,0.0,0.0);
    __m256d norm1_pd = _mm256_set_pd(0.0,0.0,0.0,0.0);
    __m256d norm2_pd = _mm256_set_pd(0.0,0.0,0.0,0.0);
    __m256d rprd1_pd = _mm256_set_pd(0.0,0.0,0.0,0.0);
    __m256d rprd2_pd = _mm256_set_pd(0.0,0.0,0.0,0.0);
    __m256d rsum1_pd = _mm256_set_pd(0.0,0.0,0.0,0.0);
    __m256d rsum2_pd = _mm256_set_pd(0.0,0.0,0.0,0.0);
    for (j=0; j<nRowsDiv8; ++j) {
      __m256d xx1_pd = _mm256_load_pd(x+xOffset+j*8);
      __m256d xx2_pd = _mm256_load_pd(x+xOffset+j*8+4);
      __m256d xy1_pd = _mm256_load_pd(x+yOffset+j*8);
      __m256d xy2_pd = _mm256_load_pd(x+yOffset+j*8+4);
      __m256d r1_pd = _mm256_load_pd(r+j*8);
      __m256d r2_pd = _mm256_load_pd(r+j*8+4);
      __m256d prod1_pd = _mm256_mul_pd(xx1_pd, xy1_pd);
      __m256d prod2_pd = _mm256_mul_pd(xx2_pd, xy2_pd);
      __m256d sqr1_pd = _mm256_mul_pd(prod1_pd, prod1_pd);
      __m256d sqr2_pd = _mm256_mul_pd(prod2_pd, prod2_pd);
      __m256d prdr1_pd = _mm256_mul_pd(r1_pd, prod1_pd);
      __m256d prdr2_pd = _mm256_mul_pd(r2_pd, prod2_pd);
      mean1_pd = _mm256_add_pd(mean1_pd,prod1_pd);
      mean2_pd = _mm256_add_pd(mean2_pd,prod2_pd);
      norm1_pd = _mm256_add_pd(norm1_pd,sqr1_pd);
      norm2_pd = _mm256_add_pd(norm2_pd,sqr2_pd);
      rprd1_pd = _mm256_add_pd(rprd1_pd,prdr1_pd);
      rprd2_pd = _mm256_add_pd(rprd2_pd,prdr2_pd);
      rsum1_pd = _mm256_add_pd(rsum1_pd,r1_pd);
      rsum2_pd = _mm256_add_pd(rsum2_pd,r2_pd);
    }
    __m256d t_mean = _mm256_add_pd(mean1_pd, mean2_pd);
    __m256d h_mean = _mm256_add_pd(t_mean, _mm256_permute2f128_pd(t_mean, t_mean, 0x1));
    _mm_store_sd(&mean, _mm_hadd_pd( _mm256_castpd256_pd128(h_mean), _mm256_castpd256_pd128(h_mean)));
    __m256d t_norm = _mm256_add_pd(norm1_pd, norm2_pd);
    __m256d h_norm = _mm256_add_pd(t_norm, _mm256_permute2f128_pd(t_norm, t_norm, 0x1));
    _mm_store_sd(&norm, _mm_hadd_pd( _mm256_castpd256_pd128(h_norm), _mm256_castpd256_pd128(h_norm)));
    __m256d t_rprd = _mm256_add_pd(rprd1_pd, rprd2_pd);
    __m256d h_rprd = _mm256_add_pd(t_rprd, _mm256_permute2f128_pd(t_rprd, t_rprd, 0x1));
    _mm_store_sd(&rprd, _mm_hadd_pd( _mm256_castpd256_pd128(h_rprd), _mm256_castpd256_pd128(h_rprd)));
    __m256d t_rsum = _mm256_add_pd(rsum1_pd, rsum2_pd);
    __m256d h_rsum = _mm256_add_pd(t_rsum, _mm256_permute2f128_pd(t_rsum, t_rsum, 0x1));
    _mm_store_sd(&rsum, _mm_hadd_pd( _mm256_castpd256_pd128(h_rsum), _mm256_castpd256_pd128(h_rsum)));
    
    for (j=nRowsDiv8*8; j<n; ++j){
      product = x[xOffset+j]*x[yOffset+j];
      mean += product;
      norm += product*product;
      rprd += r[j] * product;
      rsum += r[j];
    }
    mean /= n;
    temp = rprd - mean * rsum;
    result[i] += n2 
                  * ((contNorms[xIndices[i]-1] * contNorms[xIndices[i]-1]) 
                    + (contNorms[yIndices[i]-1] * contNorms[yIndices[i]-1])) 
                  + (norm > 0 ? (temp * temp)/(norm-n*(mean * mean)) : 0);
    result[i] = sqrt(result[i]/3)/n;
  }
#else
  if(max_alignment(x) < 32) {
    Rf_error("alignment of x is %d; need at least 32bit alignment", max_alignment(x));
  }
  if(max_alignment(r) < 32) {
    Rf_error("alignment of r is %d; need at least 32bit alignment", max_alignment(r));
  }
  int nRowsDiv2 = n/2;
#pragma omp parallel for shared(x, contNorms, r, n, p, xIndices, yIndices, result) private(i, j, xOffset, yOffset, mean, norm, temp, product)
  for (i=0; i<p; ++i){
    xOffset = ((size_t)xIndices[i] - 1)*n;
    yOffset = ((size_t)yIndices[i] - 1)*n;
    __m128d mean_pd = _mm_set_pd(0.0,0.0);
    __m128d norm_pd = _mm_set_pd(0.0,0.0);
    __m128d rprd_pd = _mm_set_pd(0.0,0.0);
    __m128d rsum_pd = _mm_set_pd(0.0,0.0);
    for (j=0; j<nRowsDiv2; ++j){
      __m128d xx_pd = _mm_load_pd(x+xOffset+j*2);
      __m128d xy_pd = _mm_load_pd(x+yOffset+j*2);
      __m128d r_pd  = _mm_load_pd(r+j*2);
      __m128d prod_pd = _mm_mul_pd(xx_pd, xy_pd);
      __m128d sqr_pd  = _mm_mul_pd(prod_pd, prod_pd);
      __m128d prdr_pd = _mm_mul_pd(r_pd, prod_pd);
      mean_pd = _mm_add_pd(mean_pd, prod_pd);
      norm_pd = _mm_add_pd(norm_pd, sqr_pd);
      rprd_pd = _mm_add_pd(rprd_pd, prdr_pd);
      rsum_pd = _mm_add_pd(rsum_pd, r_pd);
    }
    __m128d mennrm = _mm_hadd_pd(mean_pd, norm_pd);
    __m128d prdsum = _mm_hadd_pd(rprd_pd, rsum_pd);
    _mm_storel_pd(&mean, mennrm);
    _mm_storeh_pd(&norm, mennrm);
    _mm_storel_pd(&rprd, prdsum);
    _mm_storeh_pd(&rsum, prdsum);
    
    for (j=nRowsDiv2*2; j<n; ++j){
      product = x[xOffset+j]*x[yOffset+j];
      mean += product;
      norm += product*product;
      rprd += r[j] * product;
      rsum += r[j];
    }

    mean /= n;
    temp = rprd - mean * rsum;
    
    result[i] += n2 
                  * ((contNorms[xIndices[i]-1] * contNorms[xIndices[i]-1]) 
                    + (contNorms[yIndices[i]-1] * contNorms[yIndices[i]-1])) 
                  + (norm > 0 ? (temp * temp)/(norm-n*(mean * mean)) : 0);
    result[i] = sqrt(result[i]/3)/n;
  }
#endif
#pragma pomp inst end(crout_compute_norms_cont_cont)
}

SEXP R_compute_norms_cont_cont(SEXP R_x, SEXP R_contNorms, SEXP R_r, SEXP R_nRows, SEXP R_nVars, SEXP R_xIndices, SEXP R_yIndices, Rboolean verbose, SEXP R_numCores, SEXP R_result){
  PROTECT(R_x = coerceVector(R_x, REALSXP));
  PROTECT(R_contNorms = coerceVector(R_contNorms, REALSXP));
  PROTECT(R_r = coerceVector(R_r, REALSXP));
  PROTECT(R_nRows = coerceVector(R_nRows, INTSXP));
  PROTECT(R_nVars = coerceVector(R_nVars, INTSXP));
  PROTECT(R_xIndices = coerceVector(R_xIndices, INTSXP));
  PROTECT(R_yIndices = coerceVector(R_yIndices, INTSXP));
  PROTECT(R_numCores = coerceVector(R_numCores, INTSXP));
  PROTECT(R_result = coerceVector(R_result, REALSXP));
  double *restrict x = REAL(R_x);
  double *restrict contNorms = REAL(R_contNorms);
  double *restrict r = REAL(R_r);
  int *restrict nRows = INTEGER(R_nRows);
  int *restrict nVars = INTEGER(R_nVars);
  int *restrict xIndices = INTEGER(R_xIndices);
  int *restrict yIndices = INTEGER(R_yIndices);
  int *restrict numCores = INTEGER(R_numCores);
  double *restrict result = REAL(R_result);
  struct timespec timer_norms_cont_cont;
  if (verbose) {  
    timer_norms_cont_cont = timer_start();
  }
  compute_norms_cont_cont(x, contNorms, r, nRows, nVars, xIndices, yIndices, numCores, result);
  if (verbose) {
    Rprintf("---> timer (compute_norms_cont_cont): %lf [s]\n", timer_end(timer_norms_cont_cont));
  }
  UNPROTECT(9);
  return R_result;
}
