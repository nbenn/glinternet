#ifndef UTIL_HEADER
#define UTIL_HEADER

#define _POSIX_C_SOURCE 199309L

#include <time.h>

struct timespec timer_start();
double timer_end(struct timespec);

#endif