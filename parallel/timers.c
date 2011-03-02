#define _REENTRANT
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <time.h>               /* clock */
#include <string.h>             /* strrchr */
#include <math.h>               /* sqrt, M_PI, cos, sin */
#include <float.h>              /* DBL_EPSILON */


#include "timers.h"

void initialize_timer(Timer * t)
{
    t->clock_holder = 0;
    t->duration_clocks = 0;
}

void start_timer(Timer * t)
{
    t->clock_holder = clock();
}

void stop_timer(Timer * t)
{
    clock_t end_clock = clock();

    t->duration_clocks += (end_clock - t->clock_holder);
    assert(t->duration_clocks >= 0);
}

double timer_duration(const Timer * t)
{
    return ((double) t->duration_clocks) / CLOCKS_PER_SEC;
}
