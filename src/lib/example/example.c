#include "example.h"

#include <math.h>

double log_binomial_coeff(unsigned n, unsigned k) {
    if (n < k)
        return 0;
    return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}
