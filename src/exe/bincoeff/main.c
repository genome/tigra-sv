#include "example/example.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv) {
    unsigned n = 0;
    unsigned k = 0;
    char* endp = NULL;
    double r = 0.0;

    if (argc != 3) {
        fprintf(stderr, "usage: %s <n> <k>\n\n", argv[0]);
        fprintf(stderr, "    print binomial coefficient via lgamma function\n");
        return 1;
    }


    n = strtoul(argv[1], &endp, 10);
    if (endp-argv[1] != strlen(argv[1])) {
        fprintf(stderr, "Couldn't parse %s as a number\n", argv[1]);
        return 1;
    }

    k = strtoul(argv[2], &endp, 10);
    if (endp-argv[2] != strlen(argv[2])) {
        fprintf(stderr, "Couldn't parse %s as a number\n", argv[2]);
        return 1;
    }

    r = log_binomial_coeff(n, k);
    printf("log(%u choose %u): %f\n", n, k, r);
    printf("%u choose %u: %lu\n", n, k, (unsigned long)(exp(r) + 0.499));
    
    return 0;
}
