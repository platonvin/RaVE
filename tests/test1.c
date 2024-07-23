#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>
#include <pthread.h>
// #include <stdio.h>
// #include <stdlib.h>
#include "../src/count_cores.c"

static unsigned long int next = 1;
// static __m128i next_simd;
int rand_simd(void)  /* RAND_MAX assumed to be 32767. */
{
    // int a = 1103515245;
    next = next * 1103515245 + 12345;
    return (unsigned)(next);
    // return (unsigned)(next) % (1 << 16) - (1 << 15);
}

void printBitsColored(int num) {
    unsigned int mask = 1 << (sizeof(int) * 8 - 1); // Mask to isolate the highest bit
    printf("Bits: ");
    
    for (int i = 0; i < sizeof(int) * 8; ++i) {
        if (num & mask) {
            printf("\033[32m1\033[0m"); // Green for "1" bit
        } else {
            printf("\033[31m0\033[0m"); // Red for "0" bit
        }
        mask >>= 1;
    }
    printf("\n");
}

void print_int128(__int128 value) {
    char buffer[40]; // Buffer to hold the string representation
    int index = 39;  // Start at the end of the buffer
    buffer[index] = '\0'; // Null-terminate the string

    int is_negative = 0;
    if (value < 0) {
        is_negative = 1;
        value = -value;
    }

    // Convert the number to a string
    do {
        buffer[--index] = '0' + (value % 10);
        value /= 10;
    } while (value != 0);

    if (is_negative) {
        buffer[--index] = '-';
    }

    // Print the string
    printf("%s\n", &buffer[index]);
}

int main(){
    // int n = sysconf(_SC_NPROCESSORS_ONLN);
    int n = get_core_count();

    printf("%d", n);
}