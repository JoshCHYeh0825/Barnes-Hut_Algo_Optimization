#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Fast inverse square root from Quake III Arena
float fast_inverse_sqrt(float number) {
    long i;
    float x2, y;
    const float threehalfs = 1.5F;

    x2 = number * 0.5F;
    y = number;
    i = *(long *)&y;           // bit-level hack
    i = 0x5f3759df - (i >> 1); // magic number
    y = *(float *)&i;
    y = y * (threehalfs - (x2 * y * y)); // 1st Newton iteration
    // y = y * (threehalfs - (x2 * y * y));  // 2nd iteration for more accuracy

    return y;
}

// Method 1: Using powf directly
float method_powf(float d_sq, float e_sq) {
    return powf(d_sq + e_sq, 1.5f);
}

// Method 2: Manual calculation (base * sqrt(base))
float method_manual(float d_sq, float e_sq) {
    float base = d_sq + e_sq;
    return base * sqrtf(base);
}

// Method 3: Using fast inverse square root
float method_fast_invsqrt(float d_sq, float e_sq) {
    float base = d_sq + e_sq;
    float invSqrt = fast_inverse_sqrt(base);
    return base * base * invSqrt;
}

int main() {
    const int NUM_TESTS = 10000000; // Number of iterations for benchmarking
    const int NUM_VALUES = 1000;    // Number of different distance values to test
    float *test_values = malloc(NUM_VALUES * sizeof(float));

    // Generate some random test values
    srand(time(NULL));
    for (int i = 0; i < NUM_VALUES; i++) {
        // Generate values that might be typical in a simulation (between 0.1 and 100)
        test_values[i] = 0.1f + (float)rand() / RAND_MAX * 99.9f;
    }

    // Small epsilon squared value typical for softening
    const float e_sq = 0.01f;

    clock_t start, end;
    double cpu_time_used;
    float result = 0.0f; // To prevent compiler from optimizing away calculations

    // Test method 1: powf
    start = clock();
    for (int i = 0; i < NUM_TESTS; i++) {
        float d_sq = test_values[i % NUM_VALUES];
        result += method_powf(d_sq, e_sq);
    }
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Method 1 (powf):              %f seconds\n", cpu_time_used);
    printf("  Result checksum: %f\n", result);

    // Reset result for next test
    result = 0.0f;

    // Test method 2: manual calculation
    start = clock();
    for (int i = 0; i < NUM_TESTS; i++) {
        float d_sq = test_values[i % NUM_VALUES];
        result += method_manual(d_sq, e_sq);
    }
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Method 2 (manual calculation): %f seconds\n", cpu_time_used);
    printf("  Result checksum: %f\n", result);

    // Reset result for next test
    result = 0.0f;

    // Test method 3: fast inverse square root
    start = clock();
    for (int i = 0; i < NUM_TESTS; i++) {
        float d_sq = test_values[i % NUM_VALUES];
        result += method_fast_invsqrt(d_sq, e_sq);
    }
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Method 3 (fast inverse sqrt):  %f seconds\n", cpu_time_used);
    printf("  Result checksum: %f\n", result);

    // Compare accuracy
    printf("\nAccuracy comparison (sample values):\n");
    printf("%-10s %-15s %-15s %-15s %-15s\n", "d_sq", "powf", "manual", "fast_invsqrt", "max_diff");

    for (int i = 0; i < 5; i++) {
        float d_sq = test_values[i];
        float res1 = method_powf(d_sq, e_sq);
        float res2 = method_manual(d_sq, e_sq);
        float res3 = method_fast_invsqrt(d_sq, e_sq);
        float diff1 = fabsf(res1 - res2);
        float diff2 = fabsf(res1 - res3);
        float max_diff = fmaxf(diff1, diff2);

        printf("%-10.4f %-15.6f %-15.6f %-15.6f %-15.6f\n", d_sq, res1, res2, res3, max_diff);
    }

    free(test_values);
    return 0;
}