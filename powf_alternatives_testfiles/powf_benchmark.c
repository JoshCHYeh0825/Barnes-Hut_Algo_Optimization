#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Number of iterations for benchmarking
#define NUM_ITERATIONS 100000000
#define NUM_TEST_VALUES 10000

// Fast inverse square root (Quake III technique)
float fast_inv_sqrt(float x) {
    float xhalf = 0.5f * x;
    int i = *(int *)&x;
    i = 0x5f3759df - (i >> 1);
    x = *(float *)&i;
    x = x * (1.5f - xhalf * x * x); // One iteration for more accuracy
    return x;
}

// Original powf implementation
float powf_original(float base, float exponent) {
    return powf(base, exponent);
}

// Method 1: Direct multiplication
float powf_direct(float base, float exponent) {
    if (exponent != 1.5f) {
        fprintf(stderr, "powf_direct only supports exponent 1.5\n");
        return 0.0f;
    }
    return base * sqrtf(base);
}

// Method 2: Using separate square root
float powf_sqrt_reuse(float base, float exponent) {
    if (exponent != 1.5f) {
        fprintf(stderr, "powf_sqrt_reuse only supports exponent 1.5\n");
        return 0.0f;
    }
    float sqrt_val = sqrtf(base);
    return base * sqrt_val;
}

// Method 3: Using fast inverse square root
float powf_fast_inv_sqrt(float base, float exponent) {
    if (exponent != 1.5f) {
        fprintf(stderr, "powf_fast_inv_sqrt only supports exponent 1.5\n");
        return 0.0f;
    }
    float inv_sqrt = fast_inv_sqrt(base);
    return base * (1.0f / inv_sqrt);
}

// For testing accuracy
typedef struct {
    float input;
    float result_original;
    float result_direct;
    float result_sqrt_reuse;
    float result_fast_inv_sqrt;
    float error_direct;
    float error_sqrt_reuse;
    float error_fast_inv_sqrt;
} TestResult;

double get_time_ms() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1000000.0;
}

void print_benchmark_results(const char *method_name, double time_ms) {
    printf("%-25s %10.3f ms\n", method_name, time_ms);
}

int main() {
    printf("===== POWF OPTIMIZATION BENCHMARK =====\n\n");

    // Initialize random seed
    srand((unsigned int)time(NULL));

    // Generate test values (d_sq + e_sq values)
    float *test_values = (float *)malloc(NUM_TEST_VALUES * sizeof(float));
    for (int i = 0; i < NUM_TEST_VALUES; i++) {
        // Generate values that mimic realistic distance calculations (d_sq + e_sq)
        // Range from small to large values to test different scenarios
        if (i < NUM_TEST_VALUES / 3) {
            // Small values (close objects)
            test_values[i] = ((float)rand() / RAND_MAX) * 0.1f + 0.0001f;
        } else if (i < 2 * NUM_TEST_VALUES / 3) {
            // Medium values (typical distances)
            test_values[i] = ((float)rand() / RAND_MAX) * 10.0f + 0.1f;
        } else {
            // Large values (far objects)
            test_values[i] = ((float)rand() / RAND_MAX) * 1000.0f + 10.0f;
        }
    }

    // ------------ PERFORMANCE BENCHMARK ------------

    printf("Performance Benchmark (%d iterations per method):\n", NUM_ITERATIONS);
    printf("%-25s %10s\n", "Method", "Time (ms)");
    printf("--------------------------------------\n");

    // Benchmark original powf
    double start_time = get_time_ms();
    volatile float result = 0.0f; // volatile to prevent optimization
    for (int i = 0; i < NUM_ITERATIONS; i++) {
        float val = test_values[i % NUM_TEST_VALUES];
        result = powf_original(val, 1.5f);
    }
    double end_time = get_time_ms();
    print_benchmark_results("Original powf", end_time - start_time);

    // Benchmark direct multiplication
    start_time = get_time_ms();
    for (int i = 0; i < NUM_ITERATIONS; i++) {
        float val = test_values[i % NUM_TEST_VALUES];
        result = powf_direct(val, 1.5f);
    }
    end_time = get_time_ms();
    print_benchmark_results("Direct multiplication", end_time - start_time);

    // Benchmark sqrt reuse
    start_time = get_time_ms();
    for (int i = 0; i < NUM_ITERATIONS; i++) {
        float val = test_values[i % NUM_TEST_VALUES];
        result = powf_sqrt_reuse(val, 1.5f);
    }
    end_time = get_time_ms();
    print_benchmark_results("Sqrt reuse", end_time - start_time);

    // Benchmark fast inverse sqrt
    start_time = get_time_ms();
    for (int i = 0; i < NUM_ITERATIONS; i++) {
        float val = test_values[i % NUM_TEST_VALUES];
        result = powf_fast_inv_sqrt(val, 1.5f);
    }
    end_time = get_time_ms();
    print_benchmark_results("Fast inverse sqrt", end_time - start_time);

    // ------------ ACCURACY BENCHMARK ------------

    printf("\nAccuracy Benchmark:\n");
    printf("%-10s %-13s %-13s %-13s %-13s %-11s %-11s %-11s\n", "Input", "Original", "Direct", "Sqrt Reuse",
           "Fast InvSqrt", "Error Dir", "Error Sqrt", "Error Fast");
    printf("----------------------------------------------------------------------------------------\n");

    // Calculate accuracy for a subset of test values
    float max_error_direct = 0.0f;
    float max_error_sqrt_reuse = 0.0f;
    float max_error_fast_inv_sqrt = 0.0f;
    float avg_error_direct = 0.0f;
    float avg_error_sqrt_reuse = 0.0f;
    float avg_error_fast_inv_sqrt = 0.0f;

    int num_display = 10; // Number of rows to display
    int step = NUM_TEST_VALUES / num_display;

    for (int i = 0; i < NUM_TEST_VALUES; i++) {
        float val = test_values[i];
        float res_original = powf_original(val, 1.5f);
        float res_direct = powf_direct(val, 1.5f);
        float res_sqrt_reuse = powf_sqrt_reuse(val, 1.5f);
        float res_fast_inv_sqrt = powf_fast_inv_sqrt(val, 1.5f);

        float err_direct = fabsf((res_direct - res_original) / res_original);
        float err_sqrt_reuse = fabsf((res_sqrt_reuse - res_original) / res_original);
        float err_fast_inv_sqrt = fabsf((res_fast_inv_sqrt - res_original) / res_original);

        // Update max errors
        max_error_direct = fmaxf(max_error_direct, err_direct);
        max_error_sqrt_reuse = fmaxf(max_error_sqrt_reuse, err_sqrt_reuse);
        max_error_fast_inv_sqrt = fmaxf(max_error_fast_inv_sqrt, err_fast_inv_sqrt);

        // Update average errors
        avg_error_direct += err_direct;
        avg_error_sqrt_reuse += err_sqrt_reuse;
        avg_error_fast_inv_sqrt += err_fast_inv_sqrt;

        // Display some sample values
        if (i % step == 0) {
            printf("%-10.4f %-13.6f %-13.6f %-13.6f %-13.6f %-11.6f %-11.6f %-11.6f\n", val, res_original, res_direct,
                   res_sqrt_reuse, res_fast_inv_sqrt, err_direct * 100.0f, err_sqrt_reuse * 100.0f,
                   err_fast_inv_sqrt * 100.0f);
        }
    }

    // Calculate average errors
    avg_error_direct /= NUM_TEST_VALUES;
    avg_error_sqrt_reuse /= NUM_TEST_VALUES;
    avg_error_fast_inv_sqrt /= NUM_TEST_VALUES;

    printf("\nError Summary:\n");
    printf("%-20s %-15s %-15s\n", "Method", "Max Error (%)", "Avg Error (%)");
    printf("----------------------------------------------\n");
    printf("%-20s %-15.6f %-15.6f\n", "Direct multiplication", max_error_direct * 100.0f, avg_error_direct * 100.0f);
    printf("%-20s %-15.6f %-15.6f\n", "Sqrt reuse", max_error_sqrt_reuse * 100.0f, avg_error_sqrt_reuse * 100.0f);
    printf("%-20s %-15.6f %-15.6f\n", "Fast inverse sqrt", max_error_fast_inv_sqrt * 100.0f,
           avg_error_fast_inv_sqrt * 100.0f);

    // Free memory
    free(test_values);

    return 0;
}