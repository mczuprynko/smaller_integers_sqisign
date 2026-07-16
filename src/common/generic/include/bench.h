// SPDX-License-Identifier: Apache-2.0
#ifndef BENCH_H__
#define BENCH_H__

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <execinfo.h>
#if defined(__APPLE__)
#include "bench_macos.h"
#endif

#if defined(TARGET_ARM) || defined(TARGET_S390X) || defined(NO_CYCLE_COUNTER)
#define BENCH_UNIT0 "nanoseconds"
#define BENCH_UNIT3 "microseconds"
#define BENCH_UNIT6 "milliseconds"
#define BENCH_UNIT9 "seconds"
#else
#define BENCH_UNIT0 "cycles"
#define BENCH_UNIT3 "kilocycles"
#define BENCH_UNIT6 "megacycles"
#define BENCH_UNIT9 "gigacycles"
#endif

#define MAX_FRAMES 50
#define MAX_SAMPLES 600000
#define MAX_BIT_SAMPLES 600000
// #define RECORD_CALLSTACK

extern uint64_t lll_count;
extern unsigned long long MAX_bit;
extern unsigned long long count_calls;
extern char* frame_array[MAX_SAMPLES][MAX_FRAMES];
extern unsigned bit_samples_array[MAX_BIT_SAMPLES];

enum MeasureFunc {
    COMMIT,
    COMMIT_RndIdealGivNorm,
    COMMIT_RndEquivPrimeIdeal,
    COMMIT_IdealToIso,
    COMMIT_ITI_SuitableIdeals,
    COMMIT_ITI_FixedDegIso,
    COMMIT_ITI_Iso22Chain,
    HASH_CHALL,
    CHALL_TO_QUAT,
    RAND_AUX_IDEAL,
    ODD_ISO,
    ODD_RndIdealGivNorm,
    ODD_IdealIntersection,
    ODD_IdealToIso,
    ODD_ITI_SuitableIdeals,
    ODD_ITI_FixedDegIso,
    ODD_ITI_Iso22Chain,
    ODD_SplitAuxIso,
    EVEN_ISO,
    END_SIG,
    NUM_MEASURMENTS
};

extern const char *NamesMeasures[NUM_MEASURMENTS];

extern const int groupings[NUM_MEASURMENTS];

extern uint64_t cyc_data[NUM_MEASURMENTS];

extern uint64_t choice;
extern int FileDesc;

/*#define GET_BACKTRACE_SYMBOLS(array, curr_counts, bit)                     \
    if (curr_counts < MAX_SAMPLES)                          \
    {                                                        \
        if (0) {                                        \
            void *callstack[MAX_FRAMES];      \
            int frames = backtrace(callstack, MAX_FRAMES);               \
            char **symbols = backtrace_symbols(callstack, frames); \
            if (symbols == NULL) {                                  \
                perror("backtrace_symbols");                       \
                exit(EXIT_FAILURE);                                 \
            }                                                       \
            for (int i = 0; i < frames; i++) {                     \
                array[i] = symbols[i];                              \
            }                                                       \
            array[frames] = 0;                                      \
            free(symbols);                                          \
        }                                                           \
    }                                                               \
    if (curr_counts < MAX_BIT_SAMPLES) {                        \
        bit_samples_array[count_calls] = bit;                 \
    } */

#ifdef RECORD_CALLSTACK
    #define GET_BACKTRACE_SYMBOLS(array, curr_counts, bit)                     \
        if (curr_counts < MAX_SAMPLES)                          \
        {                                                                      \
            void *callstack[MAX_FRAMES];      \
            int frames = backtrace(callstack, MAX_FRAMES);               \
            backtrace_symbols_fd(callstack, frames, FileDesc); \
            char var_str[] = "\n";\
            int ret =  write(FileDesc, var_str, strlen(var_str)); \
            (void) ret;             \
        }                                                               \
        if (curr_counts < MAX_BIT_SAMPLES) {                        \
            bit_samples_array[count_calls] = bit;                 \
        }
#else 
    #define GET_BACKTRACE_SYMBOLS(array, curr_counts, bit) 
#endif 

#define PRINT_BACKTRACE_SYMBOLS(bit)                                                \
        void *callstack[MAX_FRAMES];                                  \
        int frames = backtrace(callstack, MAX_FRAMES);                \
        backtrace_symbols_fd(callstack, frames, STDOUT_FILENO);       \
        int ret = write(STDOUT_FILENO, "\n", 1);                        \
        (void) ret;                                                     \
        printf("Bit size: %lu\n", bit);                                \

static inline void
cpucycles_init(void) {
#if defined(__APPLE__) && defined(TARGET_ARM64)
    macos_init_rdtsc();
#endif
}

static inline uint64_t
cpucycles(void)
{
#if defined(TARGET_AMD64) || defined(TARGET_X86)
    uint32_t hi, lo;

    asm volatile("rdtsc" : "=a"(lo), "=d"(hi));
    return ((uint64_t)lo) | ((uint64_t)hi << 32);
#elif defined(TARGET_S390X)
    uint64_t tod;
    asm volatile("stckf %0\n" : "=Q"(tod) : : "cc");
    return (tod * 1000 / 4096);
#elif defined(TARGET_ARM64) && !defined(NO_CYCLE_COUNTER)
#if defined(__APPLE__)
    return macos_rdtsc();
#else
    uint64_t cycles;
    asm volatile("mrs %0, PMCCNTR_EL0" : "=r"(cycles));
    return cycles;
#endif // __APPLE__
#else
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    return (uint64_t)time.tv_sec * 1000000000 + time.tv_nsec;
#endif
}

#define BEG_MES() { uint64_t cyc_count = cpucycles();

#define END_MES(store_var_id) cyc_data[store_var_id] += cpucycles() - cyc_count;}

#define END_MES_CHOICE(store_var_id0, store_var_id1) \
    if (choice == 0) cyc_data[store_var_id0] += cpucycles() - cyc_count;\
    else cyc_data[store_var_id1] += cpucycles() - cyc_count;\
}

static inline int
CMPFUNC(const void *a, const void *b)
{
    uint64_t aa = *(uint64_t *)a, bb = *(uint64_t *)b;

    if (aa > bb)
        return +1;
    if (aa < bb)
        return -1;
    return 0;
}

static inline uint32_t
ISQRT(uint64_t x)
{
    uint32_t r = 0;
    for (ssize_t i = 31; i >= 0; --i) {
        uint32_t s = r + (1 << i);
        if ((uint64_t)s * s <= x)
            r = s;
    }
    return r;
}

static inline double
_TRUNC(uint64_t x)
{
    return x / 1000 / 1000.;
}
#define _FMT ".3lf"
#define _UNIT BENCH_UNIT6

#define BENCH_CODE_1(RUNS)                                                                         \
    {                                                                                              \
        const size_t count = (RUNS);                                                               \
        if (!count)                                                                                \
            abort();                                                                               \
        uint64_t cycles, cycles1, cycles2;                                                         \
        uint64_t cycles_list[count];                                                               \
        cycles = 0;                                                                                \
        for (size_t i = 0; i < count; ++i) {                                                       \
            cycles1 = cpucycles();

#define BENCH_CODE_2(name)                                                                         \
    cycles2 = cpucycles();                                                                         \
    cycles_list[i] = cycles2 - cycles1;                                                            \
    cycles += cycles2 - cycles1;                                                                   \
    }                                                                                              \
    qsort(cycles_list, count, sizeof(uint64_t), CMPFUNC);                                          \
    uint64_t variance = 0;                                                                         \
    for (size_t i = 0; i < count; ++i) {                                                           \
        int64_t off = cycles_list[i] - cycles / count;                                             \
        variance += off * off;                                                                     \
    }                                                                                              \
    variance /= count;                                                                             \
    printf("  %-10s", name);                                                                       \
    printf(" | average %9" _FMT " | stddev %9" _FMT,                                               \
           _TRUNC(cycles / count),                                                                 \
           _TRUNC(ISQRT(variance)));                                                               \
    printf(" | median %9" _FMT " | min %9" _FMT " | max %9" _FMT,                                  \
           _TRUNC(cycles_list[count / 2]),                                                         \
           _TRUNC(cycles_list[0]),                                                                 \
           _TRUNC(cycles_list[count - 1]));                                                        \
    printf("  (%s)\n", _UNIT);                                                                     \
    }

#define PRINT_STATS(RUNS, vals)                                                                             \
    {                                                                                                       \
    int width = 0;                                                                                          \
    for (int j = 0; j < NUM_MEASURMENTS; ++j){                                                              \
        if (width < groupings[j]) width = groupings[j];                                                     \
    }                                                                                                       \
    for (int j = 0; j < NUM_MEASURMENTS; ++j){                                                              \
        uint64_t c_vals[RUNS];                                                                              \
        memcpy(c_vals, &vals[j], sizeof(c_vals));                                                           \
        qsort(c_vals, RUNS, sizeof(uint64_t), CMPFUNC);                                                     \
        uint64_t variance = 0;                                                                              \
        uint64_t cycles = 0;                                                                                \
        for (size_t i = 0; i < RUNS; ++i) {                                                                 \
            cycles += c_vals[i];                                                                            \
        }                                                                                                   \
        for (size_t i = 0; i < RUNS; ++i) {                                                                 \
            int64_t off = c_vals[i] - cycles / RUNS;                                                        \
            variance += off * off;                                                                          \
        }                                                                                                   \
        variance /= RUNS;                                                                                   \
        if (groupings[j] == 0) {                                                                            \
            for (int i = 0; i < 150; ++i) printf("=");                                                      \
            printf("\n");                                                                                   \
        }                                                                                                   \
        printf("  ");                                                                                       \
        for (int k =0; k < 4; ++k) printf("%*s", groupings[j], "");                                         \
        printf("%-20s", NamesMeasures[j]);                                                                  \
        for (int k =0; k < 4; ++k) printf("%*s", width - groupings[j], "");                                 \
        printf(" | average %9" _FMT " | stddev %9" _FMT,                                                    \
            _TRUNC(cycles / RUNS),                                                                          \
            _TRUNC(ISQRT(variance)));                                                                       \
        printf(" | median %9" _FMT " | min %9" _FMT " | max %9" _FMT,                                       \
            _TRUNC(c_vals[RUNS / 2]),                                                                       \
            _TRUNC(c_vals[0]),                                                                              \
            _TRUNC(c_vals[RUNS - 1]));                                                                      \
        printf("  (%s)\n", _UNIT);                                                                          \
        if ((groupings[j] != 0) && ((j == (NUM_MEASURMENTS-1)) || (groupings[j] > groupings[j+1]))){        \
            for (int curr_level = groupings[j]; curr_level > groupings[j+1]; curr_level--){                 \
                uint64_t tot_variance = 0;                                                                  \
                uint64_t tot_cycles = 0;                                                                    \
                uint64_t sub_tot_cycles[RUNS];                                                              \
                memset(sub_tot_cycles, 0, sizeof(sub_tot_cycles));                                          \
                for (size_t i = 0; i < RUNS; ++i) {                                                         \
                    int jj = j;                                                                             \
                    while ((jj >= 0) && (curr_level <= groupings[jj])){                                     \
                        if (curr_level != groupings[jj]) {jj -= 1; continue;}                               \
                        sub_tot_cycles[i] += vals[jj][i];                                                   \
                        jj -= 1;                                                                            \
                    }                                                                                       \
                    tot_cycles += sub_tot_cycles[i];                                                        \
                }                                                                                           \
                for (size_t i = 0; i < RUNS; ++i) {                                                         \
                    int64_t off = sub_tot_cycles[i] - tot_cycles / RUNS;                                    \
                    tot_variance += off * off;                                                              \
                }                                                                                           \
                qsort(sub_tot_cycles, RUNS, sizeof(uint64_t), CMPFUNC);                                     \
                tot_variance /= RUNS;                                                                       \
                printf("  ");                                                                               \
                for (int k =0; k < 4; ++k) printf("%*s", curr_level, "");                                   \
                printf("%-20s", "Sub_total");                                                               \
                for (int k =0; k < 4; ++k) printf("%*s", width - curr_level, "");                           \
                printf(" | average %9" _FMT " | stddev %9" _FMT,                                            \
                    _TRUNC(tot_cycles / RUNS),                                                              \
                    _TRUNC(ISQRT(tot_variance)));                                                           \
                printf(" | median %9" _FMT " | min %9" _FMT " | max %9" _FMT,                               \
                    _TRUNC(sub_tot_cycles[RUNS / 2]),                                                       \
                    _TRUNC(sub_tot_cycles[0]),                                                              \
                    _TRUNC(sub_tot_cycles[RUNS - 1]));                                                      \
                printf("  (%s)\n", _UNIT);                                                                  \
            }                                                                                               \
        }                                                                                                   \
    }                                                                                                       \
    }
#endif
