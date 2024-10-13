#include "stdio.h"
#include <cstdint>
#include <unistd.h>
#include "time.h"
#include <cstdlib>

#define DISABLE_DIAGNOSTICS \
  _Pragma("clang diagnostic push")\
  _Pragma("clang diagnostic ignored \"-Wint-to-pointer-cast\"")

#define ENABLE_DIAGNOSTICS _Pragma("clang diagnostic pop")

#define ASSERT(Expression)                                                              \
  DISABLE_DIAGNOSTICS \
    if (!(Expression)) {                                                                \
        printf("ASSERT FAILED: %s - %s - line: %d\n", #Expression, __FILE__, __LINE__); \
        volatile int a = 0;                                                             \
        *(int*) a = 0;                                                                  \
    }\
  ENABLE_DIAGNOSTICS


#define PRINT_NL(str) printf("%s \n", str);

uint64_t _rdtsc() {
  uint64_t lo, hi = 0;
  asm volatile( "rdtsc" : "=a" (lo), "=d" (hi) );
  return (lo | (hi << 32));
}

uint64_t calibrateTimer() {

  timespec TS = {0, 100000000};
  uint64_t Accum = 0;
  for (int i = 0; i < 30; ++i) {
    timespec Rem = {};
    uint64_t St = _rdtsc();
    nanosleep(&TS, &Rem);
    Accum += _rdtsc() - St;
  }

  double Freq = Accum/3.0;
  printf("Cycles per second: %lf  \n", Freq);
  return Freq;
}

struct timer { 
  uint64_t Start;
  double Freq;
  const char* Name;
  double Iter;

  timer(double VFreq, double Iterations, const char* VName) : Freq(VFreq), Name(VName), Iter(Iterations) { Start = _rdtsc(); }
  ~timer() { 
    uint64_t Elapsed = _rdtsc() - Start;
    printf(" %s --> Elapsed cycles: %ld -- Elapsed Time: %lf \n", Name, Elapsed, Elapsed/Freq);
    printf(" %s --> Cycles per iteration: %lf \n", Name, Elapsed/(Iter));
  }
};

#define TIME(Name) timer Name##_Timer(Freq, Iterations, #Name);
#define TIME_FUNC() timer Timer(Freq, Iterations, __func__);
