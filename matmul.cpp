#include "egor_include.cpp"

struct matrix {
  union { unsigned int M; unsigned int Row; };
  union { unsigned int N; unsigned int Col; };
  float* Buffer;
};

struct arena {
  char* Start;
  unsigned int Size;
  unsigned int Used;
};

arena Arena;
double Freq;
double Iterations;

arena arenaAllocatorCreate(unsigned int SizeInBytes) {
  arena Arena = {nullptr, SizeInBytes, 0};
  Arena.Start = (char*)malloc(SizeInBytes);
  ASSERT(Arena.Start);
  return Arena;
}

void* arenaAllocatorAlloc(arena* Arena, unsigned int Size) {
  ASSERT(Arena->Size - Arena->Used >= Size);
  void* ReturnPtr = Arena->Start + Arena->Used;
  Arena->Used += Size;

  return ReturnPtr;
}


void zeroMatrix(matrix* M) {
  uint32_t Size = M->Row*M->Col;
  for (int i = 0; i < Size; ++i) { M->Buffer[i] = 0.0f; }
}

void initializeIdentity(matrix* M) {
  for (int i = 0; i < M->Row; ++i) {
    for (int j = 0; j < M->Col; ++j) {
      M->Buffer[i*M->Row + j] = (i == j) ? 1 : 0;
    }
  }
}

void initializeMatrix(matrix* M) {
  unsigned int Size = M->M*M->N;
  for (int i = 0; i < Size; ++i) M->Buffer[i] = 1;
}

matrix createMatrix(unsigned int M, unsigned int N) {
  matrix Result = {M, N, nullptr};
  Result.Buffer = (float*)arenaAllocatorAlloc(&Arena, Result.M*Result.N*sizeof(float));
  ASSERT(Result.Buffer);
  return Result; 
}

//[[gnu::noinline]]
void matmul_switched(matrix* Left, matrix* Right, matrix* Result) {
  TIME_FUNC();
  float* L = Left->Buffer;
  float* R = Right->Buffer;
  float* F = Result->Buffer;

  unsigned int Col = Left->M;
  unsigned int Common = Left->N;
  unsigned int Row = Right->N;

  for (int i = 0; i < Col; ++i) {
    for (int j = 0; j < Common; ++j) {
      for (int k = 0; k < Row; ++k) {
        F[i*Row + k] += L[i*Common + j] * R[j*Row + k]; 
      }
    }
  }
}

//[[gnu::noinline]]
void matmul_switched2(matrix* Left, matrix* Right, matrix* Result) {
  TIME_FUNC();

  float* L = Left->Buffer;
  float* R = Right->Buffer;
  float* F = Result->Buffer;

  unsigned int Col = Left->M;
  unsigned int Common = Left->N;
  unsigned int Row = Right->N;

  for (int i = 0; i < Col; ++i) {
    float* F_row = F + i*Row;
    float* L_row = L + i*Common;
    for (int j = 0; j < Common; ++j) {
      float L_val = L_row[j];
      float* R_row = R + j*Row;
#pragma unroll 4 
      for (int k = 0; k < Row; ++k) {
        F_row[k] += L_val * R_row[k]; 
      }
    }
  }
}


//[[gnu::noinline]]
void matmul_switched3(matrix* Left, matrix* Right, matrix* Result) {
  TIME_FUNC();

  float* L = Left->Buffer;
  float* R = Right->Buffer;
  float* F = Result->Buffer;

  unsigned int Col = Left->M;
  unsigned int Common = Left->N;
  unsigned int Row = Right->N;

  for (int i = 0; i < Col; ++i) {
    float* F_row = F + i*Row;
    float* L_row = L + i*Common;
    for (int j = 0; j < Common; ++j) {
      float L_val = L_row[j];
      float* R_row = R + j*Row;
#pragma unroll 4 
      for (int k = 0; k < Row; ++k) {
        float Res = L_val * R_row[k]; 
        F_row[k] += Res;
      }
    }
  }
}


void matmul_switched_micro(matrix* Left, matrix* Right, matrix* Result) {
  TIME_FUNC();

  float* L = Left->Buffer;
  float* R = Right->Buffer;
  float* F = Result->Buffer;

  unsigned int Col = Left->M;
  unsigned int Common = Left->N;
  unsigned int Row = Right->N;

#define MICRO_SIZE 512

  for (int k = 0; k < Row; k += MICRO_SIZE) {

    for (int i = 0; i < Col; ++i) {
      float Accum[MICRO_SIZE] = {};
      float* F_row = F + i*Row;
      float* L_row = L + i*Common;

      for (int j = 0; j < Common; ++j) {
        float L_val = L_row[j];
        float* R_row = R + j*Row;

        for (int micro_k = 0; micro_k < MICRO_SIZE; ++micro_k) {
          Accum[micro_k] += L_val * R_row[k + micro_k]; 
        }
      }

      for (int micro_k = 0; micro_k < MICRO_SIZE; ++micro_k) {
        F_row[k + micro_k] = Accum[micro_k];
      }
    }
  }
}


void matmul_switched_micro_256(matrix* Left, matrix* Right, matrix* Result) {
  TIME_FUNC();

  float* L = Left->Buffer;
  float* R = Right->Buffer;
  float* F = Result->Buffer;

  unsigned int Col = Left->M;
  unsigned int Common = Left->N;
  unsigned int Row = Right->N;

#define MICRO_DIM0 256
#define MICRO_DIM1 8 

  for (int k = 0; k < Row; k += MICRO_DIM0) {

    for (int i = 0; i < Col; ++i) {
      float Accum[MICRO_DIM0] = {};
      float* F_row = F + i*Row;
      float* L_row = L + i*Common;

      for (int j = 0; j < Common; j += MICRO_DIM1) {

        float L_val[MICRO_DIM1] = {};
        for (int micro_j = 0; micro_j < MICRO_DIM1; ++micro_j) {
          L_val[micro_j] = L_row[j + micro_j];
        }

        for (int micro_j = 0; micro_j < MICRO_DIM1; ++micro_j) {
          float* R_row = R + (j + micro_j)*Row;
          for (int micro_k = 0; micro_k < MICRO_DIM0; ++micro_k) {
            Accum[micro_k] += L_val[micro_j] * R_row[k + micro_k]; 
          }
        }
      }

      for (int micro_j = 0; micro_j < MICRO_DIM1; ++micro_j) {
        for (int micro_k = 0; micro_k < MICRO_DIM0; ++micro_k) {
          F_row[k + micro_k] = Accum[micro_k];
        }
      }

    }
  }
}


void matmul_vanila(matrix* Left, matrix* Right, matrix* Result) {
  TIME_FUNC();
  float* L = Left->Buffer;
  float* R = Right->Buffer;
  float* F = Result->Buffer;

  unsigned int Col = Left->M;
  unsigned int Common = Left->N;
  unsigned int Row = Right->N;

  for (int i = 0; i < Col; ++i) {
    for (int k = 0; k < Row; ++k) {
      for (int j = 0; j < Common; ++j) {
        F[i*Row + k] += L[i*Common + j] * R[j*Row + k]; 
      }
    }
  }
}


void isMatrixIdentity(matrix* M) {
  for (int i = 0; i < M->Row; ++i) {
    for (int j = 0; j < M->Col; ++j) {
      float Val = M->Buffer[i*M->Row + j];
      if ( i == j ) { ASSERT(Val); }
      else ASSERT(Val == 0.0f);
    }
  }
}

void printMatrix(matrix* M) {
  for (int i = 0; i < M->M; ++i) {
    printf("\n");
    for (int j = 0; j < M->N; ++j) {
      printf("%f ", M->Buffer[i * M->N + j]);
    }
  }
  printf("\n");
}


#define createMatmulGroup(Name, M, K, N) \
  matrix Name##_A = createMatrix(M, K);\
  matrix Name##_B = createMatrix(K, N);\
  matrix Name##_C = createMatrix(M, N);\
  initializeIdentity(&Name##_A);\
  initializeIdentity(&Name##_B);\
  //printMatrix(&Name##_B); \
  //printMatrix(&Name##_A); \

#define validateMatmulGroup(Name) \
  matmul_vanila(&Name##_A, &Name##_B, &Name##_C);\
  isMatrixIdentity(&Name##_C); \
  zeroMatrix(&Name##_C); \
  matmul_switched(&Name##_A, &Name##_B, &Name##_C);\
  isMatrixIdentity(&Name##_C); \
  zeroMatrix(&Name##_C); \
  matmul_switched2(&Name##_A, &Name##_B, &Name##_C);\
  isMatrixIdentity(&Name##_C); \
  zeroMatrix(&Name##_C); \
  matmul_switched3(&Name##_A, &Name##_B, &Name##_C);\
  isMatrixIdentity(&Name##_C); \
  zeroMatrix(&Name##_C); \
  matmul_switched_micro(&Name##_A, &Name##_B, &Name##_C);\
  isMatrixIdentity(&Name##_C); \
  zeroMatrix(&Name##_C); \
  matmul_switched_micro_256(&Name##_A, &Name##_B, &Name##_C);\
  isMatrixIdentity(&Name##_C); \
  zeroMatrix(&Name##_C); \
  PRINT_NL("==== End of validate ========");

#define testMatmulGroup(Name) \
  matmul_switched(&Name##_A, &Name##_B, &Name##_C);\
  matmul_switched(&Name##_A, &Name##_B, &Name##_C);\
  matmul_switched2(&Name##_A, &Name##_B, &Name##_C);\
  matmul_switched3(&Name##_A, &Name##_B, &Name##_C);\
  matmul_switched_micro(&Name##_A, &Name##_B, &Name##_C);\
  matmul_switched_micro_256(&Name##_A, &Name##_B, &Name##_C);\
  PRINT_NL("==== End of test ============");

int main(int argc, char** argv) {

  Freq = calibrateTimer();
  Iterations = (uint64_t)MICRO_SIZE*MICRO_SIZE*MICRO_SIZE*4*8*4; 

  TIME(main);
  Arena = arenaAllocatorCreate(2E9);
  createMatmulGroup(Test, MICRO_SIZE*1, MICRO_SIZE*1, MICRO_SIZE*1);
  validateMatmulGroup(Test);

  createMatmulGroup(Perf, MICRO_SIZE*8, MICRO_SIZE*8, MICRO_SIZE*4);
  testMatmulGroup(Perf);
}
