#include <iostream>
#include <random>
#include <chrono>
#include <vector>

struct crsMatrix {
    size_t N;   // matrix size (N x N)
    size_t NZ;  // number of non-zero elements
    std::vector<double> value;     // value array (size is NZ)
    std::vector<size_t> column;    // column number array (size is NZ)
    std::vector<size_t> rowIndex;  // row index array (size is N + 1)

    crsMatrix(size_t N, size_t NZ) : N(N), NZ(NZ) {
        value.resize(NZ, 0.0);
        column.resize(NZ, 0);
        rowIndex.resize(N + 1, 0);
    }

    friend bool operator==(const crsMatrix& m1, const crsMatrix& m2) {
        return m1.N == m2.N && m1.NZ == m2.NZ &&
            m1.value == m2.value && m1.column == m2.column && m1.rowIndex == m2.rowIndex;
    }
};

typedef void (*TransposeFuncType)(const crsMatrix&, crsMatrix&);

// matrix generator
void generateRegularCRS(size_t N, size_t cntInRow, crsMatrix& mtx);
// function to check correctness
void testCorrectness(TransposeFuncType transposeFunc);

// naive transpose implementation
__declspec(noinline)
void transposeNaive(const crsMatrix& m, crsMatrix& mT) {
    // create additional vectors
    
    for (size_t i = 0; i < m.N; i++) {
        size_t k1 = m.rowIndex[i], k2 = m.rowIndex[i + 1];
        for (size_t k = k1; k < k2; k++) {
            // push m.value[k] and i to additional vectors
        }
    }
    
    // copy arrays from additional vectors to mT
}

// optimized transpose implementation
__declspec(noinline)
void transposeOpt(const crsMatrix& m, crsMatrix& mT) {
    
    // calculate mT.rowIndex:
    // 1. calculate number of elements in every column
    // 2. calculate mT.rowIndex (use partial sums)
    // ...
    
    // next 2 lines: calculation of start indices (shift)
    for (size_t j = mT.N; j >= 1; j--)
        mT.rowIndex[j] = mT.rowIndex[j - 1];
    
    // generate mT
    for (size_t i = 0; i < m.N; i++) {
        size_t k1 = m.rowIndex[i], k2 = m.rowIndex[i + 1];
        for (size_t k = k1; k < k2; k++) {
            size_t j = m.column[k];
            // fill mT.column and mT.values
            ++mT.rowIndex[j + 1];
        }
    }
}

int main() {
    
    // choose function to test
    TransposeFuncType transposeFunc = transposeNaive;
    //TransposeFuncType transposeFunc = transposeOpt;

    // check correctness of your implementation
    testCorrectness(transposeFunc);

    // do performance experiments
    const size_t N = 100000;
    const size_t cntInRow = 500;
    const size_t NZ = N * cntInRow;
    std::cout << "Matrix size is ~" <<
        (NZ * sizeof(double) + (NZ + N + 1) * sizeof(size_t)) / 1024.0 / 1024.0 / 1024.0 <<
        " GB" << std::endl;

    crsMatrix m(N, NZ), mT(N, NZ);

    const auto start{ std::chrono::steady_clock::now() };
    
    //for (int i = 0; i < 1000; i++)
    transposeFunc(m, mT);

    const auto end{ std::chrono::steady_clock::now() };
    const std::chrono::duration<double> elapsed_seconds{ end - start };

    std::cout << "Time of transposing is " << elapsed_seconds.count() << " s" << std::endl;
    
    return 0;
}

void testCorrectness(TransposeFuncType transposeFunc) {
    crsMatrix m(6, 9);
    m.value = {1, 2, 3, 4, 8, 5, 7, 1, 6};
    m.column = {0, 4, 2, 3, 3, 5, 1, 2, 5};
    m.rowIndex = {0, 2, 4, 4, 6, 6, 9};
    
    crsMatrix mTExp(6, 9);
    mTExp.value = {1, 7, 3, 1, 4, 8, 2, 5, 6};
    mTExp.column = {0, 5, 1, 5, 1, 3, 0, 3, 5};
    mTExp.rowIndex = {0, 1, 2, 4, 6, 7, 9};
    
    crsMatrix mT(6, 9);
    transposeFunc(m, mT);
    
    if (mT == mTExp) std::cout << "OK" << std::endl;
    else std::cout << "ERROR: incorrect implementation" << std::endl;
}

void generateRegularCRS(size_t N, size_t cntInRow, crsMatrix& mtx) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> idis(0, N - 1);
    std::uniform_real_distribution<> rdis(0, 1);
    size_t NZ = cntInRow * N;
    mtx.N = N; mtx.NZ = NZ;
    mtx.value.resize(NZ);
    mtx.column.resize(NZ);
    mtx.rowIndex.resize(N + 1);
    for (size_t i = 0; i < N; i++) {
        // form column numbers in the line i
        for (size_t j = 0; j < cntInRow; j++) {
            int f = 0;
            do {
                mtx.column[i * cntInRow + j] = idis(gen);
                f = 0;
                for (size_t k = 0; k < j; k++)
                    if (mtx.column[i * cntInRow + j] == mtx.column[i * cntInRow + k])
                        f = 1;
            } while (f == 1);
        }
        // sort column numbers in the line i
        for (size_t j = 0; j < cntInRow - 1; j++)
            for (size_t k = 0; k < cntInRow - 1; k++)
                if (mtx.column[i * cntInRow + k] > mtx.column[i * cntInRow + k + 1]) {
                    size_t tmp = mtx.column[i * cntInRow + k];
                    mtx.column[i * cntInRow + k] =
                    mtx.column[i * cntInRow + k + 1];
                    mtx.column[i * cntInRow + k + 1] = tmp;
                }
    }
    // fill value array
    for (int i = 0; i < cntInRow * N; i++)
        mtx.value[i] = rdis(gen);
    // fill row index array
    mtx.rowIndex[0] = 0;
    for (int i = 1; i <= N; i++)
        mtx.rowIndex[i] = mtx.rowIndex[i + 1] + cntInRow;
}
