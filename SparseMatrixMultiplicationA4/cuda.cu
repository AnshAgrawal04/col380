#include "header.h"
#include <cuda_runtime.h>
#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstring>                       

__global__
void bcsr_mul_kernel(int nb_rows, int k,
                     const int* __restrict__ A_rowptr,
                     const int* __restrict__ A_colind,
                     const unsigned long long* __restrict__ A_vals,
                     const int* __restrict__ B_rowptr,
                     const int* __restrict__ B_colind,
                     const unsigned long long* __restrict__ B_vals,
                     const int* __restrict__ C_rowptr,
                     const int* __restrict__ C_colind,
                     unsigned long long* C_vals)
{
    int i   = blockIdx.x;            
    int tid = threadIdx.x; 
    int elems = k * k;

    extern __shared__ unsigned long long shmem[];
    unsigned long long* Ash = shmem;
    unsigned long long* Bsh = shmem + elems;

    int bi = tid / k;
    int bj = tid % k;

    for (int a = A_rowptr[i]; a < A_rowptr[i+1]; ++a) {
        int p = A_colind[a];
        const unsigned long long* Ablk = A_vals + size_t(a) * elems;
        Ash[tid] = Ablk[tid];
        __syncthreads();

        for (int b = B_rowptr[p]; b < B_rowptr[p+1]; ++b) {
            int j = B_colind[b];
            const unsigned long long* Bblk = B_vals + size_t(b) * elems;
            Bsh[tid] = Bblk[tid];
            __syncthreads();

            unsigned long long sum = 0;
            for (int kk = 0; kk < k; ++kk) {
                sum += Ash[bi*k + kk] * Bsh[kk*k + bj];
            }

            int low = C_rowptr[i];
            int high = C_rowptr[i+1] - 1;
            while (low <= high) {
                int mid = (low + high) >> 1;
                int col = C_colind[mid];
                if (col == j) { low = mid; break; }
                if (col < j) low = mid + 1;
                else         high = mid - 1;
            }
            int pos = low;

            // Atomic accumulate
            atomicAdd(&C_vals[size_t(pos)*elems + bi*k + bj], sum);
            __syncthreads();
        }
    }
}

void CSR_mul(struct BCSR& A, struct BCSR& B, struct BCSR& C, int k) {

    int nb_rows = (A.rows+(k-1)) / k;
    vector<std::unordered_map<int,bool>> rowBlocks(nb_rows);
    for (int i = 0; i < nb_rows; ++i) {
        for (int a = A.rowptr[i]; a < A.rowptr[i+1]; ++a) {
            int p = A.colind[a];
            for (int b = B.rowptr[p]; b < B.rowptr[p+1]; ++b) {
                rowBlocks[i][B.colind[b]] = true;
            }
        }
    }

    C.rows = A.rows;
    C.cols = B.cols;    
    C.rowptr  = new int[nb_rows + 1];
    C.rowptr[0]  = 0;


    for (int i = 0; i < nb_rows; ++i)
        C.rowptr[i+1] = C.rowptr[i] + int(rowBlocks[i].size());


    int nnzC = C.rowptr[nb_rows];

    C.nz = nnzC;
    // cout<<"finally c nnz "<<nnzC<<endl;
    C.colind = new int[nnzC];
    int idx = 0;

    for (int i = 0; i < nb_rows; ++i) {
        std::vector<int> cols;
        cols.reserve(rowBlocks[i].size());
        for (auto &kv : rowBlocks[i]) cols.push_back(kv.first);
        std::sort(cols.begin(), cols.end());
        for (int col : cols) C.colind[idx++] = col;
    }

    // Allocate and zero C.values
    size_t elems = size_t(k) * k;
    C.values = new unsigned long long[nnzC * elems];
    std::memset(C.values, 0, sizeof(unsigned long long) * nnzC * elems);

    // Device buffers
    int *d_A_rowptr, *d_A_colind, *d_B_rowptr, *d_B_colind, *d_C_rowptr, *d_C_colind;
    unsigned long long *d_A_vals, *d_B_vals, *d_C_vals;

    // Helper: alloc+copy
    auto alloc_and_copy = [&](const void* src, size_t bytes, void** dst) {
        cudaMalloc(dst, bytes);
        cudaMemcpy(*dst, src, bytes, cudaMemcpyHostToDevice);
    };
    int padded_b = (B.rows + k -1)/k ;
    alloc_and_copy(A.rowptr, (nb_rows+1)*sizeof(int), (void**)&d_A_rowptr);
    alloc_and_copy(A.colind, A.rowptr[nb_rows]*sizeof(int), (void**)&d_A_colind);
    alloc_and_copy(A.values, elems*A.rowptr[nb_rows]*sizeof(unsigned long long), (void**)&d_A_vals);
    alloc_and_copy(B.rowptr, (padded_b+1)*sizeof(int),   (void**)&d_B_rowptr);
    alloc_and_copy(B.colind, B.rowptr[padded_b]*sizeof(int), (void**)&d_B_colind);
    alloc_and_copy(B.values, elems*B.rowptr[padded_b]*sizeof(unsigned long long), (void**)&d_B_vals);
    alloc_and_copy(C.rowptr, (nb_rows+1)*sizeof(int),  (void**)&d_C_rowptr);
    alloc_and_copy(C.colind, nnzC*sizeof(int),  (void**)&d_C_colind);

    cudaMalloc((void**)&d_C_vals, elems*nnzC*sizeof(unsigned long long));
    cudaMemset(d_C_vals, 0, elems*nnzC*sizeof(unsigned long long));

    dim3 grid(nb_rows);
    dim3 block(k * k);
    size_t sharedMem = 2ULL * k * k * sizeof(unsigned long long);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess){
        fprintf(stderr, "Kernel launch failed: %s\n", cudaGetErrorString(err));
        cout<<"This one error"<<endl;
    }
        

    bcsr_mul_kernel<<<grid, block, sharedMem>>>(
        nb_rows, k,
        d_A_rowptr, d_A_colind, d_A_vals,
        d_B_rowptr, d_B_colind, d_B_vals,
        d_C_rowptr, d_C_colind, d_C_vals
    );
    err = cudaGetLastError();
    if (err != cudaSuccess)
        fprintf(stderr, "Kernel launch failed: %s\n", cudaGetErrorString(err));

    cudaDeviceSynchronize();

    cudaMemcpy(C.values, d_C_vals, elems*nnzC*sizeof(unsigned long long), cudaMemcpyDeviceToHost);

    cudaFree(d_A_rowptr); cudaFree(d_A_colind); cudaFree(d_A_vals);
    cudaFree(d_B_rowptr); cudaFree(d_B_colind); cudaFree(d_B_vals);
    cudaFree(d_C_rowptr); cudaFree(d_C_colind); cudaFree(d_C_vals);
}