#include <bits/stdc++.h>    
#include <iostream>
#include <mpi.h>
#include <fstream>


using namespace std;

struct BlockCSRMatrix {
    int rows, cols;               // Full matrix size (in elements, not blocks)
    int block_size;               // Size of each block (assuming square k x k blocks)
    int num_block_rows;            // Number of block rows
    int num_block_cols;            // Number of block columns
    int num_nonzeros;           // Number of non-zero blocks
    int orig_rows, orig_cols; // Original matrix size (in elements, not blocks)

    vector<uint64_t*> block_values; // Each entry points to a block (size block_size*block_size)
    vector<int> block_col_indices;  // Column index (which block column)
    vector<int> block_row_ptr;      // Starting offset of blocks for each block-row

    BlockCSRMatrix(int r, int c, int bsize) 
        : rows(r), cols(c), block_size(bsize) {
        num_block_rows = (r + bsize - 1) / bsize;
        num_block_cols = (c + bsize - 1) / bsize;
        block_row_ptr.push_back(0); // Like CSR, starts with 0
    }
};

struct BCSR {
    int rows, cols;
    int nz;
    int* rowptr;
    int* colind;
    unsigned long long* values;
};

void CSR_mul(struct BCSR &A, struct BCSR &B , struct BCSR &C, int k);
