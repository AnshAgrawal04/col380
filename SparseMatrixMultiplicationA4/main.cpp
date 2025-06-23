#include "header.h"
#include <omp.h>

bool check_args(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <folder_path>" << std::endl;
        return false;
    }
    return true;
}

void read_size_file(int rank, int& size, string folder_path , int &k , int &start, int &end, int &flag){
    ifstream infile(folder_path + "/size");
    if (!infile) {
        cerr << "Error opening file" << endl;
        return;
    }
    int n;
    infile >> n >> k;
    infile.close();

    if (n < size){
        if (rank < n){ 
            start = rank ; 
            end = rank+1;
            size = n;
            return;
        }
        else{flag = -1;return;}
    }
    int r = 0 ;
    start = 0;
    end = 0 ;
    while (r <= rank){
        int cur_td = n/size + (1?(r)<n%size:0);
        start = end;
        end += cur_td;
        r++;
    }
    cout<<"I am rank: "<<rank<<" and I am reading from "<< start<< " "<<end<<endl;
    cout << "Read n and k from file: " << n << " and " << k << endl;

}

void read_file(string folder_path , int i , int k , BCSR &input){
    ifstream infile(folder_path + "/matrix" + to_string(i+1));
    if (!infile) {
        cerr << "Error opening file" << endl;
        return;
    }
    int rows, cols, nz;
    infile >> rows >> cols >> nz;

    int padded_rows = ((rows + k - 1) / k) * k;
    // int padded_cols = ((cols + k - 1) / k) * k;
    int nb_rows = padded_rows / k;

    vector<int> block_rows(nz);
    vector<int> block_cols(nz);
    vector<unsigned long long*> block_values(nz);

    for (int j = 0; j < nz; ++j) {
        int row, col;
        infile >> row >> col;
        block_rows[j] = row / k;
        block_cols[j] = col / k;

        unsigned long long* block = new unsigned long long[k * k];
        // if (row + k > rows || col + k > cols) {
            std::fill_n(block, k * k, 0ULL); 
            int x = min(k, rows - row);
            int y = min(k, cols - col);
            for (int a = 0; a < x; ++a) {
                for (int b = 0; b < y; ++b) {
                    infile >> block[a * k + b];
                    // cout<<block[a * k + b];
                }
            }
            // cout<<endl;
        // }
        // else{
        //     for (int a = 0 ;a<k*k;a++){
        //         infile >> block[a];
        //     }
        // }

        block_values[j] = block;
    }
    infile.close();

    vector<int> row_counts(nb_rows,0);
    for (int r : block_rows) {
        ++row_counts[r];
    }
    int* rowptr = new int[nb_rows+1];
    rowptr[0] = 0;
    for (int r = 0; r < nb_rows; ++r) {
        rowptr[r + 1] = rowptr[r] + row_counts[r];
    }   
    int nnz_blocks = rowptr[nb_rows];
    cout<<"Read Matrix nnz blocks "<<nnz_blocks<<" and nz "<<nz<<endl;
    int* colind = new int[nnz_blocks];
    unsigned long long* values = new unsigned long long[nnz_blocks*k*k];

    vector<int> cursor(nb_rows);
    for (int r = 0; r < nb_rows; ++r) cursor[r] = rowptr[r];

    for (int j = 0; j < nnz_blocks; ++j) {
        int r = block_rows[j];
        int pos = cursor[r]++;
        colind[pos] = block_cols[j];
        memcpy(values+pos*(k*k),block_values[j],sizeof(unsigned long long)*(k*k));
        delete[] block_values[j];
    }

    input.rows =  rows;
    input.cols = cols;
    input.nz = nz;
    input.rowptr = rowptr;
    input.colind = colind;
    input.values = values;
}

void send_bcsr(int dest, BCSR &input, int k){
    int padded_rows = (input.rows + k - 1) / k;
    
    int packet[3];
    packet[0] = input.rows; packet[1] = input.cols; packet[2] = input.nz;
    // make the send non blpcking
    MPI_Send(packet , 3 , MPI_INT , dest , 0 , MPI_COMM_WORLD);
    MPI_Send(input.rowptr , padded_rows+1 , MPI_INT , dest , 0 , MPI_COMM_WORLD);
    MPI_Send(input.colind , input.nz , MPI_INT , dest , 0 , MPI_COMM_WORLD);
    MPI_Send(input.values , input.nz*k*k , MPI_UNSIGNED_LONG_LONG , dest , 0 , MPI_COMM_WORLD);
}

void recv_bcsr(int src, BCSR &input, int k){
    int packet[3];
    MPI_Recv(packet , 3 , MPI_INT , src , 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    input.rows = packet[0]; input.cols = packet[1]; input.nz = packet[2];
    int padded_rows = (input.rows + k - 1) / k;
    input.rowptr = new int[padded_rows+1];
    input.colind = new int[input.nz];
    input.values = new unsigned long long[input.nz*k*k];
    MPI_Recv(input.rowptr , padded_rows+1 , MPI_INT , src , 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(input.colind , input.nz , MPI_INT , src , 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(input.values , input.nz*k*k , MPI_UNSIGNED_LONG_LONG , src , 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // add mpi wait here
}


void write_final_ans(string folder_path, BCSR &final_result, int k){
    ofstream outfile("matrix");
    if (!outfile) {
        cerr << "Error opening file" << endl;
        return;
    }
    outfile << final_result.rows << " " << final_result.cols << " \n" << final_result.nz << endl;
    int nb_rows = (final_result.rows + k - 1) / k;
    cout<<"NB: "<<nb_rows<<endl;
    for (int i = 0; i < nb_rows; ++i) {
        int row0 = i * k;
        for (int ptr = final_result.rowptr[i]; ptr < final_result.rowptr[i + 1]; ++ptr) {
            int bj = final_result.colind[ptr];
            int col0 = bj * k;
            outfile << row0 << " " << col0 << endl;
            // cout<< "Printing "<<row0<<" "<<col0<<endl;
            // bool edge_row = (row0 + k > final_result.rows);
            // bool edge_col = (col0 + k > final_result.cols);

            // if (edge_row || edge_col) {
                int x = min(k, final_result.rows - row0);
                int y = min(k, final_result.cols - col0);
                for (int a = 0; a < x; ++a) {
                    for (int b = 0; b < y; ++b) {
                        outfile << final_result.values[ptr * k * k + a * k + b];
                        if (b== y-1) outfile << "\n";
                        else outfile << " ";
                    }
                }
            // } else {
            //     for (int a = 0; a < k; ++a) {
            //         for (int b = 0; b < k; ++b) {
            //             outfile << final_result.values[ptr * k * k + a * k + b];
            //             if (b == k-1) outfile << "\n";
            //             else outfile << " ";
            //         }
            //     }
            // }
        }
    }

    outfile.close();
}

int main(int argc, char* argv[]) {
    if (!check_args(argc, argv)) { return 1; }

    int rank, size;
    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int start,end; int flag = 0;
    int k ;
    string folder_path = argv[1]; 
    read_size_file(rank, size, folder_path , k , start, end, flag);

    if (flag == -1){
        MPI_Finalize();
        cout<<"Exiting since nothing can be done"<<endl;
        return 0;
    }

    vector<BCSR> inputs(end - start);
    int num_threads = omp_get_max_threads(); // Get the maximum number of threads available

    #pragma omp parallel for num_threads(num_threads)
    for (int i = start; i < end; i++) {
        read_file(folder_path, i, k, inputs[i - start]);
    }
    // Perform matrix multiplication
    int total_matrices = end - start;
    while (total_matrices > 1){
        int new_total_matrices = (total_matrices + 1) / 2;
        vector<BCSR> new_inputs(new_total_matrices);
        #pragma omp parallel for num_threads(num_threads)
        for (int i = 0; i < new_total_matrices; i++) {
            int a = i * 2;
            int b = a + 1;
            if (b >= total_matrices) {
                new_inputs[i] = inputs[a];
            } else {
                BCSR C;
                CSR_mul(inputs[a], inputs[b], C, k);
                new_inputs[i] = C;
                delete[] inputs[a].rowptr;
                delete[] inputs[a].colind;
                delete[] inputs[a].values;
                delete[] inputs[b].rowptr;
                delete[] inputs[b].colind;
                delete[] inputs[b].values;
            }
        }
        if (total_matrices % 2 == 1) {
            new_inputs[new_total_matrices - 1] = inputs[total_matrices - 1];
        }
        inputs = new_inputs;
        total_matrices = new_total_matrices;
    }
    BCSR final_result = inputs[0];

    int itersize = 1;
    
    while (itersize < size){
        if (rank%(itersize*2)){
            // send result
            int dest = rank - itersize;
            send_bcsr(dest, final_result,k);            
            break;
        }
        else{
            if (rank + itersize < size){
                // receive result
                int src = rank + itersize;
                BCSR input;
                recv_bcsr(src, input,k);
                BCSR C;
                CSR_mul(final_result, input, C, k);
                // now change the final_result
                delete[] final_result.rowptr;
                delete[] final_result.colind;
                delete[] final_result.values;
                final_result = C;
                delete[] input.rowptr;
                delete[] input.colind;
                delete[] input.values;
            }
            
        }
        itersize *= 2;
    }

    if (rank == 0){
        write_final_ans(folder_path, final_result, k);
    }
    MPI_Finalize();
    return 0;
   
}