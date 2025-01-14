using namespace std;
#include <fstream>
#include <stdio.h>
// #include <iostream>
#define rep(a, b) for (int a = 0; a < b; a++)


void matrixMultiplyIJK(double **matrix_A, double **matrix_B, double **matrix_C, int number_row1, int number_col1, int number_col2){
    rep(i,number_row1){rep(j,number_col2){rep(k,number_col1){matrix_C[i][j] += matrix_A[i][k]*matrix_B[k][j];}}}
}
void matrixMultiplyIKJ(double **matrix_A, double **matrix_B, double **matrix_C, int number_row1, int number_col1, int number_col2){
    rep(i,number_row1){rep(k,number_col1){rep(j,number_col2){matrix_C[i][j] += matrix_A[i][k]*matrix_B[k][j];}}}
}
void matrixMultiplyJIK(double **matrix_A, double **matrix_B, double **matrix_C, int number_row1, int number_col1, int number_col2){
    rep(j,number_col2){rep(i,number_row1){rep(k,number_col1){matrix_C[i][j] += matrix_A[i][k]*matrix_B[k][j];}}}
}
void matrixMultiplyJKI(double **matrix_A, double **matrix_B, double **matrix_C, int number_row1, int number_col1, int number_col2){
    rep(j,number_col2){rep(k,number_col1){rep(i,number_row1){matrix_C[i][j] += matrix_A[i][k]*matrix_B[k][j];}}}
}
void matrixMultiplyKIJ(double **matrix_A, double **matrix_B, double **matrix_C, int number_row1, int number_col1, int number_col2){
    rep(k,number_col1){rep(i,number_row1){rep(j,number_col2){matrix_C[i][j] += matrix_A[i][k]*matrix_B[k][j];}}}
}
void matrixMultiplyKJI(double **matrix_A, double **matrix_B, double **matrix_C, int number_row1, int number_col1, int number_col2){
    rep(k,number_col1){rep(j,number_col2){rep(i,number_row1){matrix_C[i][j] += matrix_A[i][k]*matrix_B[k][j];}}}
}

int main(int argc, char *argv[])
{
    if (argc != 7)
    {
        // cout << "Invalid arguments" << endl;
        return 1;
    }
    int type = atoi(argv[1]);
    int number_row1 = atoi(argv[2]);
    int number_col1 = atoi(argv[3]);
    int number_col2 = atoi(argv[4]);

    string path_input = argv[5];
    string path_output = argv[6];

    string fileA = path_input + "/mtx_A.bin";
    string fileB = path_input + "/mtx_B.bin";

    double **matrix_A = new double *[number_row1];
    FILE *fp = fopen(fileA.c_str(), "rb");
    for (int i = 0; i < number_row1; i++)
    {
        matrix_A[i] = new double[number_col1];
        fread(matrix_A[i], sizeof(double), number_col1, fp);
    }
    fclose(fp);

    double **matrix_B = new double *[number_col1];
    fp = fopen(fileB.c_str(), "rb");
    for (int i = 0; i < number_col1; i++)
    {
        matrix_B[i] = new double[number_col2];
        fread(matrix_B[i], sizeof(double), number_col2, fp);
    }
    fclose(fp);

    // make matrix C ussing the number of rows and columns and initializing the values to 0
    double **matrix_C = new double *[number_row1];
    for (int i = 0; i < number_row1; i++)
    {
        matrix_C[i] = new double[number_col2];
        for (int j = 0; j < number_col2; j++)
        {
            matrix_C[i][j] = 0;
        }
    }

    // multiply the matrices using the type
    int repeat = 1;
    while (repeat--){
    switch (type)
    {
    case 0:
        matrixMultiplyIJK(matrix_A, matrix_B, matrix_C, number_row1, number_col1, number_col2);
        break;
    case 1:
        matrixMultiplyIKJ(matrix_A, matrix_B, matrix_C, number_row1, number_col1, number_col2);
        break;
    case 2:
        matrixMultiplyJIK(matrix_A, matrix_B, matrix_C, number_row1, number_col1, number_col2);
        break;
    case 3:
        matrixMultiplyJKI(matrix_A, matrix_B, matrix_C, number_row1, number_col1, number_col2);
        break;
    case 4:
        matrixMultiplyKIJ(matrix_A, matrix_B, matrix_C, number_row1, number_col1, number_col2);
        break;
    case 5:
        matrixMultiplyKJI(matrix_A, matrix_B, matrix_C, number_row1, number_col1, number_col2);
        break;
    }
    }
    // write the output file
    fp = fopen((path_output + "/mtx_C.bin").c_str(), "wb");
    for (int i = 0; i < number_row1; i++)
    {
        fwrite(matrix_C[i], sizeof(double), number_col2, fp);
    }
    fclose(fp);


    // free the memory
    for (int i = 0; i < number_row1; i++)
    {
        delete[] matrix_A[i];
    }
    delete[] matrix_A;
    for (int i = 0; i < number_col1; i++)
    {
        delete[] matrix_B[i];
    }   
    delete[] matrix_B;
    for (int i = 0; i < number_row1; i++)
    {
        delete[] matrix_C[i];
    }
    delete[] matrix_C;
    

}