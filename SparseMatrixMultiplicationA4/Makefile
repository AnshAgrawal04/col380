CC = mpic++
GCC = g++
NVCC = nvcc
CFLAGS = -O3 -Wall -fopenmp -g
CUDAFLAGS = -Xcompiler -fopenmp -O3

MPI_SRC = main.cpp
MPI_OBJ = comms.o

CUDA_SRC = cuda.cu
CUDA_OBJ = matmul.o

IO_SRC = read_write.cpp
IO_OBJ = read_write.o

TARGET = a4

all: $(TARGET)

$(MPI_OBJ): $(MPI_SRC)
	# $(CC) $(CFLAGS) -c $(MPI_SRC) -o $(MPI_OBJ)
	# $(GCC) $(CFLAGS) -c $(IO_SRC) -o $(IO_OBJ)
	$(GCC) $(CFLAGS) -c main.cpp -o main.o

$(CUDA_OBJ): $(CUDA_SRC)
	$(NVCC) $(CUDAFLAGS) -c $(CUDA_SRC) -o $(CUDA_OBJ)

$(TARGET): $(MPI_OBJ) $(CUDA_OBJ)
	$(CC) $(CFLAGS) main.o $(CUDA_OBJ) -lcudart -L/usr/local/cuda/lib64 -o $(TARGET)

run: $(TARGET)
	mpirun -np 1 ./$(TARGET) $(FOLDERNAME)

clean:
	rm -f $(MPI_OBJ) $(CUDA_OBJ) $(IO_OBJ) $(TARGET) main.o
	rm JOB_*
	rm JOB*

