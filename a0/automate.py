import subprocess
import os
import numpy as np
output_file = "perf_output.txt"

def gen_input_files():
    for j in range(9,14):
        # make a directory for the input files and one for output filess
        os.system(f"mkdir paths/input_path_{2**j}")
        os.system(f"mkdir paths/output_path_{2**j}")
        # generate the input files
        mtx_A = np.random.random(size = (2**j, 2**j)) * 1e2 # dtype = float64
        mtx_B = np.random.random(size = (2**j, 2**j)) * 1e2 # dtype = float64

        with open(f"paths/input_path_{2**j}/mtx_A.bin", "wb") as fp:
            fp.write(mtx_A.tobytes())
        with open(f"paths/input_path_{2**j}/mtx_B.bin", "wb") as fp:
            fp.write(mtx_B.tobytes())


def find_cache_misses():
    output_file = "perf_output_cache.txt"
    with open(output_file, "w") as f:
        for j in [2**x for x in range(4,9)]:
            for i in range(6):
                for k in range(3):      
                    command = f"perf stat -e cache-misses ./main {i} {j} {j} {j} paths/input_path_{j}/ paths/output_path_{j}/"
                    f.write(f"Command: {command}\n")
                    try:
                        result = subprocess.run(command, shell=True, capture_output=True, text=True)
                        f.write(result.stdout)
                        f.write(result.stderr)
                    except Exception as e:
                        f.write(f"Error running command: {e}\n")
                    f.write("\n")
            print(f"Done with {j}x{j} matrix multiplication with {i}th algorithm")
    


def find_cache_references():
    with open(output_file, "w") as f:
        for j in [2**x for x in range(12, 14)]:
            for i in range(6):
                for k in range(3):      
                    command = f"perf stat ./main {i} {j} {j} {j} paths/input_path_{j}/ paths/output_path_{j}/"
                    f.write(f"Command: {command}\n")
                    try:
                        result = subprocess.run(command, shell=True, capture_output=True, text=True)
                        f.write(result.stdout)
                        f.write(result.stderr)
                    except Exception as e:
                        f.write(f"Error running command: {e}\n")
                    f.write("\n")
                print(f"Done with {j}x{j} matrix multiplication with {i}th algorithm")

if __name__ == "__main__":
    #gen_input_files()
    #find_cache_references()
    find_cache_misses()
    # os.system("rm -r paths")