import subprocess
import os
import numpy as np
output_file = "perf_output.txt"

def gen_input_files():
    for j in range(1,9):
        # make a directory for the input files and one for output filess
        os.system(f"mkdir paths/input_path_{1000*j}")
        os.system(f"mkdir paths/output_path_{1000*j}")
        # generate the input files
        mtx_A = np.random.random(size = (1000*j,1000*j)) * 1e2 # dtype = float64
        mtx_B = np.random.random(size = (1000*j,1000*j)) * 1e2 # dtype = float64

        with open(f"paths/input_path_{1000*j}/mtx_A.bin", "wb") as fp:
            fp.write(mtx_A.tobytes())
        with open(f"paths/input_path_{1000*j}/mtx_B.bin", "wb") as fp:
            fp.write(mtx_B.tobytes())


def find_cache_misses():
    output_file = "perf_cache_1000.txt"
    with open(output_file, "w") as f:
        for j in [1000*x for x in range(1,9)]:
            for i in range(6):
                for k in range(3):      
                    command = f"perf stat -e cycles,instructions,task-clock,branch-misses,cache-misses,instructions,cycles,cache-references,cache-misses ./main {i} {j} {j} {j} paths/input_path_{j}/ paths/output_path_{j}/"
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
    output_file = "perf_cpu_1000.txt"
    with open(output_file, "w") as f:
        for j in [1000*x for x in range(1,9)]:
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

def run_gprof():
    for j in [1000*x for x in range(1,6)]:
        for i in range(6):
            for k in range(3):
                newfile = f"cm_outputs/cache_misses_{j}_{i}_{k}.txt"
                with open(newfile, "w") as f:
                    command = f"perf stat -e cycles,instructions,task-clock,branch-misses,cache-misses,instructions,cycles,cache-references,cache-misses ./main {i} {j} {j} {j} paths/input_path_{j}/ paths/output_path_{j}/"
                    f.write(f"Command: {command}\n")
                    try:
                        result = subprocess.run(command, shell=True, capture_output=True, text=True)
                        f.write(result.stdout)
                        f.write(result.stderr)
                    except Exception as e:
                        f.write(f"Error running command: {e}\n")
                    f.write("\n")
                command = "ls -sh gmon.out"
                result = subprocess.run(command, shell=True, capture_output=True, text=True)
                
                gprof_file = f"gprof_outputs/gprof_{j}_{i}_{k}.txt"
                command = f"gprof --graph ./main  > {gprof_file}"
                os.system(command)
                print(f"Done with {j}x{j} matrix multiplication with {i}th algorithm")


if __name__ == "__main__":
    #gen_input_files()
    #find_cache_references()
    #find_cache_misses()
    # os.system("rm -r paths")
    run_gprof()