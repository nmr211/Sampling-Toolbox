import subprocess

"""
This script generates samples using Latin Hypercube Sampling (LHS)
based on fitted distributions of space weather variables.

How to use:

1. Set the number of samples:
   - Modify the value of SAMPLES (e.g., 1000, 10000, 50000)

2. Select variables:
   - Edit the VARIABLES list
   - Remove any variables you do NOT want to sample
   - Available options:
     ["bx", "by", "bz", "swvel", "swden", "al", "un", "f107"]

3. Run the script
   

Output:
- The generated NetCDF file will be saved in the same directory
- Default file name: lhs_samples.nc (unless changed in settings)
"""

SAMPLES = 10000    # You can change the number of samples that you need

VARIABLES = [
    "bx",         # You can hash out which input parameters you dont want to use
    "by",
    "bz",
    "swvel",
    "swden",
    "al",
    "un",
    "f107"
]

OUTPUT_FILE = "lhs_samples.nc"   # Name of the output file 


# ============================================================================
# RUN SCRIPT (DO NOT EDIT)
# =============================================================================

def main():
    cmd = [
        "python", "toolbox.py",
        "--samples", str(SAMPLES),
        "--output", OUTPUT_FILE,
        "--variables"
    ] + VARIABLES

    print("\nRunning sampling with configuration:")
    print(f"  Samples   : {SAMPLES}")
    print(f"  Variables  : {VARIABLES}")
    print(f"  Output     : {OUTPUT_FILE}\n")

    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()