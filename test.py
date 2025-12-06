import subprocess

files_primal = [
    r".\example_exp_iters.txt",
    r".\example_phase1.txt",
    r".\example_phase2.txt",
    r".\example_stalling.txt",
]

files = files_primal

script_path = r".\SIMPLEX.py"

for file in files:
    print(f"Running SIMPLEX.py with file {file}...\n")
    result = subprocess.run(["python", script_path, file], capture_output=True, text=True)
    print(result.stdout)
    if result.stderr:
        print("Errors:", result.stderr)
    print("="*80)
