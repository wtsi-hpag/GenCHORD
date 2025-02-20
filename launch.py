import os
import subprocess
from itertools import islice

# CONFIGURATION
ROOT = "/lustre/scratch126/casm/team267ms/sb71/project_dft_low_cov/data/devils_bam_7261/"
LOWCOV_FILE = f"lowcovindex.dat"  # Updated source file

N = 15  # Number of parallel screen sessions
COMMAND_TEMPLATE = "echo `./genchord -file " + ROOT + "{filename}.bam -output /Archives/LowQuality/{filename} -process 1'"


def read_file(filename):
    """Read the file and extract filenames."""
    with open(filename, "r") as f:
        lines = f.readlines()
    return [line.strip() for line in lines if line.strip()]


def chunk_list(data, n):
    """Yield n roughly equal parts of the data."""
    it = iter(data)
    return [list(islice(it, len(data) // n + (1 if i < len(data) % n else 0))) for i in range(n)]


def create_screen_script(index, files):
    """Create a shell script to run commands in a screen session."""
    script_name = f"screen_task_{index}.sh"
    with open(script_name, "w") as f:
        f.write("#!/bin/bash\n")
        for filename in files:
            command = COMMAND_TEMPLATE.format(filename=filename)
            f.write(f"{command} \n")
    os.chmod(script_name, 0o755)
    return script_name


def launch_screens(scripts):
    """Launch screen sessions running the generated scripts."""
    for i, script in enumerate(scripts):
        screen_name = f"cancer_task_{i}"
        subprocess.run(["screen", "-dmS", screen_name, "bash", "-c", f"./{script}; exec bash"])
        print(f"Launched screen: {screen_name} running {script}")


def main():
    file_list = read_file(LOWCOV_FILE)
    chunks = chunk_list(file_list, N)
    scripts = [create_screen_script(i, chunk) for i, chunk in enumerate(chunks)]
    launch_screens(scripts)


if __name__ == "__main__":
    main()
