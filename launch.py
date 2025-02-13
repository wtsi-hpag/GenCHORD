import os
import subprocess
from itertools import islice

# CONFIGURATION
ROOT = "../../zn1/ed15/Devil/Data_Upload/BAMs/cancer/"
SH_MOVE_FILE = f"{ROOT}sh.move"

N = 15  # Number of parallel screen sessions
COMMAND_TEMPLATE = "echo './genchord -f " + ROOT + "{filename}.bam -output Archives/HighQuality/{cancer_type}/{filename} -process 1'"  # Replace with your command

def read_file(filename):
    """Read the file and parse mv commands into (filename, cancer_type) pairs."""
    with open(filename, "r") as f:
        lines = f.readlines()
    parsed = [line.strip().split()[1:] for line in lines if line.startswith("mv ")]
    return [(file, cancer) for file, cancer in parsed]

def chunk_list(data, n):
    """Yield n roughly equal parts of the data."""
    it = iter(data)
    return [list(islice(it, len(data) // n + (1 if i < len(data) % n else 0))) for i in range(n)]

def create_screen_script(index, files):
    """Create a shell script to run commands in a screen session."""
    script_name = f"screen_task_{index}.sh"
    with open(script_name, "w") as f:
        f.write("#!/bin/bash\n")
        # f.write("cd /file/location\n")  # Change to the appropriate directory
        for filename, cancer_type in files:
            command = COMMAND_TEMPLATE.format(filename=filename, cancer_type=cancer_type)
            f.write(f"{command} \n")  # Run each command in the background
        # f.write("wait\n")  # Wait for all commands qto finish before exiting
    os.chmod(script_name, 0o755)
    return script_name

def launch_screens(scripts):
    """Launch screen sessions running the generated scripts."""
    for i, script in enumerate(scripts):
        screen_name = f"cancer_task_{i}"
        subprocess.run(["screen", "-dmS", screen_name, "bash", "-c", f"./{script}; exec bash"])
        print(f"Launched screen: {screen_name} running {script}")

def main():
    file_list = read_file(SH_MOVE_FILE)
    chunks = chunk_list(file_list, N)
    scripts = [create_screen_script(i, chunk) for i, chunk in enumerate(chunks)]
    launch_screens(scripts)

if __name__ == "__main__":
    main()
