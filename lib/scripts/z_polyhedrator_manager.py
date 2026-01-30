import os
import sys
import subprocess
from script_colors import SCRIPT_TAG, ERROR_TAG

class ZPolyhedratorManager:
    """Manages the cloning, building, and running of the z_polyhedrator project."""

    def __init__(self, repo_dir):
        """
        Initialize the manager with a configurable repository directory and URL.
        :param repo_dir: Directory where the repository should be cloned.
        :param repo_url: GitHub repository URL to clone.
        """
        self.Z_POLYHEDRAL_DIR = repo_dir
        self.TARGET_DIR = "target/release/z_polyhedrator"

    def generate_uzp(self, pattern_file, input_matrix_file, output_file_name):
        """Generate the initial UZP using z_polyhedrator, ensuring correct directory execution."""

        print(f"{SCRIPT_TAG} Generating UZP with {input_matrix_file} and {pattern_file}")

        # get RUSTUP_HOME and CARGO_HOME from environment variables if set
        rustup_home = os.getenv("RUSTUP_HOME", "/opt/rustup")
        cargo_home = os.getenv("CARGO_HOME", "/tmp/.cargo_pldi25_artifact")

        # Construct the command to execute everything in a single shell command
        command = f"""
        cd {self.Z_POLYHEDRAL_DIR} > /dev/null && \
        RUSTUP_HOME={rustup_home} CARGO_HOME={cargo_home} rustup install 1.85.0
        RUSTUP_HOME={rustup_home} CARGO_HOME={cargo_home} RUSTFLAGS="-C opt-level=3 -C target-cpu=native" cargo build --release && \
        {self.TARGET_DIR} search {pattern_file} {input_matrix_file} -w {output_file_name} && \
        cd - > /dev/null
        """

        try:
            subprocess.run(command, shell=True, check=True, executable="bash")
            print(f"{SCRIPT_TAG} UZP generation completed successfully")
        except subprocess.CalledProcessError as e:
            print(f"{ERROR_TAG} UZP generation failed: {e}")


if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python z_polyhedrator_manager.py <input_matrix_file> <output_directory>")
        sys.exit(1)
        
    # get UZP_PROJECT_ROOT from environment variable if set
    project_root = os.getenv("UZP_PROJECT_ROOT", "/uzp-artifact")

    # Pass parameters dynamically
    repo_directory = os.path.join(project_root, "z_polyhedrator")

    manager = ZPolyhedratorManager(repo_directory)

    input_matrix_file = os.path.abspath(sys.argv[1])
    output_directory = os.path.abspath(sys.argv[2])
    output_file_name = os.path.join(output_directory, os.path.basename(input_matrix_file).replace('.mtx', ''))
    manager.generate_uzp(f"{project_root}/z_polyhedrator/data/patterns.txt", input_matrix_file, output_file_name)
