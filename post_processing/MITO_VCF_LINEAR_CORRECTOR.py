import os
import tkinter as tk
from tkinter import filedialog, messagebox

# It's good practice to handle potential import errors.
try:
    from Bio import SeqIO
except ImportError:
    # If BioPython is not installed, show an error and exit.
    root = tk.Tk()
    root.withdraw()
    messagebox.showerror(
        "Missing Library",
        "BioPython is not installed. Please install it by running: pip install biopython"
    )
    exit()


# --- Core Scientific Functions ---

def create_position_map(original_fasta_file, split_position=8284):
    """
    Parses the original circular genome to create a coordinate map for linearization.
    This map is essential for converting coordinates back to the original circular reference.

    Args:
        original_fasta_file (str): Path to the original FASTA file.
        split_position (int): The position where the circular genome was split.

    Returns:
        dict: A map where {new_linearized_position: original_circular_position}.
    """
    try:
        records = list(SeqIO.parse(original_fasta_file, "fasta"))
        if not records:
            raise ValueError("The FASTA file is empty or not properly formatted.")

        genome_seq = records[0].seq
        genome_length = len(genome_seq)

        # Ensure the split position is valid.
        if not (0 < split_position < genome_length):
            raise ValueError(f"Split position {split_position} is out of bounds for genome of length {genome_length}.")

        position_map = {}
        for i in range(genome_length):
            original_pos = i + 1
            if i < split_position:
                # This part was the second half of the original genome.
                linearized_pos = i + (genome_length - split_position) + 1
            else:
                # This part was the first half of the original genome.
                linearized_pos = i - split_position + 1
            position_map[linearized_pos] = original_pos

        return position_map

    except Exception as e:
        # Propagate the error to be caught by the main function.
        raise ValueError(f"Failed to parse genome '{os.path.basename(original_fasta_file)}'.\nError: {e}")


def revert_vcf_positions(vcf_file_path, position_map):
    """
    Reads a VCF file, corrects the variant positions using the map,
    and writes to a new file with a '_corrected' suffix.

    Args:
        vcf_file_path (str): Path to the input VCF file.
        position_map (dict): The map to convert positions.
    """
    base_name = os.path.splitext(os.path.basename(vcf_file_path))[0]
    output_vcf_path = os.path.join(os.path.dirname(vcf_file_path), f"{base_name}_corrected.vcf")

    with open(vcf_file_path, "r") as f_in, open(output_vcf_path, "w") as f_out:
        for line in f_in:
            if line.startswith("#"):
                f_out.write(line)  # Write header lines directly
            else:
                parts = line.strip().split("\t")
                if len(parts) > 1 and parts[1].isdigit():
                    linearized_pos = int(parts[1])
                    # Look up the original position in the map.
                    # If not found, keep the original (as a fallback).
                    original_pos = position_map.get(linearized_pos, linearized_pos)
                    parts[1] = str(original_pos)
                    f_out.write("\t".join(parts) + "\n")
                else:
                    f_out.write(line)  # Write malformed data lines as-is


# --- Main Application Logic ---

def main():
    """
    Main function to run the VCF correction workflow.
    """
    # Create a root window but hide it. This is key to making file dialogs work well.
    root = tk.Tk()
    root.withdraw()

    # 1. Get the MODIFIED (linearized) reference genome.
    messagebox.showinfo(
        "Step 1: Select Modified Genome",
        "Please select the MODIFIED (linearized) reference genome file."
    )
    modified_fasta = filedialog.askopenfilename(
        title="Select Modified (Linearized) Reference Genome",
        filetypes=[("FASTA Files", "*.fasta *.fna *.fa"), ("All Files", "*.*")]
    )
    if not modified_fasta:
        messagebox.showwarning("Cancelled", "Operation cancelled.")
        return

    # Note: The 'modified_fasta' variable is captured as requested but not used by the
    # 'create_position_map' or 'revert_vcf_positions' functions that follow.

    # 2. Get the ORIGINAL circular reference genome to build the coordinate map.
    messagebox.showinfo(
        "Step 2: Select Original Genome",
        "Next, please select the ORIGINAL circular mitochondrial reference genome file."
    )
    original_fasta = filedialog.askopenfilename(
        title="Select Original Circular Reference Genome",
        filetypes=[("FASTA Files", "*.fasta *.fna *.fa"), ("All Files", "*.*")]
    )
    if not original_fasta:
        messagebox.showwarning("Cancelled", "Operation cancelled.")
        return

    # 3. Create the position map from the ORIGINAL genome.
    try:
        position_map = create_position_map(original_fasta)
    except Exception as e:
        messagebox.showerror("Genome Error", str(e))
        return

    # 4. Get the folder containing the VCF files.
    messagebox.showinfo(
        "Step 3: Select VCF Folder",
        "Finally, please select the folder containing the VCF files you want to correct."
    )
    vcf_folder = filedialog.askdirectory(title="Select VCF Folder")
    if not vcf_folder:
        messagebox.showwarning("Cancelled", "Operation cancelled.")
        return

    # 5. Process each VCF file.
    processed_files = []
    failed_files = []
    vcf_files_to_process = [f for f in os.listdir(vcf_folder) if f.lower().endswith(".vcf")]

    if not vcf_files_to_process:
        messagebox.showinfo("No Files Found", "No '.vcf' files were found in the selected directory.")
        return

    for filename in vcf_files_to_process:
        vcf_file_path = os.path.join(vcf_folder, filename)
        try:
            revert_vcf_positions(vcf_file_path, position_map)
            processed_files.append(filename)
        except Exception as e:
            failed_files.append(f"{filename} (Error: {e})")

    # 6. Show a final report.
    summary_message = f"Processing Complete!\n\n"
    summary_message += f"✅ Successfully processed: {len(processed_files)} files.\n"
    if processed_files:
        summary_message += f"New files have been saved with a '_corrected.vcf' suffix.\n\n"

    if failed_files:
        summary_message += f"❌ Failed to process: {len(failed_files)} files.\n\n"
        summary_message += "Failed Files:\n" + "\n".join(failed_files)

    messagebox.showinfo("Summary", summary_message)


if __name__ == "__main__":
    main()