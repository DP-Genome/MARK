import os
import subprocess
import tkinter as tk
from tkinter import filedialog


def clean_bam_files(input_folder):
    print(f"Processing files in: {input_folder}")

    # Check if folder exists
    if not os.path.exists(input_folder):
        print(f"Error: The folder '{input_folder}' does not exist.")
        return

    for filename in os.listdir(input_folder):
        if filename.endswith(".bam") and not filename.endswith("_cleaned.bam"):
            input_bam_path = os.path.join(input_folder, filename)
            output_bam_path = os.path.join(input_folder, filename.replace(".bam", "_cleaned.bam"))

            print(f"Working on {filename}...")

            # Step 1: Extract header
            header_sam_path = os.path.join(input_folder, "temp_header.sam")
            with open(header_sam_path, "w") as header_sam:
                subprocess.run(["samtools", "view", "-H", input_bam_path], stdout=header_sam)

            # Step 2: Edit the header
            new_header_sam_path = os.path.join(input_folder, "new_header.sam")
            with open(header_sam_path, "r") as infile, open(new_header_sam_path, "w") as outfile:
                subprocess.run(
                    [
                        "sed",
                        "-E",
                        "s/_linearized//g; s/\\.1//g; /^@PG.*ID:samtools/d"
                    ],
                    stdin=infile,
                    stdout=outfile
                )

            # Step 3: Replace the header
            with open(output_bam_path, "wb") as output_bam:
                subprocess.run(["samtools", "reheader", new_header_sam_path, input_bam_path],
                               stdout=output_bam)

            # Step 4: Index the cleaned BAM file
            subprocess.run(["samtools", "index", output_bam_path])

            # Clean up temporary files
            if os.path.exists(header_sam_path):
                os.remove(header_sam_path)
            if os.path.exists(new_header_sam_path):
                os.remove(new_header_sam_path)

            print(f"  -> Created: {os.path.basename(output_bam_path)}")


def select_folder_and_run():
    # Create the root window but hide it (we only want the dialog)
    root = tk.Tk()
    root.withdraw()

    # Open the directory selection dialog
    folder_path = filedialog.askdirectory(title="Select Folder Containing BAM Files")

    # Check if a folder was actually selected (user didn't click Cancel)
    if folder_path:
        print(f"Selected folder: {folder_path}")
        clean_bam_files(folder_path)
        print("All operations complete.")
    else:
        print("No folder selected. Exiting.")


if __name__ == "__main__":
    select_folder_and_run()