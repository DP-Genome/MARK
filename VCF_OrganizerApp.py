import os
import sys
import json
import datetime
import shutil
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import tkinter.scrolledtext as scrolledtext

class VCFOrganizerApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("VCF Organizer Tool v1.0")
        self.geometry("850x750")
        self.configure(padx=20, pady=20)

        # Apply a nice native-looking theme if available
        style = ttk.Style(self)
        if 'clam' in style.theme_names():
            style.theme_use('clam')
            
        style.configure('TFrame', background=self.cget('bg'))
        style.configure('TLabelframe', background=self.cget('bg'))
        style.configure('TLabelframe.Label', font=('Helvetica', 12, 'bold'))

        # Variables
        self.input_dir_var = tk.StringVar()
        self.output_dir_var = tk.StringVar()
        self.mapping_var = tk.StringVar()
        self.collection_name_var = tk.StringVar()
        
        # Default mapping file if it exists in the same directory
        default_mapping = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sample_mapping.json")
        if os.path.exists(default_mapping):
            self.mapping_var.set(default_mapping)

        self.suffix_var = tk.StringVar(value="_trimmed_snps_corrected.vcf")
        
        self.manual_keywords_var = tk.StringVar()
        self.manual_folder_var = tk.StringVar()

        self.create_widgets()

    def create_widgets(self):
        title_label = ttk.Label(self, text="Mitochondrial VCF Organizer", font=("Helvetica", 18, "bold"))
        title_label.pack(pady=(0, 20))

        # --- Directories Frame ---
        dir_frame = ttk.LabelFrame(self, text="1. Paths & Configuration", padding=15)
        dir_frame.pack(fill=tk.X, pady=(0, 15))

        # Config Grid
        ttk.Label(dir_frame, text="Input Directory (Contains Pipeline Outputs):").grid(row=0, column=0, sticky=tk.W, pady=8)
        ttk.Entry(dir_frame, textvariable=self.input_dir_var, width=55).grid(row=0, column=1, padx=10, pady=8)
        ttk.Button(dir_frame, text="Browse...", command=self.browse_input).grid(row=0, column=2, pady=8)

        ttk.Label(dir_frame, text="Destination Output Directory:").grid(row=1, column=0, sticky=tk.W, pady=8)
        ttk.Entry(dir_frame, textvariable=self.output_dir_var, width=55).grid(row=1, column=1, padx=10, pady=8)
        ttk.Button(dir_frame, text="Browse...", command=self.browse_output).grid(row=1, column=2, pady=8)

        ttk.Label(dir_frame, text="Sample Mapping JSON:").grid(row=2, column=0, sticky=tk.W, pady=8)
        ttk.Entry(dir_frame, textvariable=self.mapping_var, width=55).grid(row=2, column=1, padx=10, pady=8)
        ttk.Button(dir_frame, text="Browse...", command=self.browse_mapping).grid(row=2, column=2, pady=8)

        ttk.Label(dir_frame, text="Target File Suffix (.vcf):").grid(row=3, column=0, sticky=tk.W, pady=8)
        ttk.Entry(dir_frame, textvariable=self.suffix_var, width=55).grid(row=3, column=1, padx=10, pady=8, sticky=tk.W)

        ttk.Label(dir_frame, text="Master Folder Base Name (Optional):").grid(row=4, column=0, sticky=tk.W, pady=8)
        ttk.Entry(dir_frame, textvariable=self.collection_name_var, width=55).grid(row=4, column=1, padx=10, pady=8, sticky=tk.W)

        # --- Auto Organization Frame ---
        auto_frame = ttk.LabelFrame(self, text="2. Auto Organization", padding=15)
        auto_frame.pack(fill=tk.X, pady=(0, 15))
        
        ttk.Label(auto_frame, text="Automatically organize files into sample folders (007, 2800M, etc.) using the mapping JSON.").pack(side=tk.LEFT)
        ttk.Button(auto_frame, text="Run Auto Organize", command=self.run_auto).pack(side=tk.RIGHT)

        # --- Manual Extraction Frame ---
        manual_frame = ttk.LabelFrame(self, text="3. Manual Extraction / Custom Search", padding=15)
        manual_frame.pack(fill=tk.X, pady=(0, 15))

        manual_grid = ttk.Frame(manual_frame)
        manual_grid.pack(fill=tk.X)

        ttk.Label(manual_grid, text="Search Substrings (comma separated, e.g. 5, 7, 9):").grid(row=0, column=0, sticky=tk.W, pady=8)
        ttk.Entry(manual_grid, textvariable=self.manual_keywords_var, width=40).grid(row=0, column=1, padx=10, pady=8)

        ttk.Label(manual_grid, text="Target Folder Name:").grid(row=1, column=0, sticky=tk.W, pady=8)
        ttk.Entry(manual_grid, textvariable=self.manual_folder_var, width=40).grid(row=1, column=1, padx=10, pady=8)

        ttk.Button(manual_grid, text="Search & Extract", command=self.run_manual).grid(row=0, column=2, rowspan=2, padx=20, sticky=tk.NS)

        # --- Log Frame ---
        log_frame = ttk.LabelFrame(self, text="Logs & Output", padding=10)
        log_frame.pack(fill=tk.BOTH, expand=True)

        self.log_area = scrolledtext.ScrolledText(log_frame, state='disabled', wrap='word', font=("Consolas", 11))
        self.log_area.pack(fill=tk.BOTH, expand=True)
        
        self.log("VCF Organizer Application successfully initialized. Ready.")

    def log(self, message):
        self.log_area.config(state='normal')
        self.log_area.insert(tk.END, message + "\n")
        self.log_area.see(tk.END)
        self.log_area.config(state='disabled')
        if hasattr(self, 'update_idletasks'):
            self.update_idletasks()

    def browse_input(self):
        d = filedialog.askdirectory(title="Select Input Directory")
        if d: self.input_dir_var.set(d)

    def browse_output(self):
        d = filedialog.askdirectory(title="Select Output Directory")
        if d: self.output_dir_var.set(d)

    def browse_mapping(self):
        f = filedialog.askopenfilename(title="Select Mapping JSON", filetypes=[("JSON files", "*.json"), ("All files", "*.*")])
        if f: self.mapping_var.set(f)

    def get_vcf_files(self, input_dir, suffix):
        vcfs = []
        for root, dirs, files in os.walk(input_dir):
            for file in files:
                if file.endswith(suffix):
                    vcfs.append(os.path.join(root, file))
        return vcfs

    def get_collection_name(self, in_dir, base_name):
        now = datetime.datetime.now()
        timestamp = now.strftime("%d%B%Y-%H%M")
        
        if not base_name:
            # Try to deduce from immediate subdirectories
            try:
                subdirs = [d for d in os.listdir(in_dir) if os.path.isdir(os.path.join(in_dir, d))]
                if len(subdirs) >= 2:
                    s1 = str(subdirs[0])
                    longest = ""
                    for i in range(len(s1)):
                        for j in range(i+1, len(s1)+1):
                            sub = s1[i:j]
                            if all(sub in str(s) for s in subdirs):
                                if len(sub) > len(longest):
                                    longest = sub
                    base_name = str(longest).strip('_-')
                elif len(subdirs) == 1:
                    base_name = str(subdirs[0])
            except Exception:
                pass
                
        if not base_name or len(base_name) < 3:
            base_name = "VCF_Collection"
            
        return f"{base_name}_{timestamp}"

    def run_auto(self):
        in_dir = self.input_dir_var.get().strip()
        out_dir = self.output_dir_var.get().strip()
        map_file = self.mapping_var.get().strip()
        suffix = self.suffix_var.get().strip()

        if not in_dir or not out_dir or not map_file:
            messagebox.showerror("Validation Error", "Please configure the Input Directory, Output Directory, and Mapping File first.")
            return

        if not os.path.exists(map_file):
            messagebox.showerror("File Error", "The selected mapping JSON file cannot be found.")
            return

        self.log("\n" + "="*60)
        self.log("▶ Starting Auto Organization Process")
        self.log("="*60)
        
        # Use a thread so the UI does not freeze during intensive I/O
        threading.Thread(target=self.do_auto_organize, args=(in_dir, out_dir, map_file, suffix), daemon=True).start()

    def do_auto_organize(self, in_dir, out_dir, map_file, suffix):
        try:
            with open(map_file, 'r') as f:
                mapping = json.load(f)
            
            self.log(f"[*] Loaded mapping for {len(mapping)} biological samples.")
            vcf_files = self.get_vcf_files(in_dir, suffix)
            self.log(f"[*] Scanning input directory recursively...")
            self.log(f"[*] Found {len(vcf_files)} total files matching suffix '{suffix}'.\n")

            collection_name = self.get_collection_name(in_dir, self.collection_name_var.get().strip())
            master_out_dir = os.path.join(out_dir, collection_name)
            if not os.path.exists(master_out_dir):
                os.makedirs(master_out_dir)
            self.log(f"[+] Made Master Collection Folder: {collection_name}/")

            copied_count = 0
            for sample_name, barcodes in mapping.items():
                sample_out_dir = os.path.join(master_out_dir, sample_name)
                
                for vfile in vcf_files:
                    filename = os.path.basename(vfile)
                    # Check if the filename contains ANY of the allowed barcodes for this sample
                    if any((str(bc) in filename) for bc in barcodes):
                        if not os.path.exists(sample_out_dir):
                            os.makedirs(sample_out_dir)
                            self.log(f"[+] Made directory: {sample_name}/")
                        
                        dest = os.path.join(sample_out_dir, filename)
                        if not os.path.exists(dest):
                            try:
                                shutil.copy2(vfile, dest)
                                self.log(f"  └── 📄 Copied to {sample_name}/ : {filename}")
                                copied_count += 1
                            except Exception as e:
                                self.log(f"  └── ❌ Error copying {filename}: {e}")
                        else:
                            # Skip standard log to avoid flooding, unless debugging
                            pass
                            
            self.log(f"\n✅ Auto Organization Complete! Processed and grouped {copied_count} new files.")
        except Exception as e:
            self.log(f"\n❌ CRITICAL ERROR: {str(e)}")

    def run_manual(self):
        in_dir = self.input_dir_var.get().strip()
        out_dir = self.output_dir_var.get().strip()
        keywords = self.manual_keywords_var.get().strip()
        folder_name = self.manual_folder_var.get().strip()
        suffix = self.suffix_var.get().strip()

        if not in_dir or not out_dir:
            messagebox.showerror("Validation Error", "Please configure the Input Directory and Output Directory first.")
            return

        if not keywords or not folder_name:
            messagebox.showerror("Validation Error", "Please provide Search Substrings AND a Target Folder Name.")
            return

        # Parse keywords handling potential commas and spaces
        kw_list = [k.strip() for k in keywords.split(',') if k.strip()]
        if not kw_list:
            messagebox.showerror("Validation Error", "Invalid search substrings.")
            return

        self.log("\n" + "="*60)
        self.log(f"▶ Starting Manual Search & Extract")
        self.log(f"[*] Search Substrings : {kw_list}")
        self.log(f"[*] Target Folder     : {folder_name}")
        self.log("="*60)
        
        threading.Thread(target=self.do_manual_extract, args=(in_dir, out_dir, kw_list, folder_name, suffix), daemon=True).start()

    def do_manual_extract(self, in_dir, out_dir, kw_list, folder_name, suffix):
        try:
            self.log(f"[*] Scanning input directory recursively...")
            vcf_files = self.get_vcf_files(in_dir, suffix)
            self.log(f"[*] Found {len(vcf_files)} total files matching suffix '{suffix}'. Checking for keyword matches...\n")

            collection_name = self.get_collection_name(in_dir, self.collection_name_var.get().strip())
            master_out_dir = os.path.join(out_dir, collection_name)
            target_dir = os.path.join(master_out_dir, folder_name)
            copied_count = 0

            for vfile in vcf_files:
                filename = os.path.basename(vfile)
                # If ANY of the search substrings EXACTLY match within the filename
                if any(kw in filename for kw in kw_list):
                    if not os.path.exists(target_dir):
                        os.makedirs(target_dir)
                        self.log(f"[+] Made directory: {folder_name}/")
                    
                    dest = os.path.join(target_dir, filename)
                    if not os.path.exists(dest):
                        try:
                            shutil.copy2(vfile, dest)
                            
                            # Identify which matched for the log
                            matched_kw = [kw for kw in kw_list if kw in filename]
                            self.log(f"  └── 📄 Copied (Matched {matched_kw}): {filename}")
                            copied_count += 1
                        except Exception as e:
                            self.log(f"  └── ❌ Error copying {filename}: {e}")
                            
            self.log(f"\n✅ Manual Extract Complete! Copied {copied_count} new files into '{folder_name}'.")
        except Exception as e:
            self.log(f"\n❌ CRITICAL ERROR: {str(e)}")

if __name__ == "__main__":
    app = VCFOrganizerApp()
    app.mainloop()
