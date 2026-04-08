import os
import sys
import shutil
import subprocess
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import threading
import re
import signal

try:
    from Bio import SeqIO
except ImportError:
    messagebox.showerror("Missing Library", "Biopython is not installed.\nPlease run: pip install biopython")
    sys.exit(1)


class MitoPipelineDashboard:
    def __init__(self, root):
        self.root = root
        self.root.title("ONT/Illumina Mitochondrial Pipeline Dashboard v4.7")
        self.root.geometry("1150x1050")

        self.script_dir = os.path.dirname(os.path.abspath(__file__))

        # State management for the new STOP feature
        self.current_process = None
        self.stop_requested = False

        self.dash_defaults = {
            'MAX_LEN': '1500',
            'DISCARD_WARN_PCT': '5',
            'QS_MIN': '10',
            'SAFETY_QUAL': '20',
            'STRICT_QUAL': '60',
            'PILEUP_MAX_DEPTH': '1000000',
            'MIN_DEPTH': '15',
            'MIN_LEN': '90',
            'MIN_LEN_POST': '90',
            'MAX_LEN_POST': '300',
            'EXTRA_TRIM': '22',
            'READQ': '20',
            'UNQUAL_PCT': '40',
        }

        self.var_states = {}
        for var, default_val in self.dash_defaults.items():
            self.var_states[var] = {
                'entry_var': tk.StringVar(value=default_val),
                'override_var': tk.BooleanVar(value=False),
                'entry_widget': None,
                'chk_widget': None,
                'final_label': None
            }

        style = ttk.Style()

        self.notebook = ttk.Notebook(root)
        self.notebook.pack(expand=True, fill='both', padx=10, pady=10)

        self.tab_pipeline = ttk.Frame(self.notebook)
        self.tab_post = ttk.Frame(self.notebook)

        self.notebook.add(self.tab_pipeline, text="Step 1: Run Pipeline")
        self.notebook.add(self.tab_post, text="Step 2: Post-Processing")

        self.init_pipeline_tab()
        self.init_postproc_tab()

        self.log_frame = ttk.LabelFrame(root, text="Activity Log")
        self.log_frame.pack(fill='both', expand=True, padx=10, pady=(0, 10))
        self.log_text = scrolledtext.ScrolledText(self.log_frame, height=12, state='disabled', font=("Consolas", 9))
        self.log_text.pack(fill='both', expand=True, padx=5, pady=5)

        self.refresh_config_ui()

    def log(self, message):
        self.root.after(0, self._log_internal, message)
        print(message)

    def _log_internal(self, message):
        self.log_text.config(state='normal')
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)
        self.log_text.config(state='disabled')

    def safe_show_info(self, title, message):
        self.root.after(0, lambda: messagebox.showinfo(title, message))

    def safe_show_error(self, title, message):
        self.root.after(0, lambda: messagebox.showerror(title, message))

    def update_progress_ui(self, current, total, text_prefix="Processing", custom_text=None):
        def _update():
            if total > 0:
                self.progress_bar['maximum'] = total
                self.progress_bar['value'] = current

                if custom_text:
                    self.progress_var.set(custom_text)
                elif current == 0:
                    self.progress_var.set(f"Initializing... Found {total} samples.")
                elif current >= total:
                    if text_prefix == "Finished":
                        self.progress_var.set("Complete!")
                    else:
                        self.progress_var.set("Finishing up...")
                else:
                    display_num = current if current > 0 else 1
                    self.progress_var.set(f"{text_prefix}: Sample {display_num} of {total}")
            else:
                self.progress_bar.stop()
                self.progress_var.set(custom_text if custom_text else f"{text_prefix}...")

        self.root.after(0, _update)

    def set_ui_state(self, is_running):
        """Disables run buttons and enables stop button while active."""

        def _update():
            run_state = 'disabled' if is_running else 'normal'
            stop_state = 'normal' if is_running else 'disabled'

            self.btn_run_only.config(state=run_state)
            self.btn_run_auto.config(state=run_state)
            self.btn_stop.config(state=stop_state)

        self.root.after(0, _update)

    def count_input_files(self, folder):
        if not os.path.exists(folder): return 0

        r1_count = 0
        se_count = 0
        valid_subdirs = 0

        for f in os.listdir(folder):
            path = os.path.join(folder, f)
            if os.path.isdir(path):
                # Count if it matches known patterns or actually contains FASTQ files
                if f.startswith("barcode"):
                    valid_subdirs += 1
                elif any(file.lower().endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')) for file in os.listdir(path)):
                    valid_subdirs += 1
                continue

            f_lower = f.lower()
            if f_lower.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
                if '_r1' in f_lower or '.r1' in f_lower or re.search(r'_1\.f(?:ast)?q', f_lower):
                    r1_count += 1
                elif '_r2' in f_lower or '.r2' in f_lower or re.search(r'_2\.f(?:ast)?q', f_lower):
                    continue
                else:
                    se_count += 1

        total_samples = max(r1_count + se_count, valid_subdirs)
        return total_samples if total_samples > 0 else 1

    def get_def_path(self, filename):
        return os.path.join(self.script_dir, filename)

    def parse_script_defaults(self, script_path):
        script_defaults = {}
        if not os.path.exists(script_path): return script_defaults
        try:
            with open(script_path, 'r', encoding='utf-8') as f:
                content = f.read()
            for var in self.dash_defaults.keys():
                pattern_fallback = rf'{var}="?\$\{{{var}:-(?P<val>[^}}]+)\}}"?\s*'
                pattern_hardcoded = rf'^\s*{var}="?(?P<val>[^"\n$]+)"?\s*(?:#.*)?$'
                match_fall = re.search(pattern_fallback, content)
                if match_fall:
                    script_defaults[var] = match_fall.group('val')
                else:
                    match_hard = re.search(pattern_hardcoded, content, re.MULTILINE)
                    if match_hard:
                        script_defaults[var] = match_hard.group('val') + " (Hardcoded)"
        except Exception as e:
            self.log(f"Error parsing script defaults: {e}")
        return script_defaults

    # =========================================================================
    # TAB 1: MAIN PIPELINE
    # =========================================================================
    def init_pipeline_tab(self):
        script_frame = ttk.LabelFrame(self.tab_pipeline, text="1. Select Pipeline Script")
        script_frame.pack(fill='x', padx=10, pady=(10, 5))
        script_frame.columnconfigure(1, weight=1)

        ttk.Label(script_frame, text="Pipeline Script (.sh):").grid(row=0, column=0, sticky='e', padx=5, pady=10)
        self.script_path = tk.StringVar(value=self.get_def_path("ONT_MITO_v36.sh"))
        ttk.Entry(script_frame, textvariable=self.script_path).grid(row=0, column=1, sticky='ew', padx=5)
        btn_browse_script = ttk.Button(script_frame, text="Browse & Load", command=self.browse_script)
        btn_browse_script.grid(row=0, column=2, padx=5)

        self.config_frame = ttk.LabelFrame(self.tab_pipeline, text="2. Parameter Configuration")
        self.config_frame.pack(fill='x', padx=10, pady=5)

        btn_frame = ttk.Frame(self.config_frame)
        btn_frame.pack(fill='x', padx=5, pady=5)
        ttk.Button(btn_frame, text="Reset to Script Defaults", command=self.set_all_script_defaults).pack(side='left',
                                                                                                          padx=5)

        self.grid_frame = ttk.Frame(self.config_frame)
        self.grid_frame.pack(fill='x', padx=5, pady=5)

        io_frame = ttk.LabelFrame(self.tab_pipeline, text="3. Input & Output Selection")
        io_frame.pack(fill='x', padx=10, pady=5)
        io_frame.columnconfigure(1, weight=1)

        ttk.Label(io_frame, text="Linear Ref (FASTA):").grid(row=0, column=0, sticky='e', padx=5, pady=3)
        self.ref_path = tk.StringVar(value=self.get_def_path("linearized_mtdna.fasta"))
        ttk.Entry(io_frame, textvariable=self.ref_path).grid(row=0, column=1, sticky='ew', padx=5)
        ttk.Button(io_frame, text="Browse",
                   command=lambda: self.browse_file(self.ref_path, [("FASTA", "*.fasta *.fa")])).grid(row=0, column=2,
                                                                                                      padx=5)

        ttk.Label(io_frame, text="Regions BED:").grid(row=1, column=0, sticky='e', padx=5, pady=3)
        self.bed_path = tk.StringVar(value=self.get_def_path("linearized_regions.bed"))
        ttk.Entry(io_frame, textvariable=self.bed_path).grid(row=1, column=1, sticky='ew', padx=5)
        ttk.Button(io_frame, text="Browse", command=lambda: self.browse_file(self.bed_path, [("BED", "*.bed")])).grid(
            row=1, column=2, padx=5)

        ttk.Label(io_frame, text="Adapter File:").grid(row=2, column=0, sticky='e', padx=5, pady=3)
        self.adapter_path = tk.StringVar(value=self.get_def_path("Updated_Adapter_Primer_List_Cutadapt_cleaned.txt"))
        ttk.Entry(io_frame, textvariable=self.adapter_path).grid(row=2, column=1, sticky='ew', padx=5)
        ttk.Button(io_frame, text="Browse",
                   command=lambda: self.browse_file(self.adapter_path, [("Text", "*.txt")])).grid(row=2, column=2,
                                                                                                  padx=5)

        ttk.Label(io_frame, text="Input FASTQ Folder:").grid(row=3, column=0, sticky='e', padx=5, pady=3)
        self.input_folder = tk.StringVar()
        ttk.Entry(io_frame, textvariable=self.input_folder).grid(row=3, column=1, sticky='ew', padx=5)
        ttk.Button(io_frame, text="Browse", command=lambda: self.browse_folder(self.input_folder)).grid(row=3, column=2,
                                                                                                        padx=5)

        ttk.Label(io_frame, text="Output Base Dir:").grid(row=4, column=0, sticky='e', padx=5, pady=(5, 0))
        self.output_folder = tk.StringVar()
        ttk.Entry(io_frame, textvariable=self.output_folder).grid(row=4, column=1, sticky='ew', padx=5, pady=(5, 0))
        ttk.Button(io_frame, text="Browse", command=lambda: self.browse_folder(self.output_folder)).grid(row=4,
                                                                                                         column=2,
                                                                                                         padx=5,
                                                                                                         pady=(5, 0))

        ttk.Label(io_frame, text="(Optional: If left empty, outputs will save next to the input folder)",
                  font=("Arial", 8, "italic"), foreground="gray").grid(row=5, column=1, sticky='w', padx=5, pady=(0, 5))

        ttk.Label(io_frame, text="Custom Run Name:").grid(row=6, column=0, sticky='e', padx=5, pady=(5, 0))
        self.custom_run_name = tk.StringVar()
        ttk.Entry(io_frame, textvariable=self.custom_run_name).grid(row=6, column=1, sticky='ew', padx=5, pady=(5, 0))

        ttk.Label(io_frame, text="(Optional: If left empty, the pipeline uses its default naming convention)",
                  font=("Arial", 8, "italic"), foreground="gray").grid(row=7, column=1, sticky='w', padx=5, pady=(0, 5))

        ttk.Button(io_frame, text="Clear Fields", command=self.clear_pipeline_fields).grid(row=8, column=1, pady=10,
                                                                                           sticky='w', padx=5)

        # Run Buttons with direct references to toggle state
        run_btn_frame = ttk.Frame(self.tab_pipeline)
        run_btn_frame.pack(pady=10)

        self.btn_run_only = ttk.Button(run_btn_frame, text="RUN PIPELINE ONLY",
                                       command=lambda: self.start_thread(self.run_pipeline_wrapper))
        self.btn_run_only.pack(side='left', padx=10, ipadx=10, ipady=10)

        self.btn_run_auto = ttk.Button(run_btn_frame, text="RUN PIPELINE & POST-PROCESSING (AUTO)",
                                       command=lambda: self.start_thread(self.run_automated_cycle))
        self.btn_run_auto.pack(side='left', padx=10, ipadx=10, ipady=10)

        # NEW STOP BUTTON
        self.btn_stop = ttk.Button(run_btn_frame, text="STOP RUN", command=self.stop_pipeline, state='disabled')
        self.btn_stop.pack(side='left', padx=10, ipadx=10, ipady=10)

        self.progress_frame = ttk.LabelFrame(self.tab_pipeline, text="Progress Status")
        self.progress_frame.pack(fill='x', padx=10, pady=5)
        self.progress_bar = ttk.Progressbar(self.progress_frame, orient='horizontal', mode='determinate', length=400)
        self.progress_bar.pack(fill='x', padx=10, pady=10)
        self.progress_var = tk.StringVar(value="Idle")
        ttk.Label(self.progress_frame, textvariable=self.progress_var, font=("Arial", 10, "bold")).pack(pady=(0, 10))

    def clear_pipeline_fields(self):
        self.input_folder.set("")
        self.output_folder.set("")
        self.custom_run_name.set("")
        self.log("Pipeline input fields cleared.")

    def browse_script(self):
        path = filedialog.askopenfilename(filetypes=[("Shell Script", "*.sh")])
        if path:
            self.script_path.set(path)
            self.refresh_config_ui()
            self.log(f"Loaded parameters from: {os.path.basename(path)}")

    def refresh_config_ui(self):
        for widget in self.grid_frame.winfo_children():
            widget.destroy()

        script_defaults = self.parse_script_defaults(self.script_path.get())

        headers = ["Parameter", "Script Default", "Your Manual Value", "Override?", "Final Value to Pass"]
        for col, text in enumerate(headers):
            ttk.Label(self.grid_frame, text=text, font=("Arial", 9, "bold")).grid(row=0, column=col, padx=10, pady=5,
                                                                                  sticky='w')

        row = 1
        for var, dash_def in self.dash_defaults.items():
            script_val = script_defaults.get(var, "Not Found")
            is_hardcoded = "(Hardcoded)" in script_val
            clean_script_val = script_val.replace(" (Hardcoded)", "")

            ttk.Label(self.grid_frame, text=f"{var}:", font=("Arial", 9, "bold")).grid(row=row, column=0, sticky='e',
                                                                                       padx=5)

            color = "gray" if is_hardcoded else ("#1565c0" if clean_script_val != "Not Found" else "#9e9e9e")
            ttk.Label(self.grid_frame, text=script_val, foreground=color, font=("Arial", 9, "bold")).grid(row=row,
                                                                                                          column=1,
                                                                                                          padx=5,
                                                                                                          sticky='w')

            entry = ttk.Entry(self.grid_frame, textvariable=self.var_states[var]['entry_var'], width=10)
            entry.grid(row=row, column=2, padx=5, sticky='w')
            self.var_states[var]['entry_widget'] = entry
            self.var_states[var]['entry_var'].trace_add("write", lambda *args: self.update_inline_summary())

            chk = ttk.Checkbutton(self.grid_frame, variable=self.var_states[var]['override_var'],
                                  command=lambda v=var: [self.toggle_entry_state(v), self.update_inline_summary()])
            chk.grid(row=row, column=3, padx=5, sticky='w')
            self.var_states[var]['chk_widget'] = chk

            final_lbl = ttk.Label(self.grid_frame, text="", font=("Consolas", 10, "bold"))
            final_lbl.grid(row=row, column=4, padx=15, sticky='w')
            self.var_states[var]['final_label'] = final_lbl

            if clean_script_val == "Not Found" or is_hardcoded:
                self.var_states[var]['override_var'].set(False)
                chk.config(state='disabled')
                entry.config(state='disabled')
            else:
                self.toggle_entry_state(var)

            row += 1

        self.update_inline_summary()

    def update_inline_summary(self):
        script_defaults = self.parse_script_defaults(self.script_path.get())
        for var, state in self.var_states.items():
            s_val_raw = script_defaults.get(var, "Not Found")
            clean_s_val = s_val_raw.replace(" (Hardcoded)", "")
            val = state['entry_var'].get()

            if clean_s_val == "Not Found":
                state['final_label'].config(text="Disabled (Not found in script)", foreground="#9e9e9e")
            elif state['override_var'].get():
                state['final_label'].config(text=f"➜ {val} (Custom Override)", foreground="#c62828")
            else:
                state['final_label'].config(text=f"➜ {clean_s_val} (Script Default)", foreground="#1565c0")

    def toggle_entry_state(self, var):
        state = self.var_states[var]
        if state['override_var'].get():
            state['entry_widget'].config(state='normal')
        else:
            state['entry_widget'].config(state='disabled')

    def set_all_script_defaults(self):
        for var, state in self.var_states.items():
            if str(state['chk_widget']['state']) != 'disabled':
                state['override_var'].set(False)
                self.toggle_entry_state(var)
        self.update_inline_summary()

    # =========================================================================
    # TAB 2: POST-PROCESSING
    # =========================================================================
    def init_postproc_tab(self):
        # 1. Organize Raw Output
        frame_raw = ttk.LabelFrame(self.tab_post, text="1. Organize Raw Output")
        frame_raw.pack(fill='x', padx=10, pady=5)

        ttk.Label(frame_raw, text="Run Folder:").pack(side='left', padx=5)
        self.raw_main_dir = tk.StringVar()
        ttk.Entry(frame_raw, textvariable=self.raw_main_dir).pack(side='left', fill='x', expand=True, padx=5)
        ttk.Button(frame_raw, text="Browse", command=lambda: self.browse_folder(self.raw_main_dir)).pack(side='left',
                                                                                                         padx=5)
        ttk.Button(frame_raw, text="Run Step 1", command=lambda: self.start_thread(self.collect_raw_wrapper)).pack(
            side='left', padx=5)

        # 2. Linear VCF Correction
        frame_vcf = ttk.LabelFrame(self.tab_post, text="2. Linear VCF Correction")
        frame_vcf.pack(fill='x', padx=10, pady=5)
        frame_vcf.columnconfigure(1, weight=1)

        self.vcf_mod_ref = tk.StringVar(value=self.get_def_path("linearized_mtdna.fasta"))
        self.vcf_orig_ref = tk.StringVar(value=self.get_def_path("rCRS.fasta"))
        self.vcf_target_dir = tk.StringVar()

        ttk.Label(frame_vcf, text="Linearized Ref:").grid(row=0, column=0, sticky='e', padx=5, pady=3)
        ttk.Entry(frame_vcf, textvariable=self.vcf_mod_ref).grid(row=0, column=1, sticky='ew', padx=5)
        ttk.Button(frame_vcf, text="Browse",
                   command=lambda: self.browse_file(self.vcf_mod_ref, [("FASTA", "*.fasta *.fa")])).grid(row=0,
                                                                                                         column=2,
                                                                                                         padx=5)

        ttk.Label(frame_vcf, text="Original Circular Ref:").grid(row=1, column=0, sticky='e', padx=5, pady=3)
        ttk.Entry(frame_vcf, textvariable=self.vcf_orig_ref).grid(row=1, column=1, sticky='ew', padx=5)
        ttk.Button(frame_vcf, text="Browse",
                   command=lambda: self.browse_file(self.vcf_orig_ref, [("FASTA", "*.fasta *.fa")])).grid(row=1,
                                                                                                          column=2,
                                                                                                          padx=5)

        ttk.Label(frame_vcf, text="VCF Folder:").grid(row=2, column=0, sticky='e', padx=5, pady=3)
        ttk.Entry(frame_vcf, textvariable=self.vcf_target_dir).grid(row=2, column=1, sticky='ew', padx=5)
        ttk.Button(frame_vcf, text="Browse", command=lambda: self.browse_folder(self.vcf_target_dir)).grid(row=2,
                                                                                                           column=2,
                                                                                                           padx=5)

        ttk.Button(frame_vcf, text="Run Step 2", command=lambda: self.start_thread(self.vcf_correct_wrapper)).grid(
            row=3, column=1, pady=10)

        # 3. BAM Header Cleaning
        frame_bam = ttk.LabelFrame(self.tab_post, text="3. BAM Header Cleaning")
        frame_bam.pack(fill='x', padx=10, pady=5)

        self.bam_target_dir = tk.StringVar()
        ttk.Label(frame_bam, text="BAM Folder:").pack(side='left', padx=5)
        ttk.Entry(frame_bam, textvariable=self.bam_target_dir).pack(side='left', fill='x', expand=True, padx=5)
        ttk.Button(frame_bam, text="Browse", command=lambda: self.browse_folder(self.bam_target_dir)).pack(side='left',
                                                                                                           padx=5)
        ttk.Button(frame_bam, text="Run Step 3", command=lambda: self.start_thread(self.bam_clean_wrapper)).pack(
            side='left', padx=5)

        # 4. Final Collection
        frame_final = ttk.LabelFrame(self.tab_post, text="4. Final Collection")
        frame_final.pack(fill='x', padx=10, pady=5)

        self.final_main_dir = tk.StringVar()
        ttk.Label(frame_final, text="Output Parent Folder:").pack(side='left', padx=5)
        ttk.Entry(frame_final, textvariable=self.final_main_dir).pack(side='left', fill='x', expand=True, padx=5)
        ttk.Button(frame_final, text="Browse", command=lambda: self.browse_folder(self.final_main_dir)).pack(
            side='left', padx=5)
        ttk.Button(frame_final, text="Run Step 4", command=lambda: self.start_thread(self.collect_final_wrapper)).pack(
            side='left', padx=5)

        # Clear All Button
        ttk.Button(self.tab_post, text="Clear All Fields", command=self.clear_postproc_fields).pack(pady=10)

    def clear_postproc_fields(self):
        self.raw_main_dir.set("")
        self.vcf_target_dir.set("")
        self.bam_target_dir.set("")
        self.final_main_dir.set("")
        self.log("Post-processing directory fields cleared.")

    # =========================================================================
    # PIPELINE EXECUTION ENGINE
    # =========================================================================
    def stop_pipeline(self):
        self.stop_requested = True
        self.log("\n>>> STOP REQUESTED. Terminating processes... <<<")
        self.update_progress_ui(0, 0, custom_text="Stopping...")

        if self.current_process:
            try:
                # Kills the process group to ensure child processes (samtools, minimap2, etc) are also killed.
                os.killpg(os.getpgid(self.current_process.pid), signal.SIGTERM)
            except Exception as e:
                self.log(f"Error terminating process: {e}")

    def _handle_stop(self):
        self.update_progress_ui(1, 1, custom_text="Run Stopped.")
        self.safe_show_info("Stopped", "The process was stopped by the user.")

    def run_pipeline_wrapper(self):
        script_defaults = self.parse_script_defaults(self.script_path.get())
        overrides = []
        for v, s in self.var_states.items():
            if s['override_var'].get() and script_defaults.get(v, "Not Found").replace(" (Hardcoded)",
                                                                                       "") != "Not Found":
                overrides.append(v)

        if overrides:
            msg = "Manual Overrides Active:\n"
            for o in overrides:
                msg += f" • {o} = {self.var_states[o]['entry_var'].get()}\n"
            msg += "\nProceed with these custom settings?"
            if not messagebox.askyesno("Confirm Settings", msg):
                return

        success, wd = self._execute_pipeline()

        if self.stop_requested:
            return self._handle_stop()

        if success:
            self.update_progress_ui(1, 1, custom_text="Pipeline Execution Complete!")
            self.safe_show_info("Success", f"Pipeline Complete.\nFolder: {wd}")
            self.raw_main_dir.set(wd)
            self.final_main_dir.set(wd)
        elif wd and os.path.exists(wd):
            self.update_progress_ui(1, 1, custom_text="Pipeline Partially Failed.")
            self.safe_show_info("Partial Failure", f"Pipeline failed, but output folder exists:\n{wd}")
        else:
            self.update_progress_ui(1, 1, custom_text="Pipeline Failed.")
            self.safe_show_error("Error", "Pipeline Failed and no output folder detected.")

    def _execute_pipeline(self):
        script = self.script_path.get()
        inp = self.input_folder.get()
        out_base = self.output_folder.get().strip()

        if not os.path.exists(script) or not os.path.exists(inp):
            self.log("Error: Script or Input Folder missing.")
            return False, None

        total_files = self.count_input_files(inp)
        files_done_count = 0
        self.log(f"Detected {total_files} input samples in {inp}")
        self.update_progress_ui(0, total_files, "Starting")

        env = os.environ.copy()
        for var, state in self.var_states.items():
            if state['override_var'].get():
                env[var] = state['entry_var'].get()

        env['ref'] = self.ref_path.get()
        env['regions_bed'] = self.bed_path.get()
        if self.adapter_path.get():
            env['ADAPTER_FILE'] = self.adapter_path.get()

        exec_dir = out_base if out_base else os.path.dirname(inp)
        captured_dir = None

        try:
            self.log(f"--- Launching bash {os.path.basename(script)} ---")
            # preexec_fn=os.setsid creates a process group, critical for killing child tasks later
            self.current_process = subprocess.Popen(["bash", script, inp], cwd=exec_dir, env=env,
                                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                                    universal_newlines=True,
                                                    preexec_fn=os.setsid)

            for line in self.current_process.stdout:
                line_clean = line.strip()
                self.log(line_clean)

                if line_clean.startswith("Processing "):
                    sample_name = line_clean.replace("Processing sample: ", "").replace("Processing ", "").strip()
                    display_count = min(files_done_count + 1, total_files)
                    self.update_progress_ui(files_done_count, total_files, custom_text=f"Working on {display_count}/{total_files}: {sample_name}...")

                if line_clean.startswith("Completed: ") or line_clean.startswith("Completed ") or line_clean.startswith("Completed Diagnostic Split: "):
                    files_done_count += 1
                    display_count = min(files_done_count, total_files)
                    self.update_progress_ui(display_count, total_files, "Processing")

                if "All outputs written to:" in line_clean:
                    captured_dir = os.path.join(exec_dir, line_clean.split("written to:")[1].strip())

            self.current_process.wait()

            if self.stop_requested:
                self.log("--- PIPELINE ABORTED BY USER ---")
                return False, None

            self.update_progress_ui(total_files, total_files, "Finished")
            return (self.current_process.returncode == 0), captured_dir

        except Exception as e:
            if not self.stop_requested:
                self.log(f"Execution Error: {e}")
                self.update_progress_ui(0, 0, custom_text="Execution Error!")
            return False, None
        finally:
            self.current_process = None

    # =========================================================================
    # CORE POST-PROCESSING LOGIC
    # =========================================================================
    def _run_raw_collection(self, root):
        if not root or not os.path.exists(root): return False
        vcf_dir = os.path.join(root, "vcfs")
        bam_dir = os.path.join(root, "sorted_bams")
        os.makedirs(vcf_dir, exist_ok=True)
        os.makedirs(bam_dir, exist_ok=True)

        self.log(f"--- ORGANIZING FILES INSIDE: {os.path.basename(root)} ---")
        count = 0

        for r, d, f in os.walk(root):
            if self.stop_requested: return False
            if os.path.abspath(r) == os.path.abspath(vcf_dir) or \
                    os.path.abspath(r) == os.path.abspath(bam_dir) or \
                    "Final_Pipeline_Results" in r:
                continue

            for file in f:
                src = os.path.join(r, file)
                if file.lower().endswith(('.vcf', '.vcf.gz')):
                    self.copy_with_suffix(src, vcf_dir)
                    count += 1
                elif file.endswith(('sorted.bam', 'sorted.bam.bai', 'sorted.bai')):
                    self.copy_with_suffix(src, bam_dir)
                    count += 1

        self.log(f"Step 1: Raw files organized. Copied {count} files.")
        return True

    def _run_vcf_correction(self, folder, orig_ref):
        if not folder or not os.path.exists(orig_ref): return False
        self.log("--- STARTING VCF CORRECTION ---")

        try:
            pos_map = self.create_position_map(orig_ref, split_pos=8284)
            files = [f for f in os.listdir(folder) if f.lower().endswith('.vcf') and "corrected" not in f]
            count = 0

            for f in files:
                if self.stop_requested: return False
                input_path = os.path.join(folder, f)
                output_path = os.path.join(folder, f.replace(".vcf", "_corrected.vcf"))

                with open(input_path, "r") as fin, open(output_path, "w") as fout:
                    for line in fin:
                        if line.startswith("#"):
                            fout.write(line)
                        else:
                            parts = line.strip().split("\t")
                            if len(parts) > 1 and parts[1].isdigit():
                                lin_pos = int(parts[1])
                                orig_pos = pos_map.get(lin_pos, lin_pos)
                                parts[1] = str(orig_pos)
                                fout.write("\t".join(parts) + "\n")
                            else:
                                fout.write(line)
                count += 1

            self.log(f"Step 2: VCF positions corrected for {count} files.")
            return True
        except Exception as e:
            self.log(f"Error during VCF correction: {e}")
            return False

    def _run_bam_cleaning(self, folder):
        if not folder or not os.path.exists(folder): return False
        self.log("--- STARTING BAM CLEANING ---")

        bam_files = [f for f in os.listdir(folder) if f.endswith(".bam") and "cleaned" not in f]
        count = 0

        for f in bam_files:
            if self.stop_requested: return False
            input_path = os.path.join(folder, f)
            output_path = os.path.join(folder, f.replace(".bam", "_cleaned.bam"))

            try:
                header_sam = os.path.join(folder, "temp_header.sam")
                with open(header_sam, "w") as h:
                    subprocess.run(["samtools", "view", "-H", input_path], stdout=h, check=True)

                new_header = os.path.join(folder, "new_header.sam")
                with open(header_sam, "r") as fin, open(new_header, "w") as fout:
                    subprocess.run(["sed", "-E", "/^@PG.*ID:samtools/d"], stdin=fin,
                                   stdout=fout, check=True)

                with open(output_path, "wb") as out_bam:
                    subprocess.run(["samtools", "reheader", new_header, input_path], stdout=out_bam, check=True)

                subprocess.run(["samtools", "index", output_path], check=True)

                if os.path.exists(header_sam): os.remove(header_sam)
                if os.path.exists(new_header): os.remove(new_header)
                count += 1
            except Exception as e:
                self.log(f"Error processing {f}: {e}")

        self.log(f"Step 3: BAM headers cleaned for {count} files.")
        return True

    def _run_final_collection(self, root):
        if not root or not os.path.exists(root): return False

        final_dir = os.path.join(root, "Final_Pipeline_Results")
        final_vcf_dir = os.path.join(final_dir, "corrected_vcfs")
        final_bam_dir = os.path.join(final_dir, "cleaned_bams")

        os.makedirs(final_vcf_dir, exist_ok=True)
        os.makedirs(final_bam_dir, exist_ok=True)

        self.log(f"--- FINAL COLLECTION -> {os.path.basename(final_dir)} ---")

        src_vcf = os.path.join(root, "vcfs")
        src_bam = os.path.join(root, "sorted_bams")
        count = 0

        if os.path.exists(src_vcf):
            for f in os.listdir(src_vcf):
                if self.stop_requested: return False
                if "corrected" in f and f.endswith(('.vcf', '.vcf.gz')):
                    self.copy_with_suffix(os.path.join(src_vcf, f), final_vcf_dir)
                    count += 1

        if os.path.exists(src_bam):
            for f in os.listdir(src_bam):
                if self.stop_requested: return False
                if "cleaned" in f and f.endswith(('.bam', '.bai', '.bam.bai')):
                    self.copy_with_suffix(os.path.join(src_bam, f), final_bam_dir)
                    count += 1

        self.log(f"Step 4: Final results gathered. Copied {count} files.")
        return True

    # =========================================================================
    # UI WRAPPERS FOR POST-PROCESSING BUTTONS
    # =========================================================================
    def collect_raw_wrapper(self):
        self.update_progress_ui(0, 0, custom_text="Organizing Raw Files...")
        if self._run_raw_collection(self.raw_main_dir.get()):
            self.update_progress_ui(1, 1, custom_text="Raw Files Organized!")
            self.safe_show_info("Done", "Raw files collected.")
        elif self.stop_requested:
            self._handle_stop()

    def vcf_correct_wrapper(self):
        self.update_progress_ui(0, 0, custom_text="Correcting VCFs...")
        if self._run_vcf_correction(self.vcf_target_dir.get(), self.vcf_orig_ref.get()):
            self.update_progress_ui(1, 1, custom_text="VCF Correction Complete!")
            self.safe_show_info("Done", "VCF correction complete.")
        elif self.stop_requested:
            self._handle_stop()

    def bam_clean_wrapper(self):
        self.update_progress_ui(0, 0, custom_text="Cleaning BAM Headers...")
        if self._run_bam_cleaning(self.bam_target_dir.get()):
            self.update_progress_ui(1, 1, custom_text="BAM Headers Cleaned!")
            self.safe_show_info("Done", "BAM headers cleaned.")
        elif self.stop_requested:
            self._handle_stop()

    def collect_final_wrapper(self):
        self.update_progress_ui(0, 0, custom_text="Gathering Final Results...")
        if self._run_final_collection(self.final_main_dir.get()):
            self.update_progress_ui(1, 1, custom_text="Final Collection Complete!")
            self.safe_show_info("Done", "Final collection complete. Files moved to Final_Pipeline_Results.")
        elif self.stop_requested:
            self._handle_stop()

    # =========================================================================
    # FULL AUTO CYCLE
    # =========================================================================
    def run_automated_cycle(self):
        self.log("=== AUTO-CYCLE START ===")
        success, folder = self._execute_pipeline()

        if self.stop_requested:
            return self._handle_stop()

        if success and folder:
            self.update_progress_ui(0, 4, custom_text="Post-Processing: Organizing Raw Files (Step 1/4)")
            self.raw_main_dir.set(folder)
            self._run_raw_collection(folder)
            if self.stop_requested: return self._handle_stop()

            self.update_progress_ui(1, 4, custom_text="Post-Processing: Correcting VCFs (Step 2/4)")
            vcf_dir = os.path.join(folder, "vcfs")
            self.vcf_target_dir.set(vcf_dir)
            self._run_vcf_correction(vcf_dir, self.vcf_orig_ref.get())
            if self.stop_requested: return self._handle_stop()

            self.update_progress_ui(2, 4, custom_text="Post-Processing: Cleaning BAMs (Step 3/4)")
            bam_dir = os.path.join(folder, "sorted_bams")
            self.bam_target_dir.set(bam_dir)
            self._run_bam_cleaning(bam_dir)
            if self.stop_requested: return self._handle_stop()

            self.update_progress_ui(3, 4, custom_text="Post-Processing: Gathering Results (Step 4/4)")
            self.final_main_dir.set(folder)
            self._run_final_collection(folder)
            if self.stop_requested: return self._handle_stop()

            self.update_progress_ui(4, 4, custom_text="All Processing Complete!")
            self.safe_show_info("Auto Cycle", "Full Pipeline + Post-Processing Finished Successfully!")
        elif folder:
            self.update_progress_ui(1, 1, custom_text="Auto Cycle Partially Failed.")
            self.safe_show_info("Auto Cycle Failed",
                                "Pipeline crashed, but an output folder was generated. Post-processing aborted.")
        else:
            self.update_progress_ui(1, 1, custom_text="Auto Cycle Failed.")

    # =========================================================================
    # HELPERS
    # =========================================================================
    def copy_with_suffix(self, src, dest_dir):
        base = os.path.basename(src)
        name, ext = os.path.splitext(base)
        candidate = os.path.join(dest_dir, base)
        if os.path.exists(candidate):
            i = 1
            while True:
                new_name = f"{name}_{i}{ext}"
                candidate = os.path.join(dest_dir, new_name)
                if not os.path.exists(candidate):
                    break
                i += 1
        shutil.copy2(src, candidate)
        self.log(f"Copied: {base}")

    def create_position_map(self, fasta, split_pos=8284):
        rec = list(SeqIO.parse(fasta, "fasta"))[0]
        g_len = len(rec.seq)
        mapping = {}
        for i in range(g_len):
            orig = i + 1
            lin = (i + (g_len - split_pos) + 1) if i < split_pos else (i - split_pos + 1)
            mapping[lin] = orig
        return mapping

    def browse_file(self, var, filetypes=[("All Files", "*.*")]):
        p = filedialog.askopenfilename(filetypes=filetypes)
        if p: var.set(p)

    def browse_folder(self, var):
        p = filedialog.askdirectory()
        if p: var.set(p)

    def start_thread(self, func):
        def wrapper():
            self.set_ui_state(is_running=True)
            self.stop_requested = False
            try:
                func()
            finally:
                self.set_ui_state(is_running=False)

        threading.Thread(target=wrapper, daemon=True).start()


if __name__ == "__main__":
    root = tk.Tk()
    app = MitoPipelineDashboard(root)
    root.mainloop()