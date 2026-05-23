#!/usr/bin/env python3
import os
import gzip
import re
import datetime
import threading
import queue

import numpy as np
import pandas as pd
# Matplotlib imports deferred to plotting functions to prevent GUI hang on initial startup.

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

# =============================================================================
#  CONSTANTS & CONFIG
# =============================================================================

# Standard mitochondrial regions (rCRS coordinates)
MITO_REGIONS = [
    {"name": "HV1", "start": 16024, "end": 16383, "color": "#ff9999"},
    {"name": "HV2", "start": 57, "end": 372, "color": "#66b3ff"},
    {"name": "HV3", "start": 438, "end": 574, "color": "#99ff99"},
    {"name": "Coding", "start": 575, "end": 16023, "color": "#e0e0e0"}
]

# Homopolymer Regions (Inclusive Check)
HOMOPOLYMER_REGIONS = [
    {"start": 285, "end": 291, "name": "HP_A_6"},
    {"start": 302, "end": 309, "name": "HP_C_7"},
    {"start": 310, "end": 315, "name": "HP_C_5"},
    {"start": 451, "end": 455, "name": "HP_T_4"},
    {"start": 455, "end": 459, "name": "HP_C_4"},
    {"start": 493, "end": 498, "name": "HP_C_5"},
    {"start": 533, "end": 537, "name": "HP_C_4"},
    {"start": 540, "end": 544, "name": "HP_C_4"},
    {"start": 556, "end": 560, "name": "HP_C_4"},
    {"start": 567, "end": 573, "name": "HP_C_6"},
    {"start": 16032, "end": 16036, "name": "HP_G_4"},
    {"start": 16161, "end": 16166, "name": "HP_A_5"},
    {"start": 16179, "end": 16183, "name": "HP_A_4"},
    {"start": 16183, "end": 16188, "name": "HP_C_5"},
    {"start": 16189, "end": 16193, "name": "HP_C_4"},
    {"start": 16258, "end": 16262, "name": "HP_C_4"},
    {"start": 16362, "end": 16366, "name": "HP_C_4"},
]


# =============================================================================
#  BACKEND LOGIC
# =============================================================================

def open_text_maybe_gz(path: str):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def infer_type(ref: str, alt: str) -> str:
    if alt is None or alt.startswith("<") or alt == "*" or alt == ".":
        return "OTHER"
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    return "INDEL"


def parse_info_field(info: str) -> dict:
    d = {}
    if not info or info == ".":
        return d
    parts = info.split(";")
    for p in parts:
        if "=" in p:
            k, v = p.split("=", 1)
            d[k] = v
        else:
            d[p] = True
    return d


def safe_float(x):
    try:
        return float(x)
    except Exception:
        return np.nan


def extract_dp(info_dict):
    """Try to find depth in common INFO keys."""
    for k in ["DP", "TotalDepth", "DEPTH", "dp"]:
        if k in info_dict:
            return safe_float(info_dict[k])
    return np.nan


def check_homopolymer(pos: int) -> str:
    for hp in HOMOPOLYMER_REGIONS:
        if hp["start"] <= pos <= hp["end"]:
            return f"Yes ({hp['name']})"
    return "No"


def parse_vcf_variants(vcf_path: str):
    variants = []
    try:
        with open_text_maybe_gz(vcf_path) as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                line = line.rstrip("\n")
                cols = line.split("\t")
                if len(cols) < 8:
                    continue
                chrom, pos, vid, ref, alts, qual, flt, info = cols[:8]

                try:
                    pos_i = int(pos)
                except Exception:
                    continue

                info_dict = parse_info_field(info)
                qual_f = safe_float(qual)
                dp_val = extract_dp(info_dict)

                hp_status = check_homopolymer(pos_i)

                for alt in alts.split(","):
                    alt = alt.strip()
                    vtype = infer_type(ref, alt)
                    variants.append({
                        "chrom": chrom,
                        "pos": pos_i,
                        "ref": ref,
                        "alt": alt,
                        "qual": qual_f,
                        "dp": dp_val,
                        "type": vtype,
                        "info": info_dict,
                        "in_hp": hp_status
                    })
    except Exception as e:
        print(f"Error parsing {vcf_path}: {e}")
        return []
    return variants


def extract_af(info_dict: dict):
    for k in ("AF", "VAF", "AlleleFreq", "ALLELE_FREQ", "allele_freq"):
        if k in info_dict:
            v = str(info_dict[k])
            parts = re.split(r"[,\|]", v)
            vals = [safe_float(p) for p in parts]
            vals = [x for x in vals if np.isfinite(x)]
            return max(vals) if vals else np.nan
    return np.nan


def build_matrix(vcf_paths):
    all_records = []
    file_labels = []

    for p in vcf_paths:
        label = os.path.basename(p)
        file_labels.append(label)
        var_list = parse_vcf_variants(p)
        for v in var_list:
            all_records.append({
                "chrom": v["chrom"],
                "pos": v["pos"],
                "ref": v["ref"],
                "alt": v["alt"],
                "type": v["type"],
                "qual": v["qual"],
                "dp": v["dp"],
                "af": extract_af(v["info"]),
                "in_hp": v["in_hp"],
                "source_file": label,
            })

    if not all_records:
        raise ValueError("No variant records parsed.")

    raw_df = pd.DataFrame(all_records)

    # 1. Presence Matrix
    presence_matrix = (
        raw_df.assign(present=1)
        .pivot_table(index=["chrom", "pos"], columns="source_file", values="present", aggfunc="max", fill_value=0)
    ).reset_index()

    # 2. Depth Matrix (for plotting numbers)
    depth_matrix = (
        raw_df.pivot_table(index=["chrom", "pos"], columns="source_file", values="dp", aggfunc="max")
    )

    # 3. Quality Matrix (for the third plot)
    quality_matrix = (
        raw_df.pivot_table(index=["chrom", "pos"], columns="source_file", values="qual", aggfunc="max")
    )

    # 4. Averages per File (for X-Axis labels)
    avg_depths = raw_df.groupby("source_file")["dp"].mean().to_dict()
    avg_quals = raw_df.groupby("source_file")["qual"].mean().to_dict()

    def agg_meta(x):
        u = sorted(list(set(x)))
        return ",".join(u) if len(u) > 1 else u[0]

    meta_df = raw_df.groupby(["chrom", "pos"]).agg({
        "ref": "first",
        "alt": agg_meta,
        "type": "first",
        "in_hp": "first"
    }).reset_index()

    collapsed_df = pd.merge(meta_df, presence_matrix, on=["chrom", "pos"])

    file_cols_actual = [c for c in collapsed_df.columns if c in file_labels]
    row_sums = collapsed_df[file_cols_actual].sum(axis=1)
    n_files = len(file_labels)

    status_list = []
    for count in row_sums:
        if count == n_files:
            status_list.append("Consensus (All)")
        elif count == 1:
            status_list.append("Unique (Singleton)")
        else:
            status_list.append("Discordant (Some)")

    collapsed_df["status"] = status_list
    collapsed_df = collapsed_df.sort_values(by=["chrom", "pos"], ascending=[True, True])

    per_file = []
    for label in file_labels:
        if label in collapsed_df.columns:
            sub = collapsed_df[collapsed_df[label] == 1]
            per_file.append({
                "file": label,
                "collapsed_sites": len(sub),
                "unique_sites": len(sub[sub["status"] == "Unique (Singleton)"]),
                "consensus_sites": len(sub[sub["status"] == "Consensus (All)"]),
                "avg_variant_depth": avg_depths.get(label, np.nan),
                "avg_variant_qual": avg_quals.get(label, np.nan)
            })
    per_file_df = pd.DataFrame(per_file).sort_values("file")

    return raw_df, collapsed_df, per_file_df, file_labels, depth_matrix, avg_depths, quality_matrix, avg_quals


# =============================================================================
#  GUI CLASS
# =============================================================================

class VCFCompareGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("VCF Analyzer - Scrollable Table")
        self.geometry("1400x850")

        self.vcf_paths = []
        self.output_dir = tk.StringVar()

        self.variants_df = None
        self.raw_df = None
        self.file_labels = []

        self.sort_col = None
        self.sort_desc = False

        self.msg_queue = queue.Queue()
        self._build_ui()
        self.check_queue()

    def _build_ui(self):
        # Top Config
        frame_top = ttk.LabelFrame(self, text="Configuration")
        frame_top.pack(side=tk.TOP, fill=tk.X, padx=10, pady=5)

        btn_frame = ttk.Frame(frame_top)
        btn_frame.pack(fill=tk.X, padx=5, pady=5)
        ttk.Button(btn_frame, text="+ Add Files", command=self.add_files).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="+ Add Folder", command=self.add_folder).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text="Clear List", command=self.clear_list).pack(side=tk.LEFT, padx=5)

        out_frame = ttk.Frame(frame_top)
        out_frame.pack(fill=tk.X, padx=5, pady=5)
        ttk.Label(out_frame, text="Output Directory:").pack(side=tk.LEFT, padx=5)
        ttk.Entry(out_frame, textvariable=self.output_dir, width=60).pack(side=tk.LEFT, padx=5)
        ttk.Button(out_frame, text="Browse", command=self.browse_output).pack(side=tk.LEFT, padx=5)

        # Middle
        mid_frame = ttk.Frame(self)
        mid_frame.pack(side=tk.TOP, fill=tk.BOTH, padx=10, pady=5, expand=False)

        list_frame = ttk.LabelFrame(mid_frame, text="Selected Files")
        list_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 5))

        self.listbox = tk.Listbox(list_frame, height=8, selectmode=tk.EXTENDED)
        self.listbox.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        ttk.Button(list_frame, text="Remove Selected", command=self.remove_selected_files).pack(anchor="e", padx=5,
                                                                                                pady=2)

        action_frame = ttk.Frame(mid_frame)
        action_frame.pack(side=tk.RIGHT, fill=tk.Y, padx=(5, 0))

        self.btn_run = ttk.Button(action_frame, text="RUN ANALYSIS", command=self.run_thread, state="normal")
        self.btn_run.pack(fill=tk.X, pady=10, ipady=10)

        self.progress = ttk.Progressbar(action_frame, mode='indeterminate', length=200)

        # Bottom Table
        frame_main = ttk.Frame(self)
        frame_main.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        self.tree = ttk.Treeview(frame_main, show="headings")
        ys = ttk.Scrollbar(frame_main, orient="vertical", command=self.tree.yview)
        xs = ttk.Scrollbar(frame_main, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=ys.set, xscrollcommand=xs.set)

        self.tree.grid(row=0, column=0, sticky='nsew')
        ys.grid(row=0, column=1, sticky='ns')
        xs.grid(row=1, column=0, sticky='ew')

        frame_main.grid_rowconfigure(0, weight=1)
        frame_main.grid_columnconfigure(0, weight=1)

        self.tree.tag_configure("Consensus (All)", background="#c2f0c2")
        self.tree.tag_configure("Discordant (Some)", background="#fff4cc")
        self.tree.tag_configure("Unique (Singleton)", background="#ffcccc")

        legend = ttk.Frame(self)
        legend.pack(fill=tk.X, padx=10, pady=5)
        tk.Label(legend, text=" ■ Consensus ", bg="#c2f0c2").pack(side=tk.LEFT, padx=2)
        tk.Label(legend, text=" ■ Discordant ", bg="#fff4cc").pack(side=tk.LEFT, padx=2)
        tk.Label(legend, text=" ■ Unique ", bg="#ffcccc").pack(side=tk.LEFT, padx=2)
        ttk.Label(legend, text="(Click headers to sort)").pack(side=tk.RIGHT, padx=10)

        self.status = tk.StringVar(value="Ready.")
        ttk.Label(self, textvariable=self.status, relief=tk.SUNKEN, anchor="w").pack(side=tk.BOTTOM, fill=tk.X)

    # =========================================================================
    #  LOGIC
    # =========================================================================

    def browse_output(self):
        d = filedialog.askdirectory()
        if d: self.output_dir.set(d)

    def add_files(self):
        paths = filedialog.askopenfilenames(parent=self, title="Select VCFs")
        if not paths: return

        if not self.output_dir.get():
            self.output_dir.set(os.path.dirname(paths[0]))

        for p in paths:
            if p not in self.vcf_paths:
                self.vcf_paths.append(p)
                self.listbox.insert(tk.END, os.path.basename(p))
        self.status.set(f"Total files: {len(self.vcf_paths)}")

    def add_folder(self):
        folder = filedialog.askdirectory(title="Select Folder")
        if not folder: return

        if not self.output_dir.get():
            self.output_dir.set(os.path.dirname(folder))

        found = []
        for root, _, files in os.walk(folder):
            for fn in files:
                if fn.lower().endswith(('.vcf', '.vcf.gz')):
                    found.append(os.path.join(root, fn))

        for p in sorted(found):
            if p not in self.vcf_paths:
                self.vcf_paths.append(p)
                self.listbox.insert(tk.END, os.path.basename(p))
        self.status.set(f"Total files: {len(self.vcf_paths)}")

    def clear_list(self):
        self.vcf_paths = []
        self.listbox.delete(0, tk.END)

    def remove_selected_files(self):
        selection = self.listbox.curselection()
        if not selection:
            return
        indices = sorted(list(selection), reverse=True)
        for i in indices:
            self.listbox.delete(i)
            if i < len(self.vcf_paths):
                self.vcf_paths.pop(i)
        self.status.set(f"Removed items. Remaining: {len(self.vcf_paths)}")

    def run_thread(self):
        if len(self.vcf_paths) < 2:
            messagebox.showerror("Error", "Select at least 2 VCF files.")
            return

        if not self.output_dir.get():
            messagebox.showerror("Error", "Select an Output Directory.")
            return

        self.btn_run.config(state="disabled")
        self.progress.pack(side=tk.BOTTOM, fill=tk.X, pady=5)
        self.progress.start(10)
        self.status.set("Processing...")

        t = threading.Thread(target=self._process_data)
        t.start()

    def _process_data(self):
        try:
            # Unpack the new return values (quality_matrix, avg_quals)
            (raw_df, collapsed_df, perf_df, labels,
             depth_df, avg_depths,
             qual_df, avg_quals) = build_matrix(self.vcf_paths)

            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            res_dir = os.path.join(self.output_dir.get(), f"VCF_Results_{timestamp}")
            os.makedirs(res_dir, exist_ok=True)

            summary_csv_df = collapsed_df.copy()
            for col in labels:
                summary_csv_df[col] = summary_csv_df[col].replace({1: "Present", 0: "Absent"})
            summary_csv_df.to_csv(os.path.join(res_dir, "variants_summary_matrix.csv"), index=False)
            raw_df.to_csv(os.path.join(res_dir, "variants_detailed_raw.csv"), index=False)
            perf_df.to_csv(os.path.join(res_dir, "file_statistics.csv"), index=False)

            # 1. Standard Presence Plot
            self._generate_plots(res_dir, collapsed_df, labels, depth_df, avg_depths)

            # 2. Problematic/Detail Plot
            self._generate_problematic_plot(res_dir, collapsed_df, labels, depth_df, avg_depths)

            # 3. New Quality Heatmap Plot
            self._generate_quality_plot(res_dir, collapsed_df, labels, qual_df, avg_quals)

            self.msg_queue.put(("DONE", (raw_df, collapsed_df, labels, res_dir)))

        except Exception as e:
            self.msg_queue.put(("ERROR", str(e)))

    def _generate_plots(self, folder, df, file_cols, depth_df, avg_depths):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        plot_df = df.sort_values(by=["chrom", "pos"])
        if not depth_df.index.names == ["chrom", "pos"]:
            depth_lookup = depth_df
        else:
            depth_lookup = depth_df

        matrix = plot_df[file_cols].values
        status_list = plot_df["status"].values
        pos_list = plot_df["pos"].values
        chrom_list = plot_df["chrom"].values

        rows, cols = matrix.shape
        rgb_grid = np.ones((rows, cols, 3))

        c_green = np.array([0.76, 0.94, 0.76])
        c_orange = np.array([1.0, 0.96, 0.8])
        c_red = np.array([1.0, 0.8, 0.8])

        for r in range(rows):
            stat = status_list[r]
            if "Consensus" in stat:
                color = c_green
            elif "Unique" in stat:
                color = c_red
            else:
                color = c_orange

            for c in range(cols):
                if matrix[r, c] == 1:
                    rgb_grid[r, c] = color

        fig_h = max(12, rows * 0.15)
        if fig_h > 120: fig_h = 120

        fig = plt.figure(figsize=(15, fig_h), constrained_layout=True)
        gs = fig.add_gridspec(1, 2, width_ratios=[0.3, 10], wspace=0.01)

        ax_regions = fig.add_subplot(gs[0])
        ax_matrix = fig.add_subplot(gs[1])

        ax_matrix.imshow(rgb_grid, aspect='auto', interpolation='nearest', origin='upper')
        ax_matrix.set_title("Variant Presence by Position (Depth Overlay)")

        x_labels_new = []
        for fc in file_cols:
            ad = avg_depths.get(fc, 0)
            if np.isnan(ad): ad = 0
            x_labels_new.append(f"{fc}\n(Avg DP: {int(ad)})")

        ax_matrix.set_xlabel("VCF Files")
        ax_matrix.set_xticks(range(len(file_cols)))
        ax_matrix.set_xticklabels(x_labels_new, rotation=90, fontsize=9)

        for r in range(rows):
            p_chrom = chrom_list[r]
            p_pos = pos_list[r]

            for c, fcol in enumerate(file_cols):
                if matrix[r, c] == 1:
                    try:
                        val = depth_lookup.loc[(p_chrom, p_pos), fcol]
                        if pd.notnull(val):
                            txt = str(int(val))
                            ax_matrix.text(c, r, txt, ha='center', va='center',
                                           fontsize=7, color='black')
                    except KeyError:
                        pass

        step = 1 if rows < 200 else int(rows / 100)
        y_ticks = range(0, rows, step)
        y_labels = [str(pos_list[i]) for i in y_ticks]

        ax_matrix.set_yticks(y_ticks)
        ax_matrix.set_yticklabels(y_labels, fontsize=8)

        patches = [
            mpatches.Patch(color=c_green, label='Consensus'),
            mpatches.Patch(color=c_orange, label='Discordant'),
            mpatches.Patch(color=c_red, label='Unique'),
            mpatches.Patch(color=[1, 1, 1], label='Absent', edgecolor='gray')
        ]
        ax_matrix.legend(handles=patches, loc='upper right', bbox_to_anchor=(1.15, 1))

        ax_regions.set_ylim(rows - 0.5, -0.5)
        ax_regions.set_xlim(0, 1)
        ax_regions.axis('off')

        for i, pos in enumerate(pos_list):
            reg_color = "#ffffff"
            for mr in MITO_REGIONS:
                if mr["start"] <= pos <= mr["end"]:
                    reg_color = mr["color"]
                    break
            rect = mpatches.Rectangle((0, i - 0.5), 1, 1, color=reg_color, ec=None)
            ax_regions.add_patch(rect)

        current_region = None
        block_start = 0

        def label_block(r_name, start_idx, end_idx):
            if not r_name: return
            center = (start_idx + end_idx) / 2
            ax_regions.text(0.5, center, r_name, ha='center', va='center',
                            rotation=90, fontsize=8, fontweight='bold', color='black')

        for i, pos in enumerate(pos_list):
            r_name = None
            for mr in MITO_REGIONS:
                if mr["start"] <= pos <= mr["end"]:
                    r_name = mr["name"]
                    break

            if r_name != current_region:
                label_block(current_region, block_start, i - 1)
                current_region = r_name
                block_start = i

        label_block(current_region, block_start, rows - 1)

        plt.savefig(os.path.join(folder, "plot_presence_by_position.png"), dpi=150, bbox_inches='tight')
        plt.close()

    def _generate_problematic_plot(self, folder, df, file_cols, depth_df, avg_depths):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        plot_df = df.sort_values(by=["chrom", "pos"])
        if not depth_df.index.names == ["chrom", "pos"]:
            depth_lookup = depth_df
        else:
            depth_lookup = depth_df

        matrix = plot_df[file_cols].values
        pos_list = plot_df["pos"].values
        chrom_list = plot_df["chrom"].values
        type_list = plot_df["type"].values
        hp_list = plot_df["in_hp"].values

        rows, cols = matrix.shape
        rgb_grid = np.ones((rows, cols, 3))

        c_indel = np.array([0.82, 0.13, 0.56])  # Magenta
        c_hp_snp = np.array([1.0, 0.65, 0.0])  # Orange
        c_clean = np.array([0.68, 0.85, 0.9])  # Light Blue

        y_labels_custom = []

        for r in range(rows):
            v_type = type_list[r]
            v_hp = hp_list[r]

            row_color = c_clean
            if v_type == "INDEL":
                row_color = c_indel
            elif v_hp != "No":
                row_color = c_hp_snp

            for c in range(cols):
                if matrix[r, c] == 1:
                    rgb_grid[r, c] = row_color

            label_text = str(pos_list[r])
            if v_hp != "No":
                hp_name = v_hp.replace("Yes (", "").replace(")", "")
                label_text += f" ({hp_name})"
            elif v_type == "INDEL":
                label_text += " (INDEL)"
            y_labels_custom.append(label_text)

        fig_h = max(12, rows * 0.15)
        if fig_h > 120: fig_h = 120

        fig = plt.figure(figsize=(16, fig_h), constrained_layout=True)
        gs = fig.add_gridspec(1, 3, width_ratios=[0.3, 0.3, 10], wspace=0.01)

        ax_mito = fig.add_subplot(gs[0])
        ax_hp = fig.add_subplot(gs[1])
        ax_matrix = fig.add_subplot(gs[2])

        ax_matrix.imshow(rgb_grid, aspect='auto', interpolation='nearest', origin='upper')
        ax_matrix.set_title("Problematic Variants (Indels & Homopolymers)")

        x_labels_new = []
        for fc in file_cols:
            ad = avg_depths.get(fc, 0)
            if np.isnan(ad): ad = 0
            x_labels_new.append(f"{fc}\n(Avg DP: {int(ad)})")

        ax_matrix.set_xlabel("VCF Files")
        ax_matrix.set_xticks(range(len(file_cols)))
        ax_matrix.set_xticklabels(x_labels_new, rotation=90, fontsize=9)

        for r in range(rows):
            p_chrom = chrom_list[r]
            p_pos = pos_list[r]
            v_type = type_list[r]

            txt_color = 'black'
            if v_type == "INDEL":
                txt_color = 'white'

            for c, fcol in enumerate(file_cols):
                if matrix[r, c] == 1:
                    try:
                        val = depth_lookup.loc[(p_chrom, p_pos), fcol]
                        if pd.notnull(val):
                            txt = str(int(val))
                            ax_matrix.text(c, r, txt, ha='center', va='center',
                                           fontsize=7, color=txt_color)
                    except KeyError:
                        pass

        step = 1 if rows < 200 else int(rows / 100)
        y_ticks = range(0, rows, step)
        y_labels_subset = [y_labels_custom[i] for i in y_ticks]

        ax_matrix.set_yticks(y_ticks)
        ax_matrix.set_yticklabels(y_labels_subset, fontsize=8)

        patches = [
            mpatches.Patch(color=c_indel, label='INDEL'),
            mpatches.Patch(color=c_hp_snp, label='SNP in Homopolymer'),
            mpatches.Patch(color=c_clean, label='Clean SNP'),
            mpatches.Patch(color=[1, 1, 1], label='Absent', edgecolor='gray')
        ]
        ax_matrix.legend(handles=patches, loc='upper right', bbox_to_anchor=(1.15, 1))

        # --- Draw Mito Regions (Left Track) ---
        ax_mito.set_ylim(rows - 0.5, -0.5)
        ax_mito.set_xlim(0, 1)
        ax_mito.axis('off')

        current_region = None
        block_start = 0

        def label_block(ax, r_name, start_idx, end_idx):
            if not r_name: return
            center = (start_idx + end_idx) / 2
            ax.text(0.5, center, r_name, ha='center', va='center',
                    rotation=90, fontsize=8, fontweight='bold', color='black')

        for i, pos in enumerate(pos_list):
            reg_color = "#ffffff"
            r_name = None
            for mr in MITO_REGIONS:
                if mr["start"] <= pos <= mr["end"]:
                    reg_color = mr["color"]
                    r_name = mr["name"]
                    break

            rect = mpatches.Rectangle((0, i - 0.5), 1, 1, color=reg_color, ec=None)
            ax_mito.add_patch(rect)

            if r_name != current_region:
                label_block(ax_mito, current_region, block_start, i - 1)
                current_region = r_name
                block_start = i
        label_block(ax_mito, current_region, block_start, rows - 1)

        # --- Draw HP Track (Middle Track) ---
        ax_hp.set_ylim(rows - 0.5, -0.5)
        ax_hp.set_xlim(0, 1)
        ax_hp.axis('off')

        for i, val in enumerate(hp_list):
            color = "#eeeeee"
            if val != "No":
                color = "#444444"
            rect = mpatches.Rectangle((0, i - 0.5), 1, 1, color=color, ec=None)
            ax_hp.add_patch(rect)

        ax_hp.set_title("HP", fontsize=8, rotation=90)

        plt.savefig(os.path.join(folder, "plot_problematic_variants.png"), dpi=150, bbox_inches='tight')
        plt.close()

    def _generate_quality_plot(self, folder, df, file_cols, qual_df, avg_quals):
        """
        Generates a heatmap where COLOR intensity = Quality Score.
        TEXT in cell = Quality Score.
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import matplotlib.colors as mcolors
        import matplotlib.cm as cm

        plot_df = df.sort_values(by=["chrom", "pos"])

        # Quality lookup
        if not qual_df.index.names == ["chrom", "pos"]:
            qual_lookup = qual_df
        else:
            qual_lookup = qual_df

        matrix = plot_df[file_cols].values  # 1/0 presence
        pos_list = plot_df["pos"].values
        chrom_list = plot_df["chrom"].values

        rows, cols = matrix.shape

        # Calculate global min/max quality for color scaling
        all_vals = qual_lookup.values.flatten()
        all_vals = all_vals[~np.isnan(all_vals)]
        if len(all_vals) > 0:
            q_min, q_max = np.min(all_vals), np.max(all_vals)
        else:
            q_min, q_max = 0, 100  # Fallback

        # Use a colormap (e.g., Viridis, Plasma, or YlGnBu)
        cmap = plt.get_cmap('YlGnBu')
        norm = mcolors.Normalize(vmin=q_min, vmax=q_max)

        rgb_grid = np.ones((rows, cols, 3))  # Default White

        for r in range(rows):
            p_chrom = chrom_list[r]
            p_pos = pos_list[r]

            for c, fcol in enumerate(file_cols):
                if matrix[r, c] == 1:
                    try:
                        val = qual_lookup.loc[(p_chrom, p_pos), fcol]
                        if pd.notnull(val):
                            # Get color from cmap
                            rgba = cmap(norm(val))
                            rgb_grid[r, c] = rgba[:3]  # RGB only
                    except KeyError:
                        pass

        fig_h = max(12, rows * 0.15)
        if fig_h > 120: fig_h = 120

        fig = plt.figure(figsize=(15, fig_h), constrained_layout=True)
        gs = fig.add_gridspec(1, 2, width_ratios=[0.3, 10], wspace=0.01)

        ax_regions = fig.add_subplot(gs[0])
        ax_matrix = fig.add_subplot(gs[1])

        ax_matrix.imshow(rgb_grid, aspect='auto', interpolation='nearest', origin='upper')
        ax_matrix.set_title("Variant Quality Heatmap (Color & Number = Quality)")

        # X-Labels with Avg Quality
        x_labels_new = []
        for fc in file_cols:
            aq = avg_quals.get(fc, 0)
            if np.isnan(aq): aq = 0
            x_labels_new.append(f"{fc}\n(Avg Qual: {int(aq)})")

        ax_matrix.set_xlabel("VCF Files")
        ax_matrix.set_xticks(range(len(file_cols)))
        ax_matrix.set_xticklabels(x_labels_new, rotation=90, fontsize=9)

        # Draw Quality Numbers in Cells
        for r in range(rows):
            p_chrom = chrom_list[r]
            p_pos = pos_list[r]

            for c, fcol in enumerate(file_cols):
                if matrix[r, c] == 1:
                    try:
                        val = qual_lookup.loc[(p_chrom, p_pos), fcol]
                        if pd.notnull(val):
                            txt = str(int(val))
                            # Determine text color based on background intensity
                            # Dark background -> White text, Light -> Black
                            # Simple heuristic: if norm(val) > 0.5 -> White
                            text_color = "white" if norm(val) > 0.5 else "black"

                            ax_matrix.text(c, r, txt, ha='center', va='center',
                                           fontsize=7, color=text_color)
                    except KeyError:
                        pass

        step = 1 if rows < 200 else int(rows / 100)
        y_ticks = range(0, rows, step)
        y_labels = [str(pos_list[i]) for i in y_ticks]

        ax_matrix.set_yticks(y_ticks)
        ax_matrix.set_yticklabels(y_labels, fontsize=8)

        # Add Colorbar for reference
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax_matrix, fraction=0.046, pad=0.04)
        cbar.set_label('Quality Score')

        # --- Draw Mito Regions (Left Track) ---
        ax_regions.set_ylim(rows - 0.5, -0.5)
        ax_regions.set_xlim(0, 1)
        ax_regions.axis('off')

        current_region = None
        block_start = 0

        def label_block(ax, r_name, start_idx, end_idx):
            if not r_name: return
            center = (start_idx + end_idx) / 2
            ax.text(0.5, center, r_name, ha='center', va='center',
                    rotation=90, fontsize=8, fontweight='bold', color='black')

        for i, pos in enumerate(pos_list):
            reg_color = "#ffffff"
            r_name = None
            for mr in MITO_REGIONS:
                if mr["start"] <= pos <= mr["end"]:
                    reg_color = mr["color"]
                    r_name = mr["name"]
                    break

            rect = mpatches.Rectangle((0, i - 0.5), 1, 1, color=reg_color, ec=None)
            ax_regions.add_patch(rect)

            if r_name != current_region:
                label_block(ax_regions, current_region, block_start, i - 1)
                current_region = r_name
                block_start = i
        label_block(ax_regions, current_region, block_start, rows - 1)

        plt.savefig(os.path.join(folder, "plot_quality_heatmap.png"), dpi=150, bbox_inches='tight')
        plt.close()

    def check_queue(self):
        try:
            while True:
                msg, data = self.msg_queue.get_nowait()
                self.progress.stop()
                self.progress.pack_forget()
                self.btn_run.config(state="normal")

                if msg == "DONE":
                    raw, collapsed, labels, folder = data
                    self.raw_df = raw
                    self.variants_df = collapsed
                    self.file_labels = labels
                    self.refresh_table()
                    self.status.set(f"Saved to: {folder}")
                    messagebox.showinfo("Success", f"Analysis Complete!\nOutput:\n{folder}")

                elif msg == "ERROR":
                    self.status.set("Failed.")
                    messagebox.showerror("Error", data)
        except queue.Empty:
            pass
        finally:
            self.after(100, self.check_queue)

    def sort_by_col(self, col):
        if self.variants_df is None: return
        if self.sort_col == col:
            self.sort_desc = not self.sort_desc
        else:
            self.sort_col = col
            self.sort_desc = False
        self.variants_df = self.variants_df.sort_values(by=col, ascending=not self.sort_desc)
        self.refresh_table()

    def refresh_table(self):
        if self.variants_df is None: return
        self.tree.delete(*self.tree.get_children())

        static_cols = ["status", "chrom", "pos", "ref", "alt", "type", "in_hp"]
        cols = static_cols + self.file_labels

        self.tree["columns"] = cols

        for c in cols:
            self.tree.heading(c, text=c, command=lambda _col=c: self.sort_by_col(_col))
            if c == "status":
                w = 150
            elif c == "in_hp":
                w = 100
            elif c in self.file_labels:
                w = 120
            else:
                w = 60

            self.tree.column(c, width=w, anchor="center", stretch=False)

        df_show = self.variants_df.head(4000)

        for _, row in df_show.iterrows():
            vals = []
            for c in static_cols:
                vals.append(row[c])
            for c in self.file_labels:
                val = row[c]
                vals.append("✔" if val == 1 else "")
            stat = row["status"]
            self.tree.insert("", tk.END, values=vals, tags=(stat,))


if __name__ == "__main__":
    app = VCFCompareGUI()
    app.mainloop()