#!/usr/bin/env python3
"""
Bulk Sequence File Analyzer (BAM/FASTQ) with HTML Dashboard
----------------------------------------------------------
- Analyzes all BAM/FASTQ(.gz) files in an input folder.
- Produces per-file plots + detailed reports + consolidated CSV.
- Produces an additional run-level summary plot (averages across all files).
- Optionally generates an HTML dashboard for quick visual scanning and comparison.

Notes:
- Requires: pysam, numpy, matplotlib (plus tkinter which ships with most Python distros).
- Designed to work offline (dashboard is self-contained; no external JS/CSS CDNs).
"""

import os
import glob
import gzip
import csv
import threading
from collections import defaultdict
from datetime import datetime

import numpy as np
import pysam
import matplotlib
matplotlib.use('Agg', force=True)
import matplotlib.pyplot as plt

import tkinter as tk
from tkinter import filedialog, messagebox, font
import tkinter.ttk as ttk


# -----------------------------
# Plot styling helpers
# -----------------------------
def _apply_plot_style():
    # Use built-in Matplotlib styles only (no seaborn dependency).
    # Falls back cleanly if style name does not exist.
    for style_name in ("seaborn-v0_8-whitegrid", "seaborn-whitegrid", "ggplot"):
        try:
            plt.style.use(style_name)
            break
        except Exception:
            pass
    matplotlib.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 220,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "axes.titlepad": 10,
    })


def _safe_base_filename(file_path: str) -> str:
    base = os.path.basename(file_path)
    # Strip common FASTQ extensions robustly
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".bam", ".gz"):
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    return base


def _savefig(path: str):
    plt.savefig(path, bbox_inches="tight")
    plt.close()


# -----------------------------
# HTML dashboard helpers
# -----------------------------
def _html_escape(s: str) -> str:
    return (
        str(s)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&#39;")
    )


def _write_dashboard_html(out_dir: str, run_name: str, rows: list, assets: dict):
    """
    Create a single self-contained dashboard HTML (offline), referencing local PNGs/TXTs.
    """
    index_path = os.path.join(out_dir, f"{run_name}_Dashboard.html")

    # Lightweight CSS + JS: search + sort + image preview.
    css = r"""
    :root {
      --bg: #0b0f19;
      --panel: #121a2b;
      --panel2: #0f1524;
      --text: #e8eefc;
      --muted: #a8b3cf;
      --accent: #6ea8fe;
      --good: #50c878;
      --warn: #ffcc66;
      --bad: #ff6b6b;
      --border: rgba(255,255,255,0.08);
      --shadow: 0 10px 30px rgba(0,0,0,0.35);
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      padding: 24px;
      font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Arial, "Apple Color Emoji","Segoe UI Emoji";
      background: radial-gradient(1200px 800px at 15% 10%, #1b2a55 0%, var(--bg) 55%);
      color: var(--text);
    }
    h1 { margin: 0 0 8px; font-size: 22px; }
    .sub { color: var(--muted); margin-bottom: 18px; }
    .grid {
      display: grid;
      grid-template-columns: 1fr;
      gap: 14px;
      max-width: 1400px;
      margin: 0 auto;
    }
    .card {
      background: linear-gradient(180deg, var(--panel) 0%, var(--panel2) 100%);
      border: 1px solid var(--border);
      border-radius: 14px;
      box-shadow: var(--shadow);
      overflow: hidden;
    }
    .card-header {
      padding: 14px 16px;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      border-bottom: 1px solid var(--border);
    }
    .pill {
      display: inline-flex;
      align-items: center;
      gap: 8px;
      padding: 6px 10px;
      border-radius: 999px;
      font-size: 12px;
      color: var(--muted);
      border: 1px solid var(--border);
      background: rgba(255,255,255,0.03);
    }
    .controls {
      display: flex;
      gap: 10px;
      flex-wrap: wrap;
      padding: 14px 16px;
      align-items: center;
    }
    input[type="text"]{
      padding: 10px 12px;
      border-radius: 10px;
      border: 1px solid var(--border);
      background: rgba(255,255,255,0.03);
      color: var(--text);
      min-width: 320px;
      outline: none;
    }
    input[type="text"]::placeholder{ color: rgba(232,238,252,0.45); }
    .btn {
      padding: 9px 12px;
      border-radius: 10px;
      border: 1px solid var(--border);
      background: rgba(110,168,254,0.12);
      color: var(--text);
      cursor: pointer;
      text-decoration: none;
      display: inline-flex;
      gap: 8px;
      align-items: center;
      white-space: nowrap;
    }
    .btn:hover { background: rgba(110,168,254,0.18); }
    .btn.secondary { background: rgba(255,255,255,0.04); }
    .table-wrap { overflow: auto; }
    table { width: 100%; border-collapse: collapse; }
    thead th {
      text-align: left;
      padding: 12px 12px;
      font-size: 12px;
      color: var(--muted);
      position: sticky;
      top: 0;
      background: rgba(15,21,36,0.95);
      border-bottom: 1px solid var(--border);
      cursor: pointer;
      user-select: none;
    }
    tbody td {
      padding: 12px 12px;
      border-bottom: 1px solid var(--border);
      font-size: 13px;
      vertical-align: top;
    }
    tbody tr:hover { background: rgba(255,255,255,0.03); }
    .mono { font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", monospace; }
    .tag {
      display:inline-flex;
      padding: 3px 8px;
      border-radius: 999px;
      border: 1px solid var(--border);
      font-size: 12px;
      color: var(--muted);
      background: rgba(255,255,255,0.03);
    }
    .kpi { color: var(--text); font-weight: 600; }
    .muted { color: var(--muted); }
    .small { font-size: 12px; color: var(--muted); }
    .preview {
      display:none;
      position: fixed;
      inset: 0;
      background: rgba(0,0,0,0.78);
      z-index: 50;
      padding: 30px;
    }
    .preview-inner {
      max-width: 1200px;
      margin: 0 auto;
      background: #0f1524;
      border: 1px solid var(--border);
      border-radius: 14px;
      overflow: hidden;
      box-shadow: var(--shadow);
    }
    .preview-head {
      padding: 12px 14px;
      display:flex;
      justify-content: space-between;
      align-items:center;
      border-bottom: 1px solid var(--border);
    }
    .preview img { width: 100%; height: auto; display:block; }
    .x { cursor:pointer; padding: 8px 10px; border-radius: 10px; border:1px solid var(--border); background: rgba(255,255,255,0.04); }
    """

    js = r"""
    let sortState = { col: null, asc: true };

    function getCellText(td){ return (td.getAttribute('data-sort') || td.innerText || '').trim(); }

    function sortTable(colIndex){
      const table = document.getElementById('dataTable');
      const tbody = table.tBodies[0];
      const rows = Array.from(tbody.rows);

      const asc = (sortState.col === colIndex) ? !sortState.asc : true;
      sortState = { col: colIndex, asc: asc };

      rows.sort((a,b)=>{
        const ta = getCellText(a.cells[colIndex]);
        const tb = getCellText(b.cells[colIndex]);

        const na = parseFloat(ta);
        const nb = parseFloat(tb);
        const aIsNum = !Number.isNaN(na) && ta !== '';
        const bIsNum = !Number.isNaN(nb) && tb !== '';

        let cmp = 0;
        if(aIsNum && bIsNum){ cmp = na - nb; }
        else { cmp = ta.localeCompare(tb); }

        return asc ? cmp : -cmp;
      });

      rows.forEach(r=>tbody.appendChild(r));
      updateSortIndicators(colIndex, asc);
    }

    function updateSortIndicators(colIndex, asc){
      const ths = document.querySelectorAll('thead th');
      ths.forEach((th,i)=>{
        const base = th.getAttribute('data-base');
        th.innerText = base + (i===colIndex ? (asc ? ' ▲' : ' ▼') : '');
      });
    }

    function filterTable(){
      const q = document.getElementById('searchBox').value.toLowerCase();
      const rows = document.querySelectorAll('#dataTable tbody tr');
      rows.forEach(r=>{
        const text = r.innerText.toLowerCase();
        r.style.display = text.includes(q) ? '' : 'none';
      });
    }

    function openPreview(imgPath, title){
      const modal = document.getElementById('preview');
      document.getElementById('previewTitle').innerText = title || 'Preview';
      const img = document.getElementById('previewImg');
      img.src = imgPath;
      modal.style.display = 'block';
    }
    function closePreview(){
      document.getElementById('preview').style.display = 'none';
      document.getElementById('previewImg').src = '';
    }
    """

    # Build run header elements.
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    summary_links = []
    if assets.get("run_summary_plot"):
        summary_links.append(
            f'<a class="btn" href="{_html_escape(assets["run_summary_plot"])}" target="_blank" rel="noopener">Open run averages plot</a>'
        )
    if assets.get("consolidated_csv"):
        summary_links.append(
            f'<a class="btn secondary" href="{_html_escape(assets["consolidated_csv"])}" target="_blank" rel="noopener">Open consolidated CSV</a>'
        )

    # Construct table rows.
    tr_html = []
    for r in rows:
        # r contains: Filename, Type, Total Reads, Avg Read Length, Avg Base Quality, % Mapped, GC Content %, plots_png, detailed_report_txt
        sample = _html_escape(r.get("Sample", ""))
        raw_reads = r.get("Raw Reads", "N/A")
        aln_reads = r.get("Aligned Reads", "N/A")
        pct_m = r.get("% Mapped", "N/A")
        
        lens = f"{r.get('Raw Avg Len','N/A')} &rarr; {r.get('Prep Avg Len','N/A')} &rarr; {r.get('BAM Avg Len','N/A')}"
        quals = f"{r.get('Raw Avg Q','N/A')} &rarr; {r.get('Prep Avg Q','N/A')} &rarr; {r.get('BAM Avg Q','N/A')}"

        plots = r.get("Plots", "")

        plot_buttons = []
        if plots and os.path.exists(os.path.join(out_dir, plots)):
            plot_buttons.append(
                f'<a class="btn" href="{_html_escape(plots)}">Open plots</a>'
            )
            plot_buttons.append(
                f'<button class="btn secondary" onclick="openPreview(\'{_html_escape(plots)}\', \'{sample} plots\')">Preview</button>'
            )
        else:
            plot_buttons.append('<span class="small muted">No plot</span>')

        tr_html.append(
            "<tr>"
            f"<td class='mono kpi'>{sample}</td>"
            f"<td data-sort='{raw_reads}'>{raw_reads}</td>"
            f"<td data-sort='{aln_reads}'>{aln_reads}</td>"
            f"<td data-sort='{pct_m}'>{pct_m}</td>"
            f"<td>{lens}</td>"
            f"<td>{quals}</td>"
            f"<td>{''.join(plot_buttons)}</td>"
            "</tr>"
        )

    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>{_html_escape(run_name)} Dashboard</title>
<style>{css}</style>
</head>
<body>
  <div class="grid">
    <div>
      <h1>{_html_escape(run_name)}: v32 Sample Progression Dashboard</h1>
      <div class="sub">Generated: {timestamp} | Folder: <span class="mono">{_html_escape(out_dir)}</span></div>
    </div>

    <div class="card">
      <div class="card-header">
        <div class="pill">Run-level summaries</div>
        <div style="display:flex; gap:10px; flex-wrap:wrap;">{''.join(summary_links)}</div>
      </div>
      <div class="controls">
        <input id="searchBox" type="text" placeholder="Search samples, file names, metrics..." oninput="filterTable()"/>
        <span class="small">Tip: click column headers to sort.</span>
      </div>
      <div class="table-wrap">
        <table id="dataTable">
          <thead>
            <tr>
              <th data-base="Sample" onclick="sortTable(0)">Sample</th>
              <th data-base="Raw Reads" onclick="sortTable(1)">Raw Reads</th>
              <th data-base="Aligned Reads" onclick="sortTable(2)">Aligned Reads</th>
              <th data-base="% Mapped" onclick="sortTable(3)">% Mapped</th>
              <th data-base="Lengths (Raw/Prep/BAM)" onclick="sortTable(4)">Lengths (Raw/Prep/BAM)</th>
              <th data-base="Qualities (Raw/Prep/BAM)" onclick="sortTable(5)">Qualities (Raw/Prep/BAM)</th>
              <th data-base="Plots">Progression Plots</th>
            </tr>
          </thead>
          <tbody>
            {''.join(tr_html)}
          </tbody>
        </table>
      </div>
    </div>
  </div>

  <div id="preview" class="preview" onclick="closePreview()">
    <div class="preview-inner" onclick="event.stopPropagation()">
      <div class="preview-head">
        <div id="previewTitle" class="mono">Preview</div>
        <div class="x" onclick="closePreview()">Close</div>
      </div>
      <img id="previewImg" alt="preview"/>
    </div>
  </div>

<script>{js}</script>
</body>
</html>
"""

    with open(index_path, "w", encoding="utf-8") as f:
        f.write(html)

    return index_path


# -----------------------------
# Main Tkinter app
# -----------------------------
class SequenceAnalyzerApp:
    def __init__(self, root):
        _apply_plot_style()

        self.root = root
        self.root.title("Mito Pipeline v32 Output Analyzer (Raw/Prep/BAM)")
        self.root.geometry("860x700")

        style = ttk.Style()
        if "clam" in style.theme_names():
            style.theme_use("clam")

        self.default_font = font.nametofont("TkDefaultFont")
        self.default_font.configure(size=11)
        self.log_font = font.Font(family="Consolas", size=10)

        # Variables
        self.input_path = tk.StringVar()
        self.output_path = tk.StringVar()
        self.run_name = tk.StringVar()

        self.make_dashboard = tk.BooleanVar(value=True)

        self.setup_ui()

    def setup_ui(self):
        main_frame = ttk.Frame(self.root, padding="15")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # --- 1. Data Selection Frame ---
        folder_frame = ttk.LabelFrame(main_frame, text=" 1. Data Selection & Routing ", padding="15")
        folder_frame.pack(fill=tk.X, pady=(0, 15))

        row1 = ttk.Frame(folder_frame)
        row1.pack(fill=tk.X, pady=(0, 10))
        ttk.Label(row1, text="Input Folder:", font=("", 11, "bold"), width=15).pack(side=tk.LEFT)
        ttk.Entry(row1, textvariable=self.input_path, state="readonly").pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 10))
        ttk.Button(row1, text="Browse...", command=self.browse_input).pack(side=tk.LEFT)

        row2 = ttk.Frame(folder_frame)
        row2.pack(fill=tk.X, pady=(0, 10))
        ttk.Label(row2, text="Output Folder:", font=("", 11, "bold"), width=15).pack(side=tk.LEFT)
        ttk.Entry(row2, textvariable=self.output_path, state="readonly").pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 10))
        ttk.Button(row2, text="Browse...", command=self.browse_output).pack(side=tk.LEFT)

        row3 = ttk.Frame(folder_frame)
        row3.pack(fill=tk.X, pady=(0, 10))
        ttk.Label(row3, text="Run Name:", font=("", 11, "bold"), width=15).pack(side=tk.LEFT)
        ttk.Entry(row3, textvariable=self.run_name).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 10))
        ttk.Button(row3, text="Auto", command=self.autofill_run_name).pack(side=tk.LEFT)

        row4 = ttk.Frame(folder_frame)
        row4.pack(fill=tk.X, pady=(5, 0))
        ttk.Checkbutton(
            row4,
            text="Generate HTML dashboard (offline, in output folder)",
            variable=self.make_dashboard,
        ).pack(side=tk.LEFT)

        # --- 2. Control Frame ---
        control_frame = ttk.LabelFrame(main_frame, text=" 2. Run Control ", padding="15")
        control_frame.pack(fill=tk.X, pady=(0, 15))

        self.run_button = ttk.Button(control_frame, text="Run Analysis", command=self.run_analysis)
        self.run_button.pack(side=tk.LEFT)

        self.progress = ttk.Progressbar(control_frame, orient="horizontal", length=520, mode="determinate")
        self.progress.pack(side=tk.LEFT, padx=(15, 0), fill=tk.X, expand=True)

        # --- 3. Log Frame ---
        log_frame = ttk.LabelFrame(main_frame, text=" 3. Log ", padding="10")
        log_frame.pack(fill=tk.BOTH, expand=True)

        self.log_text = tk.Text(log_frame, height=18, wrap="word", font=self.log_font)
        self.log_text.pack(fill=tk.BOTH, expand=True)

    def log(self, msg: str):
        self.log_text.insert(tk.END, msg + "\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()

    def browse_input(self):
        p = filedialog.askdirectory(title="Select input folder containing BAM/FASTQ files")
        if p:
            self.input_path.set(p)

    def browse_output(self):
        p = filedialog.askdirectory(title="Select output folder")
        if p:
            self.output_path.set(p)

    def autofill_run_name(self):
        now = datetime.now().strftime("%Y%m%d_%H%M")
        base = os.path.basename(self.input_path.get().rstrip("/")) if self.input_path.get() else "Run"
        self.run_name.set(f"{base}_{now}")

    def run_analysis(self):
        input_dir = self.input_path.get().strip()
        out_dir = self.output_path.get().strip()
        run = self.run_name.get().strip()

        if not input_dir or not os.path.isdir(input_dir):
            messagebox.showerror("Missing input", "Select a valid input folder.")
            return
        if not out_dir or not os.path.isdir(out_dir):
            messagebox.showerror("Missing output", "Select a valid output folder.")
            return
        if not run:
            messagebox.showerror("Missing run name", "Enter a run name (or click Auto).")
            return

        self.run_button.config(state="disabled")
        self.progress["value"] = 0
        self.log("\nStarting analysis...")

        thread = threading.Thread(target=self.run_analysis_pipeline, args=(input_dir, out_dir, run))
        thread.daemon = True
        thread.start()


    def run_analysis_pipeline(self, input_dir, out_dir, run_name):
        try:
            final_out_dir = os.path.join(out_dir, f"{run_name}_Analysis_Results")
            os.makedirs(final_out_dir, exist_ok=True)
            self.log(f"Routing outputs to: {final_out_dir}")

            # Find sample directories in the v32 pipeline output folder
            sample_dirs = []
            for item in os.listdir(input_dir):
                item_path = os.path.join(input_dir, item)
                if os.path.isdir(item_path):
                    # Ignore aggregate pipeline folders
                    if item in ("sorted_bams", "vcfs", "Final_Pipeline_Results", "fastp_reports", "multiqc_data"):
                        continue
                    
                    # Check if it looks like a sample dir
                    if item.startswith("barcode") or glob.glob(os.path.join(item_path, "*sorted.bam")):
                        sample_dirs.append(item_path)

            if not sample_dirs:
                self.log("No sample directories found in the selected folder. Looking for 'barcode*' or folders with '*sorted.bam'.")
                return

            self.log(f"Found {len(sample_dirs)} sample directories to process.")
            self.progress["maximum"] = len(sample_dirs)

            dashboard_rows = []
            tsv_data = []

            for i, d in enumerate(sorted(sample_dirs)):
                sample = os.path.basename(d)
                self.log(f"Processing sample: {sample}...")

                try:
                    # Enforce strict priority using `or` instead of `+` and sorting by size
                    raw_cand = (
                        glob.glob(os.path.join(d, "*_lenLE1500.fastq")) or 
                        glob.glob(os.path.join(d, "*_lenLE1500_R1.fastq")) or 
                        glob.glob(os.path.join(d, "*_lenLE3000_R1.fastq")) or 
                        glob.glob(os.path.join(d, "*_qsGE10.fastq")) or 
                        glob.glob(os.path.join(d, "*_merged.fastq"))
                    )
                    prep_cand = (
                        glob.glob(os.path.join(d, "*_trim5.fastq")) or 
                        glob.glob(os.path.join(d, "*_qsGE10_trim5.fastq"))
                    )
                    bam_cand = glob.glob(os.path.join(d, "*trimmed_sorted.bam")) or glob.glob(os.path.join(d, "*initial_sorted.bam"))

                    raw_f = raw_cand[0] if raw_cand else None
                    prep_f = prep_cand[0] if prep_cand else None
                    bam_f = bam_cand[0] if bam_cand else None

                    raw_metrics = self.get_fastq_metrics(raw_f) if raw_f else None
                    prep_metrics = self.get_fastq_metrics(prep_f) if prep_f else None
                    bam_metrics = self.get_bam_metrics(bam_f) if bam_f else None

                    plots_img = self.plot_sample_progression(sample, raw_metrics, prep_metrics, bam_metrics, final_out_dir)

                    def get_len_bins(m):
                        if not m: return {"<90": 0, "90-300": 0, ">300": 0}
                        bins = {"<90": 0, "90-300": 0, ">300": 0}
                        for L, count in m.get("len_dist", {}).items():
                            if L < 90: bins["<90"] += count
                            elif 90 <= L <= 300: bins["90-300"] += count
                            else: bins[">300"] += count
                        return bins
                    
                    rb = get_len_bins(raw_metrics)
                    pb = get_len_bins(prep_metrics)
                    bb = get_len_bins(bam_metrics)

                    tsv_data.append({
                        "Sample": sample,
                        "Raw_<90": rb["<90"], "Raw_90-300": rb["90-300"], "Raw_>300": rb[">300"],
                        "Prep_<90": pb["<90"], "Prep_90-300": pb["90-300"], "Prep_>300": pb[">300"],
                        "BAM_<90": bb["<90"], "BAM_90-300": bb["90-300"], "BAM_>300": bb[">300"]
                    })

                    dashboard_rows.append({
                        "Sample": sample,
                        "Raw Reads": raw_metrics["reads"] if raw_metrics else "N/A",
                        "Aligned Reads": bam_metrics["mapped"] if bam_metrics else "N/A",
                        "% Mapped": f"{bam_metrics['perc_mapped']:.2f}%" if bam_metrics else "N/A",
                        "Raw Avg Len": f"{raw_metrics['avg_len']:.1f}" if raw_metrics else "N/A",
                        "Prep Avg Len": f"{prep_metrics['avg_len']:.1f}" if prep_metrics else "N/A",
                        "BAM Avg Len": f"{bam_metrics['avg_len']:.1f}" if bam_metrics else "N/A",
                        "Raw Avg Q": f"{raw_metrics['avg_q']:.1f}" if raw_metrics else "N/A",
                        "Prep Avg Q": f"{prep_metrics['avg_q']:.1f}" if prep_metrics else "N/A",
                        "BAM Avg Q": f"{bam_metrics['avg_q']:.1f}" if bam_metrics else "N/A",
                        "Plots": plots_img
                    })

                    self.log(f"Successfully finished sample: {sample}")
                except Exception as e:
                    self.log(f"ERROR processing sample {sample}: {e}")

                self.progress["value"] = i + 1

            if not dashboard_rows:
                self.log("No samples were successfully processed.")
                return

            if tsv_data:
                tsv_path = os.path.join(final_out_dir, f"{run_name}_Length_Distribution_Summary.tsv")
                with open(tsv_path, 'w') as f:
                    headers = list(tsv_data[0].keys())
                    f.write("\t".join(headers) + "\n")
                    for row in tsv_data:
                        f.write("\t".join(str(row[h]) for h in headers) + "\n")
                self.log(f"Length distribution TSV generated: {tsv_path}")

            # HTML dashboard
            if self.make_dashboard.get():
                dashboard_path = _write_dashboard_html(final_out_dir, run_name, dashboard_rows, {})
                self.log(f"HTML dashboard generated: {dashboard_path}")

            self.log("=== ALL PROCESSING COMPLETE ===")

        finally:
            self.run_button.config(state="normal")

    def get_fastq_metrics(self, path):
        m = {"reads": 0, "bases": 0, "len_dist": defaultdict(int), "qual_dist": defaultdict(int)}
        open_f = gzip.open if path.endswith(".gz") else open
        with open_f(path, "rt", errors="ignore") as f:
            for i, l in enumerate(f):
                mod = i % 4
                if mod == 1:
                    seq = l.strip()
                    m["reads"] += 1
                    m["bases"] += len(seq)
                    m["len_dist"][len(seq)] += 1
                elif mod == 3:
                    for c in l.strip():
                        m["qual_dist"][ord(c)-33] += 1
        m['avg_len'] = m["bases"] / m["reads"] if m["reads"] else 0
        total_q = sum(k*v for k,v in m["qual_dist"].items())
        m['avg_q'] = total_q / m["bases"] if m["bases"] else 0
        return m

    def get_bam_metrics(self, path):
        m = {"reads": 0, "mapped": 0, "len_dist": defaultdict(int), "qual_dist": defaultdict(int)}
        qual_sum = 0
        qual_count = 0
        with pysam.AlignmentFile(path, "rb") as bam:
            for r in bam:
                m["reads"] += 1 # Total reads in BAM
                if not r.is_unmapped and not r.is_secondary and not r.is_supplementary:
                    m["mapped"] += 1
                    if r.query_length:
                        m["len_dist"][r.query_length] += 1
                    if m["mapped"] % 10 == 0 and r.query_qualities:
                        for q in r.query_qualities:
                            m["qual_dist"][q] += 1
                            qual_sum += q
                            qual_count += 1
        tl = sum(l*c for l,c in m["len_dist"].items())
        tr = sum(m["len_dist"].values())
        m['avg_len'] = tl / tr if tr else 0
        m['perc_mapped'] = (m["mapped"] / m["reads"] * 100) if m["reads"] else 0
        m['avg_q'] = qual_sum / qual_count if qual_count > 0 else 0
        return m

    def plot_sample_progression(self, sample, raw, prep, bam, out_dir):
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)
        
        def get_dist(m, key):
            if not m: return [], []
            items = sorted(m[key].items())
            return [k for k,_ in items], [v for _,v in items]

        # 1. Read Length Distribution
        rl_x, rl_y = get_dist(raw, "len_dist")
        pl_x, pl_y = get_dist(prep, "len_dist")
        bl_x, bl_y = get_dist(bam, "len_dist")

        if rl_x: axes[0].plot(rl_x, rl_y, label="Raw FASTQ", marker='o', markersize=4, linewidth=2, alpha=0.8, color="#ff6b6b")
        if pl_x: axes[0].plot(pl_x, pl_y, label="Preprocessed FASTQ", marker='s', markersize=4, linewidth=2, alpha=0.8, color="#ffcc66")
        if bl_x: axes[0].plot(bl_x, bl_y, label="Aligned BAM", marker='^', markersize=4, linewidth=2, alpha=0.8, color="#50c878")
        axes[0].set_xlabel("Read Length (bp)")
        axes[0].set_ylabel("Count")
        axes[0].set_title("Read Length Progression")
        axes[0].legend()
        if any(max(y or [0]) > 5000 for y in [rl_y, pl_y, bl_y]):
            axes[0].set_yscale("log")

        # 2. Base Quality Distribution
        rq_x, rq_y = get_dist(raw, "qual_dist")
        pq_x, pq_y = get_dist(prep, "qual_dist")
        bq_x, bq_y = get_dist(bam, "qual_dist")

        if rq_x: axes[1].plot(rq_x, rq_y, label="Raw FASTQ", marker='o', markersize=4, linewidth=2, alpha=0.8, color="#ff6b6b")
        if pq_x: axes[1].plot(pq_x, pq_y, label="Preprocessed FASTQ", marker='s', markersize=4, linewidth=2, alpha=0.8, color="#ffcc66")
        if bq_x: axes[1].plot(bq_x, bq_y, label="Aligned BAM", marker='^', markersize=4, linewidth=2, alpha=0.8, color="#50c878")
        axes[1].set_xlabel("Phred Quality Score")
        axes[1].set_ylabel("Count")
        axes[1].set_title("Base Quality Progression")
        axes[1].legend()

        out_path = os.path.join(out_dir, f"{sample}_Progression_Plots.png")
        _savefig(out_path)
        return os.path.basename(out_path)

if __name__ == "__main__":
    root = tk.Tk()
    app = SequenceAnalyzerApp(root)
    root.mainloop()
