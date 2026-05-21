#!/usr/bin/env python3
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


# =============================================================================
# THEME — colors pulled from the MARK logo
# =============================================================================
class Theme:
    # Brand palette (sampled from logo)
    NAVY       = "#1B2D5C"   # Primary brand / "M" and "R"
    NAVY_DARK  = "#142346"   # Hover / pressed
    NAVY_LIGHT = "#2A4180"
    TEAL       = "#1D9CB0"   # "K" / accents
    TEAL_DARK  = "#157C8C"
    ORANGE     = "#E8742C"   # "A" / stop / call-to-action
    ORANGE_DARK= "#C45E1F"
    GREEN      = "#3FA34D"   # Success / accents
    PURPLE     = "#7A4FB5"   # Accent

    # Neutrals
    BG         = "#F4F6FA"   # App background
    CARD       = "#FFFFFF"   # Card / frame background
    BORDER     = "#D5DCE5"
    TEXT       = "#1C2333"
    TEXT_MUTED = "#6B7280"
    LOG_BG     = "#0E1726"
    LOG_FG     = "#E6EAF2"

    # Status colors
    SCRIPT_DEF = "#1565C0"   # Blue for "script default" labels
    OVERRIDE   = "#C62828"   # Red for "custom override" labels
    DISABLED   = "#9E9E9E"


class MitoPipelineDashboard:
    def __init__(self, root):
        self.root = root
        self.root.title("MARK — Mitochondrial Amplicon Resolving Kit  v5.2-global-scrollfix")
        self.root.geometry("1180x880")
        self.root.minsize(900, 620)
        self.root.configure(bg=Theme.BG)

        self.script_dir = os.path.dirname(os.path.abspath(__file__))
        self.logo_image = None  # keep a reference so Tk doesn't GC it

        # State management for STOP feature
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

        # Apply theme BEFORE building widgets
        self._apply_theme()

        # Build the modern header (logo + title bar)
        self._build_header()

        # Main content container
        content = tk.Frame(root, bg=Theme.BG)
        content.pack(fill='both', expand=True, padx=14, pady=(0, 10))

        # Activity Log — packed FIRST at the bottom with a fixed height. Packing it
        # before the notebook (and using side='bottom') ensures the notebook gets
        # all the remaining vertical space, so the scrollable tabs always have
        # room to render and the log doesn't squeeze them.
        log_wrap = tk.Frame(content, bg=Theme.BG)
        log_wrap.pack(side='bottom', fill='x', expand=False, pady=(8, 0))

        log_header = tk.Frame(log_wrap, bg=Theme.NAVY, height=32)
        log_header.pack(fill='x')
        log_header.pack_propagate(False)
        tk.Label(log_header, text="●  ACTIVITY LOG",
                 bg=Theme.NAVY, fg="#FFFFFF",
                 font=("Segoe UI", 9, "bold")).pack(side='left', padx=14)
        tk.Label(log_header, text="real-time pipeline output",
                 bg=Theme.NAVY, fg="#9DB0CC",
                 font=("Segoe UI", 8, "italic")).pack(side='left', padx=4)

        self.log_text = scrolledtext.ScrolledText(
            log_wrap, height=4, state='disabled',
            font=("Consolas", 9),
            bg=Theme.LOG_BG, fg=Theme.LOG_FG,
            insertbackground=Theme.LOG_FG,
            relief='flat', borderwidth=0,
            padx=10, pady=6
        )
        self.log_text.pack(fill='x', expand=False)

        # Notebook — packed AFTER the log (in code) but it will fill the remaining
        # top area because the log was placed at side='bottom'.
        self.notebook = ttk.Notebook(content, style="Brand.TNotebook")
        self.notebook.pack(expand=True, fill='both', pady=(8, 0))

        self.tab_pipeline = ttk.Frame(self.notebook, style="Card.TFrame")
        self.tab_post = ttk.Frame(self.notebook, style="Card.TFrame")

        self.notebook.add(self.tab_pipeline, text="  Step 1 · Run Pipeline  ")
        self.notebook.add(self.tab_post, text="  Step 2 · Post-Processing  ")

        self.init_pipeline_tab()
        self.init_postproc_tab()

        # When the user switches tabs, force a scrollregion recompute on every
        # scrollable canvas — geometry only finalizes once a tab is shown.
        def _on_tab_changed(_event=None):
            for c in getattr(self, "_scroll_canvases", []):
                try:
                    c.update_idletasks()
                    bbox = c.bbox('all')
                    if bbox:
                        c.configure(scrollregion=(0, 0, bbox[2], bbox[3] + 4))
                except tk.TclError:
                    pass
        self.notebook.bind('<<NotebookTabChanged>>', _on_tab_changed)

        # Surface header/logo status now that the log widget exists.
        if getattr(self, "_logo_status", None):
            self.log(self._logo_status)

        self.refresh_config_ui()

    # =========================================================================
    # THEME / STYLE
    # =========================================================================
    def _apply_theme(self):
        style = ttk.Style()
        # 'clam' gives us the most styling control across platforms
        try:
            style.theme_use('clam')
        except tk.TclError:
            pass

        # --- Base frames ---
        style.configure("TFrame", background=Theme.BG)
        style.configure("Card.TFrame", background=Theme.CARD)
        style.configure("Header.TFrame", background=Theme.NAVY)

        # --- Labels ---
        style.configure("TLabel", background=Theme.BG, foreground=Theme.TEXT,
                        font=("Segoe UI", 10))
        style.configure("Card.TLabel", background=Theme.CARD, foreground=Theme.TEXT,
                        font=("Segoe UI", 10))
        style.configure("Muted.TLabel", background=Theme.CARD, foreground=Theme.TEXT_MUTED,
                        font=("Segoe UI", 9, "italic"))
        style.configure("FieldLabel.TLabel", background=Theme.CARD, foreground=Theme.TEXT,
                        font=("Segoe UI", 10))
        style.configure("Bold.TLabel", background=Theme.CARD, foreground=Theme.NAVY,
                        font=("Segoe UI", 10, "bold"))

        # --- LabelFrame (cards) ---
        style.configure("Card.TLabelframe", background=Theme.CARD,
                        bordercolor=Theme.BORDER, relief='solid', borderwidth=1,
                        padding=12)
        style.configure("Card.TLabelframe.Label", background=Theme.CARD,
                        foreground=Theme.NAVY, font=("Segoe UI", 10, "bold"))

        # --- Notebook ---
        style.configure("Brand.TNotebook", background=Theme.BG, borderwidth=0,
                        tabmargins=[2, 6, 2, 0])
        style.configure("Brand.TNotebook.Tab",
                        padding=[18, 9],
                        font=("Segoe UI", 10, "bold"),
                        background="#DCE3EE",
                        foreground=Theme.NAVY,
                        borderwidth=0)
        style.map("Brand.TNotebook.Tab",
                  background=[("selected", Theme.CARD),
                              ("active", "#E8EDF5")],
                  foreground=[("selected", Theme.NAVY),
                              ("active", Theme.NAVY)])

        # --- Entries ---
        style.configure("TEntry",
                        fieldbackground=Theme.CARD,
                        bordercolor=Theme.BORDER,
                        lightcolor=Theme.BORDER,
                        darkcolor=Theme.BORDER,
                        foreground=Theme.TEXT,
                        padding=5)
        style.map("TEntry",
                  bordercolor=[("focus", Theme.TEAL)],
                  lightcolor=[("focus", Theme.TEAL)],
                  darkcolor=[("focus", Theme.TEAL)])

        # --- Default Button (secondary / browse) ---
        style.configure("TButton",
                        font=("Segoe UI", 9),
                        background="#E8EDF5",
                        foreground=Theme.NAVY,
                        bordercolor=Theme.BORDER,
                        focusthickness=0,
                        padding=(10, 5),
                        relief='flat')
        style.map("TButton",
                  background=[("active", "#D5DEEC"), ("pressed", "#C2CDE0")],
                  foreground=[("active", Theme.NAVY), ("disabled", Theme.DISABLED)])

        # --- Primary action button (navy) ---
        style.configure("Primary.TButton",
                        font=("Segoe UI", 10, "bold"),
                        background=Theme.NAVY,
                        foreground="#FFFFFF",
                        padding=(18, 11),
                        relief='flat',
                        focusthickness=0,
                        borderwidth=0)
        style.map("Primary.TButton",
                  background=[("active", Theme.NAVY_DARK),
                              ("pressed", Theme.NAVY_DARK),
                              ("disabled", "#B6BFD0")],
                  foreground=[("disabled", "#FFFFFF")])

        # --- Secondary action (teal) ---
        style.configure("Accent.TButton",
                        font=("Segoe UI", 10, "bold"),
                        background=Theme.TEAL,
                        foreground="#FFFFFF",
                        padding=(18, 11),
                        relief='flat',
                        focusthickness=0,
                        borderwidth=0)
        style.map("Accent.TButton",
                  background=[("active", Theme.TEAL_DARK),
                              ("pressed", Theme.TEAL_DARK),
                              ("disabled", "#B7D5DC")],
                  foreground=[("disabled", "#FFFFFF")])

        # --- Stop button (orange) ---
        style.configure("Danger.TButton",
                        font=("Segoe UI", 10, "bold"),
                        background=Theme.ORANGE,
                        foreground="#FFFFFF",
                        padding=(18, 11),
                        relief='flat',
                        focusthickness=0,
                        borderwidth=0)
        style.map("Danger.TButton",
                  background=[("active", Theme.ORANGE_DARK),
                              ("pressed", Theme.ORANGE_DARK),
                              ("disabled", "#E8C8B5")],
                  foreground=[("disabled", "#FFFFFF")])

        # --- Step button (small accent, used in post-processing) ---
        style.configure("Step.TButton",
                        font=("Segoe UI", 9, "bold"),
                        background=Theme.GREEN,
                        foreground="#FFFFFF",
                        padding=(12, 6),
                        relief='flat',
                        focusthickness=0,
                        borderwidth=0)
        style.map("Step.TButton",
                  background=[("active", "#358A41"),
                              ("pressed", "#358A41")])

        # --- Checkbutton ---
        style.configure("TCheckbutton",
                        background=Theme.CARD,
                        foreground=Theme.TEXT,
                        focuscolor=Theme.CARD)
        style.map("TCheckbutton",
                  background=[("active", Theme.CARD)])

        # --- Progress bar ---
        style.configure("Brand.Horizontal.TProgressbar",
                        troughcolor="#E2E7F0",
                        bordercolor="#E2E7F0",
                        background=Theme.TEAL,
                        lightcolor=Theme.TEAL,
                        darkcolor=Theme.TEAL,
                        thickness=18)

    # =========================================================================
    # HEADER (logo + title bar)
    # =========================================================================
    def _build_header(self):
        header = tk.Frame(self.root, bg=Theme.NAVY, height=110)
        header.pack(fill='x', side='top')
        header.pack_propagate(False)

        # Find the logo. Search several locations and a few common filenames so the
        # user doesn't have to rename anything.
        candidate_names = [
            "mark_logo.png", "MARK_logo.png", "mark.png", "MARK.png", "logo.png",
            "F0B5C749-F261-4A31-B164-D580FF584E28.png",
        ]
        candidate_dirs = [
            self.script_dir,
            os.getcwd(),
            os.path.join(self.script_dir, "assets"),
            os.path.join(self.script_dir, "resources"),
        ]
        logo_path = None
        for d in candidate_dirs:
            for name in candidate_names:
                p = os.path.join(d, name)
                if os.path.exists(p):
                    logo_path = p
                    break
            if logo_path:
                break

        logo_loaded = False
        load_error = None

        if logo_path:
            try:
                # Prefer Pillow — clean resizing and PNG transparency support.
                try:
                    from PIL import Image, ImageTk
                    img = Image.open(logo_path).convert("RGBA")
                    target_h = 90
                    ratio = target_h / img.height
                    target_w = int(img.width * ratio)
                    img = img.resize((target_w, target_h), Image.LANCZOS)
                    bg = Image.new("RGBA", img.size, Theme.NAVY)
                    bg.paste(img, (0, 0), img)
                    self.logo_image = ImageTk.PhotoImage(bg)
                except ImportError:
                    # Fallback: native Tk. Note: tk.PhotoImage only handles PNGs
                    # on Tk 8.6+; some platforms still fail on certain PNGs.
                    self.logo_image = tk.PhotoImage(file=logo_path)
                    w = self.logo_image.width()
                    if w > 280:
                        factor = max(2, w // 260)
                        self.logo_image = self.logo_image.subsample(factor, factor)

                tk.Label(header, image=self.logo_image, bg=Theme.NAVY,
                         borderwidth=0).pack(side='left', padx=(18, 10), pady=8)
                logo_loaded = True
                # Defer the log call until after self.log_text exists — store and emit later.
                self._logo_status = f"Logo loaded from: {logo_path}"
            except Exception as e:
                load_error = f"{type(e).__name__}: {e}"
                self._logo_status = (
                    f"Logo found at {logo_path} but failed to load ({load_error}). "
                    "Tip: pip install Pillow"
                )
        else:
            self._logo_status = (
                "Logo not found. Place 'mark_logo.png' next to this script "
                f"(searched: {self.script_dir})."
            )

        # Text block on the right of the logo
        text_block = tk.Frame(header, bg=Theme.NAVY)
        text_block.pack(side='left', padx=(6 if logo_loaded else 22, 0), pady=14)

        # Brand wordmark (text fallback if no logo)
        if not logo_loaded:
            brand = tk.Frame(text_block, bg=Theme.NAVY)
            brand.pack(anchor='w')
            tk.Label(brand, text="M", bg=Theme.NAVY, fg="#FFFFFF",
                     font=("Segoe UI", 28, "bold")).pack(side='left')
            tk.Label(brand, text="A", bg=Theme.NAVY, fg=Theme.ORANGE,
                     font=("Segoe UI", 28, "bold")).pack(side='left')
            tk.Label(brand, text="R", bg=Theme.NAVY, fg="#FFFFFF",
                     font=("Segoe UI", 28, "bold")).pack(side='left')
            tk.Label(brand, text="K", bg=Theme.NAVY, fg=Theme.TEAL,
                     font=("Segoe UI", 28, "bold")).pack(side='left')

        tk.Label(text_block,
                 text="Mitochondrial Amplicon Resolving Kit",
                 bg=Theme.NAVY, fg="#FFFFFF",
                 font=("Segoe UI", 14, "bold")).pack(anchor='w')
        tk.Label(text_block,
                 text="ONT / Illumina mitochondrial sequencing pipeline   ·   v5.0",
                 bg=Theme.NAVY, fg="#9DB0CC",
                 font=("Segoe UI", 9)).pack(anchor='w', pady=(2, 0))

        # Right-side status pill
        right = tk.Frame(header, bg=Theme.NAVY)
        right.pack(side='right', padx=22)
        tk.Label(right, text="Dashboard",
                 bg=Theme.NAVY, fg="#9DB0CC",
                 font=("Segoe UI", 8, "bold")).pack(anchor='e')
        tk.Label(right, text="● Ready",
                 bg=Theme.NAVY, fg=Theme.GREEN,
                 font=("Segoe UI", 10, "bold")).pack(anchor='e', pady=(2, 0))

        # Thin accent stripe under the header
        stripe = tk.Frame(self.root, bg=Theme.TEAL, height=3)
        stripe.pack(fill='x', side='top')

    # =========================================================================
    # LOGGING / PROGRESS / UI STATE
    # =========================================================================
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
    # SCROLLABLE TAB HELPER
    # =========================================================================
    def _make_scrollable_tab(self, parent):
        """
        Wrap a tab's contents in a vertically scrollable canvas.

        This version forces wheel/trackpad scrolling through a custom Tk bindtag
        that is inserted ahead of normal widget/class bindings. That is important
        on macOS because ttk widgets can consume <MouseWheel> before root-level
        or bind_all handlers ever see the event.
        """
        container = tk.Frame(parent, bg=Theme.CARD)
        container.pack(fill='both', expand=True)

        canvas = tk.Canvas(container, bg=Theme.CARD, highlightthickness=0, borderwidth=0)
        vscroll = ttk.Scrollbar(container, orient='vertical', command=canvas.yview)
        canvas.configure(yscrollcommand=vscroll.set)

        canvas.pack(side='left', fill='both', expand=True)
        vscroll.pack(side='right', fill='y')

        inner = tk.Frame(canvas, bg=Theme.CARD)
        inner_id = canvas.create_window((0, 0), window=inner, anchor='nw')

        canvas._is_scroll_canvas = True
        inner._parent_scroll_canvas = canvas

        if not hasattr(self, "_scroll_canvases"):
            self._scroll_canvases = []
        self._scroll_canvases.append(canvas)

        def _activate_canvas(_event=None, c=canvas):
            self._active_wheel_canvas = c

        canvas.bind('<Enter>', _activate_canvas, add='+')
        inner.bind('<Enter>', _activate_canvas, add='+')

        def _on_canvas_configure(event):
            try:
                canvas.itemconfigure(inner_id, width=event.width)
                self.root.after_idle(_update_scrollregion)
            except tk.TclError:
                pass
        canvas.bind('<Configure>', _on_canvas_configure, add='+')

        def _update_scrollregion(_event=None):
            try:
                inner.update_idletasks()
                bbox = canvas.bbox('all')
                if bbox:
                    canvas.configure(
                        scrollregion=(0, 0, max(bbox[2], canvas.winfo_width()), bbox[3] + 80)
                    )
            except tk.TclError:
                pass

        def _apply_scroll_bindtags(widget=None):
            """Attach the dashboard scroll bindtag to every non-text child widget."""
            if widget is None:
                widget = inner
            try:
                cls = widget.winfo_class()
                if cls not in ('Text', 'ScrolledText'):
                    tags = list(widget.bindtags())
                    tag = 'MARKDashboardScroll'
                    if tag in tags:
                        tags.remove(tag)
                    widget.bindtags((tag, *tags))
            except tk.TclError:
                return
            for child in widget.winfo_children():
                _apply_scroll_bindtags(child)

        def _refresh_everything(_event=None):
            _update_scrollregion()
            self.root.after_idle(_apply_scroll_bindtags)

        inner.bind('<Configure>', _refresh_everything, add='+')
        canvas.bind('<Map>', lambda e: self.root.after_idle(_refresh_everything), add='+')
        canvas.bind('<Visibility>', lambda e: self.root.after_idle(_refresh_everything), add='+')

        if not getattr(self, "_wheel_handler_installed", False):
            self._install_global_wheel_handler()
            self._wheel_handler_installed = True

        # Apply several times because some controls are created after the tab
        # frame itself exists.
        for delay in (50, 250, 750, 1500, 3000):
            self.root.after(delay, _refresh_everything)

        inner.refresh_scroll_bindings = _refresh_everything
        return inner

    def _install_global_wheel_handler(self):
        """
        Wheel/trackpad scrolling for all scrollable dashboard tabs.

        Two layers are used:
        1. A custom bindtag added to child widgets inside scrollable tabs. This
           fires before widget/class bindings and prevents ttk widgets from
           swallowing trackpad events.
        2. Root-level fallbacks for events that originate from empty canvas space
           or unusual platform-specific widgets.
        """
        self._wheel_accum = 0.0
        self._active_wheel_canvas = None

        def _active_scroll_canvas():
            try:
                active_tab = self.notebook.nametowidget(self.notebook.select())
            except Exception:
                active_tab = None
            for c in getattr(self, "_scroll_canvases", []):
                try:
                    if active_tab is None or str(c).startswith(str(active_tab)):
                        return c
                except Exception:
                    continue
            return None

        def _can_scroll(canvas):
            try:
                bbox = canvas.bbox('all')
                if bbox is None:
                    return False
                return (bbox[3] - bbox[1]) > canvas.winfo_height()
            except tk.TclError:
                return False

        def _scroll(canvas, units):
            try:
                if canvas is None or units == 0 or not _can_scroll(canvas):
                    return "break"
                canvas.yview_scroll(units, 'units')
                return "break"
            except tk.TclError:
                return "break"

        def _wheel_units(event):
            delta = getattr(event, 'delta', 0)
            if delta == 0:
                return 0

            # Windows wheels often use +/-120. macOS trackpads often send small
            # integer deltas. Convert both into canvas units and accumulate the
            # small trackpad deltas so a light gesture still moves the page.
            if abs(delta) >= 120:
                return -int(delta / 120)
            if abs(delta) >= 2:
                return -1 if delta > 0 else 1

            self._wheel_accum += -delta
            units = int(self._wheel_accum)
            if units != 0:
                self._wheel_accum -= units
            return units

        def _target_canvas(event=None):
            # Always scroll the visible dashboard tab. This is the behaviour the
            # UI needs: when the pointer is anywhere on the page, the page scrolls.
            return _active_scroll_canvas() or getattr(self, "_active_wheel_canvas", None)

        def _on_mousewheel(event):
            return _scroll(_target_canvas(event), _wheel_units(event))

        def _on_linux_up(event):
            return _scroll(_target_canvas(event), -3)

        def _on_linux_down(event):
            return _scroll(_target_canvas(event), 3)

        def _key_scroll(units):
            return _scroll(_active_scroll_canvas(), units)

        def _safe_key(units):
            def handler(_event):
                focused = self.root.focus_get()
                if focused is not None:
                    try:
                        cls = focused.winfo_class()
                        if cls in ('Entry', 'TEntry', 'Text', 'ScrolledText', 'TCombobox'):
                            return
                    except tk.TclError:
                        pass
                return _key_scroll(units)
            return handler

        # Bind to the custom tag that is inserted into every child widget inside
        # each scrollable tab. This is the key fix for macOS trackpads.
        self.root.bind_class('MARKDashboardScroll', '<MouseWheel>', _on_mousewheel)
        self.root.bind_class('MARKDashboardScroll', '<Button-4>', _on_linux_up)
        self.root.bind_class('MARKDashboardScroll', '<Button-5>', _on_linux_down)

        # Fallbacks for empty canvas/background areas and platform edge cases.
        self.root.bind_all('<MouseWheel>', _on_mousewheel, add='+')
        self.root.bind_all('<Button-4>', _on_linux_up, add='+')
        self.root.bind_all('<Button-5>', _on_linux_down, add='+')

        self.root.bind_all('<Prior>', _safe_key(-5), add='+')
        self.root.bind_all('<Next>', _safe_key(5), add='+')
        self.root.bind_all('<Home>', _safe_key(-999), add='+')
        self.root.bind_all('<End>', _safe_key(999), add='+')

    # =========================================================================
    # TAB 1: MAIN PIPELINE
    # =========================================================================
    def init_pipeline_tab(self):
        # Scrollable container — content can exceed window height.
        scroll_inner = self._make_scrollable_tab(self.tab_pipeline)
        self._pipeline_scroll_inner = scroll_inner  # for re-applying scroll tags
        # Padded host frame for the actual sections.
        outer = tk.Frame(scroll_inner, bg=Theme.CARD)
        outer.pack(fill='both', expand=True, padx=14, pady=14)

        # -- 1. Pipeline script ----
        script_frame = ttk.LabelFrame(outer, text="  1. Pipeline Script  ",
                                      style="Card.TLabelframe")
        script_frame.pack(fill='x', pady=(0, 10))
        script_frame.columnconfigure(1, weight=1)

        ttk.Label(script_frame, text="Pipeline Script (.sh):",
                  style="FieldLabel.TLabel").grid(row=0, column=0, sticky='e', padx=6, pady=8)
        self.script_path = tk.StringVar(value=self.get_def_path("MARK.sh"))
        ttk.Entry(script_frame, textvariable=self.script_path).grid(row=0, column=1, sticky='ew', padx=6)
        ttk.Button(script_frame, text="Browse & Load",
                   command=self.browse_script).grid(row=0, column=2, padx=6)

        # -- 2. Parameter configuration ----
        self.config_frame = ttk.LabelFrame(outer, text="  2. Parameter Configuration  ",
                                           style="Card.TLabelframe")
        self.config_frame.pack(fill='x', pady=(0, 10))

        btn_frame = tk.Frame(self.config_frame, bg=Theme.CARD)
        btn_frame.pack(fill='x', pady=(0, 6))
        ttk.Button(btn_frame, text="Reset to Script Defaults",
                   command=self.set_all_script_defaults).pack(side='left')
        ttk.Label(btn_frame,
                  text="Tick a row to override the script's default value.",
                  style="Muted.TLabel").pack(side='left', padx=14)

        self.grid_frame = tk.Frame(self.config_frame, bg=Theme.CARD)
        self.grid_frame.pack(fill='x', pady=(4, 0))

        # -- 3. Input & output ----
        io_frame = ttk.LabelFrame(outer, text="  3. Input & Output Selection  ",
                                  style="Card.TLabelframe")
        io_frame.pack(fill='x', pady=(0, 10))
        io_frame.columnconfigure(1, weight=1)

        def row_entry(parent, r, label, var, browse_cmd, kind='file', filetypes=None):
            ttk.Label(parent, text=label, style="FieldLabel.TLabel").grid(
                row=r, column=0, sticky='e', padx=6, pady=5)
            ttk.Entry(parent, textvariable=var).grid(row=r, column=1, sticky='ew', padx=6, pady=5)
            ttk.Button(parent, text="Browse", command=browse_cmd).grid(
                row=r, column=2, padx=6, pady=5)

        self.ref_path = tk.StringVar(value=self.get_def_path("linearized_mtdna.fasta"))
        row_entry(io_frame, 0, "Linear Ref (FASTA):", self.ref_path,
                  lambda: self.browse_file(self.ref_path, [("FASTA", "*.fasta *.fa")]))

        self.bed_path = tk.StringVar(value=self.get_def_path("linearized_regions.bed"))
        row_entry(io_frame, 1, "Regions BED:", self.bed_path,
                  lambda: self.browse_file(self.bed_path, [("BED", "*.bed")]))

        self.adapter_path = tk.StringVar(value=self.get_def_path("Updated_Adapter_Primer_List_Cutadapt_cleaned.txt"))
        row_entry(io_frame, 2, "Adapter File:", self.adapter_path,
                  lambda: self.browse_file(self.adapter_path, [("Text", "*.txt")]))

        self.input_folder = tk.StringVar()
        row_entry(io_frame, 3, "Input FASTQ Folder:", self.input_folder,
                  lambda: self.browse_folder(self.input_folder))

        self.output_folder = tk.StringVar()
        row_entry(io_frame, 4, "Output Base Dir:", self.output_folder,
                  lambda: self.browse_folder(self.output_folder))
        ttk.Label(io_frame,
                  text="Optional — if left empty, outputs are saved next to the input folder.",
                  style="Muted.TLabel").grid(row=5, column=1, sticky='w', padx=6, pady=(0, 6))

        ttk.Label(io_frame, text="Custom Run Name:",
                  style="FieldLabel.TLabel").grid(row=6, column=0, sticky='e', padx=6, pady=5)
        self.custom_run_name = tk.StringVar()
        ttk.Entry(io_frame, textvariable=self.custom_run_name).grid(
            row=6, column=1, sticky='ew', padx=6, pady=5)
        ttk.Label(io_frame,
                  text="Optional — if left empty, the pipeline uses its default naming convention.",
                  style="Muted.TLabel").grid(row=7, column=1, sticky='w', padx=6, pady=(0, 6))

        ttk.Button(io_frame, text="Clear Fields",
                   command=self.clear_pipeline_fields).grid(row=8, column=1, pady=(8, 0), sticky='w', padx=6)

        # -- Run buttons (prominent, colored) ----
        run_btn_frame = tk.Frame(outer, bg=Theme.CARD)
        run_btn_frame.pack(pady=(6, 8))

        self.btn_run_only = ttk.Button(run_btn_frame, text="▶  RUN PIPELINE ONLY",
                                       style="Primary.TButton",
                                       command=lambda: self.start_thread(self.run_pipeline_wrapper))
        self.btn_run_only.pack(side='left', padx=8)

        self.btn_run_auto = ttk.Button(run_btn_frame, text="⚡  RUN PIPELINE + POST-PROCESSING (AUTO)",
                                       style="Accent.TButton",
                                       command=lambda: self.start_thread(self.run_automated_cycle))
        self.btn_run_auto.pack(side='left', padx=8)

        self.btn_stop = ttk.Button(run_btn_frame, text="■  STOP RUN",
                                   style="Danger.TButton",
                                   command=self.stop_pipeline, state='disabled')
        self.btn_stop.pack(side='left', padx=8)

        # -- Progress ----
        self.progress_frame = ttk.LabelFrame(outer, text="  Progress  ",
                                             style="Card.TLabelframe")
        self.progress_frame.pack(fill='x', pady=(2, 0))
        self.progress_bar = ttk.Progressbar(self.progress_frame, orient='horizontal',
                                            mode='determinate', length=400,
                                            style="Brand.Horizontal.TProgressbar")
        self.progress_bar.pack(fill='x', padx=4, pady=(4, 8))
        self.progress_var = tk.StringVar(value="Idle")
        ttk.Label(self.progress_frame, textvariable=self.progress_var,
                  style="Bold.TLabel").pack(pady=(0, 4))

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
            tk.Label(self.grid_frame, text=text,
                     font=("Segoe UI", 9, "bold"),
                     bg=Theme.CARD, fg=Theme.NAVY).grid(row=0, column=col, padx=10, pady=(4, 8), sticky='w')

        # subtle divider
        sep = tk.Frame(self.grid_frame, bg=Theme.BORDER, height=1)
        sep.grid(row=1, column=0, columnspan=5, sticky='ew', pady=(0, 6))

        row = 2
        for var, dash_def in self.dash_defaults.items():
            script_val = script_defaults.get(var, "Not Found")
            is_hardcoded = "(Hardcoded)" in script_val
            clean_script_val = script_val.replace(" (Hardcoded)", "")

            # zebra striping for readability
            row_bg = Theme.CARD if (row % 2 == 0) else "#FAFBFD"

            # parameter name
            tk.Label(self.grid_frame, text=f"{var}",
                     font=("Segoe UI", 9, "bold"),
                     bg=row_bg, fg=Theme.TEXT).grid(row=row, column=0, sticky='ew', padx=6, pady=2, ipady=3)

            color = Theme.DISABLED if is_hardcoded else (Theme.SCRIPT_DEF if clean_script_val != "Not Found" else Theme.DISABLED)
            tk.Label(self.grid_frame, text=script_val,
                     foreground=color, font=("Segoe UI", 9, "bold"),
                     bg=row_bg).grid(row=row, column=1, padx=6, sticky='ew', ipady=3)

            entry = ttk.Entry(self.grid_frame,
                              textvariable=self.var_states[var]['entry_var'], width=12)
            entry.grid(row=row, column=2, padx=6, sticky='w', pady=2)
            self.var_states[var]['entry_widget'] = entry
            self.var_states[var]['entry_var'].trace_add("write", lambda *args: self.update_inline_summary())

            chk = ttk.Checkbutton(self.grid_frame, variable=self.var_states[var]['override_var'],
                                  command=lambda v=var: [self.toggle_entry_state(v), self.update_inline_summary()])
            chk.grid(row=row, column=3, padx=6, sticky='w')
            self.var_states[var]['chk_widget'] = chk

            final_lbl = tk.Label(self.grid_frame, text="",
                                 font=("Consolas", 10, "bold"),
                                 bg=row_bg)
            final_lbl.grid(row=row, column=4, padx=14, sticky='ew', ipady=3)
            self.var_states[var]['final_label'] = final_lbl

            # fill row background gaps
            self.grid_frame.grid_columnconfigure(0, weight=0)
            self.grid_frame.grid_columnconfigure(4, weight=1)

            if clean_script_val == "Not Found" or is_hardcoded:
                self.var_states[var]['override_var'].set(False)
                chk.config(state='disabled')
                entry.config(state='disabled')
            else:
                self.toggle_entry_state(var)

            row += 1

        self.update_inline_summary()

        # Re-apply scroll bindtags to the newly-built parameter grid widgets so
        # trackpad/wheel scrolling works when the cursor is over them.
        if hasattr(self, "_pipeline_scroll_inner") and hasattr(
            self._pipeline_scroll_inner, "refresh_scroll_bindings"
        ):
            self._pipeline_scroll_inner.refresh_scroll_bindings()

    def update_inline_summary(self):
        script_defaults = self.parse_script_defaults(self.script_path.get())
        for var, state in self.var_states.items():
            s_val_raw = script_defaults.get(var, "Not Found")
            clean_s_val = s_val_raw.replace(" (Hardcoded)", "")
            val = state['entry_var'].get()

            if clean_s_val == "Not Found":
                state['final_label'].config(text="Disabled (Not found in script)", foreground=Theme.DISABLED)
            elif state['override_var'].get():
                state['final_label'].config(text=f"➜ {val}  (Custom Override)", foreground=Theme.OVERRIDE)
            else:
                state['final_label'].config(text=f"➜ {clean_s_val}  (Script Default)", foreground=Theme.SCRIPT_DEF)

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
        scroll_inner = self._make_scrollable_tab(self.tab_post)
        outer = tk.Frame(scroll_inner, bg=Theme.CARD)
        outer.pack(fill='both', expand=True, padx=14, pady=14)

        intro = tk.Label(outer,
                         text="Run any individual step, or use AUTO from Step 1 to run the full chain.",
                         bg=Theme.CARD, fg=Theme.TEXT_MUTED,
                         font=("Segoe UI", 9, "italic"))
        intro.pack(anchor='w', pady=(0, 8))

        # 1. Organize Raw Output
        frame_raw = ttk.LabelFrame(outer, text="  1. Organize Raw Output  ",
                                   style="Card.TLabelframe")
        frame_raw.pack(fill='x', pady=(0, 8))
        frame_raw.columnconfigure(1, weight=1)

        ttk.Label(frame_raw, text="Run Folder:", style="FieldLabel.TLabel").grid(
            row=0, column=0, sticky='e', padx=6, pady=6)
        self.raw_main_dir = tk.StringVar()
        ttk.Entry(frame_raw, textvariable=self.raw_main_dir).grid(
            row=0, column=1, sticky='ew', padx=6, pady=6)
        ttk.Button(frame_raw, text="Browse",
                   command=lambda: self.browse_folder(self.raw_main_dir)).grid(row=0, column=2, padx=6)
        ttk.Button(frame_raw, text="Run Step 1", style="Step.TButton",
                   command=lambda: self.start_thread(self.collect_raw_wrapper)).grid(
            row=0, column=3, padx=(10, 6))

        # 2. Linear VCF Correction
        frame_vcf = ttk.LabelFrame(outer, text="  2. Linear VCF Correction  ",
                                   style="Card.TLabelframe")
        frame_vcf.pack(fill='x', pady=(0, 8))
        frame_vcf.columnconfigure(1, weight=1)

        self.vcf_mod_ref = tk.StringVar(value=self.get_def_path("linearized_mtdna.fasta"))
        self.vcf_orig_ref = tk.StringVar(value=self.get_def_path("rCRS.fasta"))
        self.vcf_target_dir = tk.StringVar()

        ttk.Label(frame_vcf, text="Linearized Ref:", style="FieldLabel.TLabel").grid(
            row=0, column=0, sticky='e', padx=6, pady=4)
        ttk.Entry(frame_vcf, textvariable=self.vcf_mod_ref).grid(row=0, column=1, sticky='ew', padx=6)
        ttk.Button(frame_vcf, text="Browse",
                   command=lambda: self.browse_file(self.vcf_mod_ref, [("FASTA", "*.fasta *.fa")])).grid(
            row=0, column=2, padx=6)

        ttk.Label(frame_vcf, text="Original Circular Ref:", style="FieldLabel.TLabel").grid(
            row=1, column=0, sticky='e', padx=6, pady=4)
        ttk.Entry(frame_vcf, textvariable=self.vcf_orig_ref).grid(row=1, column=1, sticky='ew', padx=6)
        ttk.Button(frame_vcf, text="Browse",
                   command=lambda: self.browse_file(self.vcf_orig_ref, [("FASTA", "*.fasta *.fa")])).grid(
            row=1, column=2, padx=6)

        ttk.Label(frame_vcf, text="VCF Folder:", style="FieldLabel.TLabel").grid(
            row=2, column=0, sticky='e', padx=6, pady=4)
        ttk.Entry(frame_vcf, textvariable=self.vcf_target_dir).grid(row=2, column=1, sticky='ew', padx=6)
        ttk.Button(frame_vcf, text="Browse",
                   command=lambda: self.browse_folder(self.vcf_target_dir)).grid(row=2, column=2, padx=6)

        ttk.Button(frame_vcf, text="Run Step 2", style="Step.TButton",
                   command=lambda: self.start_thread(self.vcf_correct_wrapper)).grid(
            row=3, column=1, pady=(8, 2))

        # 3. BAM Header Cleaning
        frame_bam = ttk.LabelFrame(outer, text="  3. BAM Header Cleaning  ",
                                   style="Card.TLabelframe")
        frame_bam.pack(fill='x', pady=(0, 8))
        frame_bam.columnconfigure(1, weight=1)

        self.bam_target_dir = tk.StringVar()
        ttk.Label(frame_bam, text="BAM Folder:", style="FieldLabel.TLabel").grid(
            row=0, column=0, sticky='e', padx=6, pady=6)
        ttk.Entry(frame_bam, textvariable=self.bam_target_dir).grid(
            row=0, column=1, sticky='ew', padx=6, pady=6)
        ttk.Button(frame_bam, text="Browse",
                   command=lambda: self.browse_folder(self.bam_target_dir)).grid(row=0, column=2, padx=6)
        ttk.Button(frame_bam, text="Run Step 3", style="Step.TButton",
                   command=lambda: self.start_thread(self.bam_clean_wrapper)).grid(
            row=0, column=3, padx=(10, 6))

        # 4. Final Collection
        frame_final = ttk.LabelFrame(outer, text="  4. Final Collection  ",
                                     style="Card.TLabelframe")
        frame_final.pack(fill='x', pady=(0, 8))
        frame_final.columnconfigure(1, weight=1)

        self.final_main_dir = tk.StringVar()
        ttk.Label(frame_final, text="Output Parent Folder:", style="FieldLabel.TLabel").grid(
            row=0, column=0, sticky='e', padx=6, pady=6)
        ttk.Entry(frame_final, textvariable=self.final_main_dir).grid(
            row=0, column=1, sticky='ew', padx=6, pady=6)
        ttk.Button(frame_final, text="Browse",
                   command=lambda: self.browse_folder(self.final_main_dir)).grid(row=0, column=2, padx=6)
        ttk.Button(frame_final, text="Run Step 4", style="Step.TButton",
                   command=lambda: self.start_thread(self.collect_final_wrapper)).grid(
            row=0, column=3, padx=(10, 6))

        # Clear All Button
        ttk.Button(outer, text="Clear All Fields",
                   command=self.clear_postproc_fields).pack(pady=(8, 0))

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
                    self.update_progress_ui(files_done_count, total_files,
                                            custom_text=f"Working on {display_count}/{total_files}: {sample_name}...")

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
