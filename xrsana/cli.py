from __future__ import annotations

import cmd
import json
import os
import sys
import textwrap
import traceback
from pathlib import Path

# ── Infer project root before startup (xrsana/cli.py ≡ <root>/xrsana/cli.py) ──
_THIS_FILE = Path(__file__).resolve()
_PROJECT_ROOT = _THIS_FILE.parent.parent  # parent of xrsana/

# Prefer environment variables; otherwise infer from cli.py location
PROJECT_ROOT = Path(os.environ.get("XRSA_ROOT", _PROJECT_ROOT)).resolve()
DATA_ROOT    = Path(os.environ.get("XRSA_DATA", PROJECT_ROOT / "ex_space")).resolve()
RESOURCE_DIR = Path(os.environ.get("XRSA_RESOURCES", PROJECT_ROOT / "xrsana" / "resources" / "data")).resolve()

PLANNING_DIR = DATA_ROOT / "planning"
ANALYSIS_DIR = DATA_ROOT / "analysis"
PERFORM_DIR  = DATA_ROOT / "performing"  # fixed: perfoming -> performing

MPLCONFIGDIR = os.environ.get("MPLCONFIGDIR", "/tmp/xrsana-mpl")

# ── Lightweight welcome banner ───────────────────────────────────────────────
BANNER = r"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║   ██╗  ██╗██████╗ ███████╗ █████╗ ███╗   ██╗ █████╗                          ║
║   ╚██╗██╔╝██╔══██╗██╔════╝██╔══██╗████╗  ██║██╔══██╗                         ║
║    ╚███╔╝ ██████╔╝███████╗███████║██╔██╗ ██║███████║                         ║
║    ██╔██╗ ██╔══██╗╚════██║██╔══██║██║╚██╗██║██╔══██║                         ║
║   ██╔╝ ██╗██║  ██║███████║██║  ██║██║ ╚████║██║  ██║                         ║
║   ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝                         ║
║                                                                              ║
║              X-Ray Raman Scattering Interactive Toolkit                      ║
║                                                                              ║
║      XRStools-based platform for HEPS ID33 beamline analysis                 ║
║                                                                              ║
║               Planning  →  Performing  →  Analysis                           ║
║                     Integrated Scientific Workflow                           ║
║                                                                              ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  • Spectrometer Geometry & Analyzer Configuration                            ║
║  • Scan Planning / Batch Processing / Data Reduction                         ║
║  • Momentum Transfer & Energy Calibration                                    ║
║  • Interactive Visualization and Spectral Analysis                           ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""


def _print(msg: str = "", **kw) -> None:
    print(msg, **kw)


def _info(msg: str) -> None:
    _print(f"  [i] {msg}")


def _warn(msg: str) -> None:
    _print(f"  [!] {msg}")


def _ok(msg: str) -> None:
    _print(f"  [✓] {msg}")


def _hr() -> None:
    _print("-" * 62)


def _list_files(directory: Path, pattern: str = "*") -> list[Path]:
    if not directory.exists():
        return []
    return sorted(directory.glob(pattern))


# ══════════════════════════════════════════════════════════════════════════════
#  1) Pre-experiment Planning Shell
# ══════════════════════════════════════════════════════════════════════════════

class PlanningShell(cmd.Cmd):
    intro = textwrap.dedent("""\
        ┌────────────────────────────────────────────────────┐
        │  Planning Mode — XRS Spectrum Prediction           │
        │  list     list available input files (.inp)        │
        │  run <f>  run prediction (default: prediction.inp) │
        │  web      launch web parameter UI                  │
        │  save <o> save result to a text file               │
        │  back     return to main menu                      │
        └────────────────────────────────────────────────────┘
    """)
    prompt = "xrsana/planning> "

    def __init__(self) -> None:
        super().__init__()
        self._last_result = None
        self._last_filename: str | None = None

    # ── Command implementations ──────────────────────────────────────────

    def do_list(self, _arg: str) -> None:
        """List .inp files under ex_space/planning/input/"""
        inp_dir = PLANNING_DIR / "input"
        files = _list_files(inp_dir, "*.inp")
        if not files:
            _warn(f"No .inp files found (search path: {inp_dir})")
            return
        _info(f"Available input files ({len(files)}):")
        for f in files:
            _print(f"      {f.name}")

    def do_run(self, arg: str) -> None:
        """Run prediction: run [filename] (default: prediction.inp)"""
        from xrsana.xrs_prediction import (
            get_all_input, beam, sample, analyzer, detector,
            compton_profiles, thomson, absolute_cross_section
        )

        filename = arg.strip() or "prediction.inp"
        inp_path = (PLANNING_DIR / "input" / filename).resolve()

        if not inp_path.exists():
            _warn(f"File not found: {inp_path}")
            _info("Use 'list' to see available files")
            return

        try:
            _info(f"Reading input file: {inp_path}")
            inp = get_all_input(str(inp_path))

            # If analyzer database_dir is empty or relative (./ or .), point to resource dir automatically
            db_dir = inp.get("analyzer", {}).get("database_dir", "")
            if db_dir in ("", "./", "."):
                inp["analyzer"]["database_dir"] = str(RESOURCE_DIR.parent)
                _info(f"Auto-set analyzer database dir: {RESOURCE_DIR.parent}")

            _info("Building objects...")
            beam_obj = beam(
                inp["beam"]["i0_intensity"],
                inp["beam"]["beam_height"],
                inp["beam"]["beam_width"],
            )
            sample_obj = sample(
                inp["sample"]["chem_formulas"],
                inp["sample"]["concentrations"],
                inp["sample"]["densities"],
                inp["sample"]["angle_tth"],
                inp["sample"]["sample_thickness"],
                shape="sphere",
                molar_masses=inp["sample"]["molar_masses"],
            )
            analyzer_obj = analyzer(
                inp["analyzer"]["material"],
                inp["analyzer"]["hkl"],
                inp["analyzer"]["mask_d"],
                inp["analyzer"]["bend_r"],
                inp["analyzer"]["energy_resolution"],
                inp["analyzer"].get("diced", False),
                inp["analyzer"].get("thickness", 500.0),
                inp["analyzer"]["database_dir"],
            )
            detector_obj = detector(
                inp["detector"]["energy"],
                inp["detector"]["thickness"],
                inp["detector"]["material"],
                inp["detector"]["pixel_size"],
            )
            cp_obj = compton_profiles(
                sample_obj,
                inp["compton_profiles"]["eloss_range"],
                inp["compton_profiles"]["E0"],
            )
            thomson_obj = thomson(
                cp_obj.get_energy_in_keV(),
                cp_obj.get_E0(),
                cp_obj.get_tth(),
            )

            self._last_result = absolute_cross_section(
                beam_obj, sample_obj, analyzer_obj, detector_obj, thomson_obj, cp_obj
            )
            self._last_filename = filename

            self._last_result.plot_abs_cross_section()
            import matplotlib.pyplot as plt
            plt.show(block=False)
            _ok(f"Prediction finished (input: {filename})")
        except Exception as exc:
            _warn(f"Prediction failed: {exc}")
            traceback.print_exc()

    def do_web(self, _arg: str) -> None:
        """Launch prediction web UI (http://127.0.0.1:8765)"""
        from xrsana.xrs_prediction import serve_web

        _info("Starting prediction web UI...")
        _info("Open in browser: http://127.0.0.1:8765")
        _info("Press Ctrl-C to stop the server")
        try:
            serve_web("127.0.0.1", 8765)
        except KeyboardInterrupt:
            _print()

    def do_save(self, arg: str) -> None:
        """Save the latest prediction result: save [output_file]"""
        if self._last_result is None:
            _warn("Run 'run' first to generate a prediction")
            return
        output = arg.strip() or "prediction_result.txt"

        # Save to current working directory (explicit by design)
        out_path = Path(output).resolve()
        self._last_result.save_txt(
            str(out_path),
            header="energy_loss_eV absolute_counts_per_sec",
        )
        _ok(f"Saved to: {out_path}")

    def do_back(self, _arg: str) -> bool:
        """Return to main menu"""
        return True

    def do_quit(self, _arg: str) -> bool:
        return True

    # Aliases
    do_q = do_back
    do_exit = do_back
    do_ls = do_list

    def emptyline(self) -> None:
        pass  # Do not repeat the last command on empty input

    def default(self, line: str) -> None:
        _warn(f"Unknown command: {line}. Type '?' for help.")


# ══════════════════════════════════════════════════════════════════════════════
#  2) Experiment Performing Shell (placeholder)
# ══════════════════════════════════════════════════════════════════════════════

class PerformingShell(cmd.Cmd):
    intro = textwrap.dedent("""\
        ┌───────────────────────────────────────────────────┐
        │  Performing Mode — Acquisition & Management       │
        │  (Work in progress)                               │
        │                                                   │
        │  Planned: scan control / live monitor / archiving │
        │  back     return to main menu                     │
        └───────────────────────────────────────────────────┘
    """)
    prompt = "xrsana/performing> "

    def do_back(self, _arg: str) -> bool:
        return True

    do_q = do_back
    do_exit = do_back
    do_quit = do_back

    def emptyline(self) -> None:
        pass

    def default(self, line: str) -> None:
        if line.startswith("#"):
            return
        _info("Performing module is under construction. Type 'back' to return.")


# ══════════════════════════════════════════════════════════════════════════════
#  3) Data Analysis Shell
# ══════════════════════════════════════════════════════════════════════════════

class AnalysisShell(cmd.Cmd):
    intro = textwrap.dedent("""\
        ┌────────────────────────────────────────────────────────┐
        │  Analysis Mode — Background / Core-edge                │
        │                                                        │
        │  list                  list available datasets         │
        │  load <dir>            load a dataset                  │
        │  info                  show current dataset info       │
        │  elastic               remove elastic peak             │
        │  stray                 remove stray background         │
        │  polycore <q> [p] [c]  polynomial + core subtraction   │
        │  extract <q>           extract valence spectrum        │
        │  average <q1,...>      average multiple q-channels     │
        │  plot                  quick plot                      │
        │  save <path>           save results                    │
        │  run <script>          run an analysis script          │
        │  script <name>         create a script template        │
        │  back                  return to main menu             │
        └────────────────────────────────────────────────────────┘
    """)
    prompt = "xrsana/analysis> "

    def __init__(self) -> None:
        super().__init__()
        self._data = None
        self._processor = None
        self._result_dir: Path | None = None

    # ── Command implementations ──────────────────────────────────────────

    def do_list(self, _arg: str) -> None:
        """List available datasets under analysis/data/"""
        data_dir = ANALYSIS_DIR / "data"
        if not data_dir.exists():
            _warn(f"Data directory not found: {data_dir}")
            return
        dirs = sorted([d for d in data_dir.iterdir() if d.is_dir()])
        if not dirs:
            _warn("No datasets found")
            return
        _info("Available datasets:")
        for d in dirs:
            sub = sorted([s.name for s in d.iterdir() if s.is_dir()])
            _print(f"  {d.name}" + (f"  (scan: {', '.join(sub)})" if sub else ""))

    def do_info(self, _arg: str) -> None:
        """Show current loaded dataset info"""
        if self._data is None:
            _warn("No dataset loaded. Use load <dataset> first.")
            return
        _hr()
        _print(f"  Path         : {self._data.path}")
        _print(f"  Incident E0  : {self._data.E0} keV")
        _print(f"  Scattering tth: {self._data.tth}")
        _print(f"  q values     : {self._data.q}")
        _print(f"  Crystal/analyzer: {self._data.key}")
        _print(f"  Selected ROIs: {len(self._data.selected_rois)}")
        _print(f"  Points (eloss): {len(self._data.eloss)}")
        _hr()

    def do_load(self, arg: str) -> None:
        """Load dataset: load <dataset_name> [q_min q_max]"""
        from xrsana import xrs_read

        parts = arg.strip().split()
        if not parts:
            _warn("Usage: load <dataset_name> [q_min q_max]")
            _info("Use 'list' to see available datasets")
            return

        dataset = parts[0]
        q_range = None
        if len(parts) >= 3:
            try:
                q_range = (float(parts[1]), float(parts[2]))
            except ValueError:
                _warn("Invalid q-range format")
                return

        base = ANALYSIS_DIR / "data" / dataset
        if not base.exists():
            _warn(f"Dataset directory not found: {base}")
            _info("Use 'list' to see available datasets")
            return

        # Prefer a sub-folder if present (typical scan layout)
        sub_dirs = [d for d in base.iterdir() if d.is_dir()]
        if sub_dirs:
            data_dir = sub_dirs[0]
            _info(f"Using subdirectory: {data_dir.name}")
        else:
            data_dir = base

        try:
            _info(f"Loading: {data_dir}")
            self._data = xrs_read.read_heps_id33(
                str(data_dir),
                q_range=q_range,
            )
            self._processor = None

            import xrsana.xrs_process as xp
            elements = input("  Element symbols (e.g. Ho, comma-separated): ").strip()
            stoich = input("  Stoichiometry (comma-separated, e.g. 1): ").strip()
            edges_raw = input("  Core edges (e.g. Ho:N4, comma-separated): ").strip()

            if not elements:
                elements = "Ho"
            if not stoich:
                stoich = "1"
            if not edges_raw:
                edges_raw = "Ho:N4"

            elem_list = [e.strip() for e in elements.split(",")]
            stoich_list = [float(s.strip()) for s in stoich.split(",")]

            edges_dict: dict[str, list[str]] = {}
            for item in edges_raw.split(","):
                item = item.strip()
                if ":" in item:
                    k, v = item.split(":", 1)
                    edges_dict[k.strip()] = [e.strip() for e in v.split(";")]

            self._processor = xp.XRSProcess(self._data, elem_list, stoich_list, edges_dict)

            # Create result directory
            self._result_dir = ANALYSIS_DIR / "result" / dataset
            self._result_dir.mkdir(parents=True, exist_ok=True)

            _ok(f"Loaded: {dataset}")
            self.do_info("")
        except Exception as exc:
            _warn(f"Load failed: {exc}")
            traceback.print_exc()

    def do_elastic(self, _arg: str) -> None:
        """Remove elastic peak"""
        if self._processor is None:
            _warn("Load data first")
            return
        _info("Removing elastic peak...")
        try:
            self._processor.xrs_remove_elastic()
            _ok("Elastic peak removed")
        except Exception as exc:
            _warn(f"Failed: {exc}")
            traceback.print_exc()

    def do_stray(self, _arg: str) -> None:
        """Remove stray background"""
        if self._processor is None:
            _warn("Load data first")
            return
        _info("Removing stray background...")
        try:
            self._processor.xrs_remove_stray_background()
            _ok("Stray background removed")
        except Exception as exc:
            _warn(f"Failed: {exc}")
            traceback.print_exc()

    def do_polycore(self, arg: str) -> None:
        """Polynomial + core subtraction: polycore <q_index> [poly_low poly_high] [core_low core_high]"""
        if self._processor is None:
            _warn("Load data first")
            return
        parts = arg.strip().split()
        if not parts:
            _warn("Usage: polycore <q_index> [poly_low poly_high] [core_low core_high]")
            _info("Example: polycore 1 200 600 150 170")
            return
        try:
            whichq = int(parts[0])
            polyregion = (float(parts[1]), float(parts[2])) if len(parts) >= 3 else (200, 600)
            coreregion = (float(parts[3]), float(parts[4])) if len(parts) >= 5 else (150, 170)
            _info(f"Running polycore: q={whichq}, poly={polyregion}, core={coreregion}")
            self._processor.xrs_remove_poly_core(
                whichq=whichq,
                polyregion=polyregion,
                coreregion=coreregion,
                polyorder=2,
                plot=True,
            )
            _ok("Polynomial + core subtraction done")
        except Exception as exc:
            _warn(f"Failed: {exc}")
            traceback.print_exc()

    def do_extract(self, arg: str) -> None:
        """Extract valence spectrum: extract <q_index>"""
        if self._processor is None:
            _warn("Load data first")
            return
        parts = arg.strip().split()
        if not parts:
            _warn("Usage: extract <q_index>")
            return
        try:
            whichq = int(parts[0])
            _info(f"Extracting valence spectrum: q={whichq}")
            self._processor.extractval(whichq=whichq, make_plots=False)
            _ok("Valence spectrum extracted")
        except Exception as exc:
            _warn(f"Failed: {exc}")
            traceback.print_exc()

    def do_average(self, arg: str) -> None:
        """Average multiple q-channels: average <q1,q2,q3,...>"""
        if self._processor is None:
            _warn("Load data first")
            return
        parts = arg.strip()
        if not parts:
            _warn("Usage: average <q1,q2,q3,...>")
            _info("Example: average 0,2,3,4,5")
            return
        try:
            whichq = [int(x.strip()) for x in parts.split(",")]
            _info(f"Averaging q-channels: {whichq}")
            result = self._processor.averageqs(whichq=whichq, return_result=True)
            _ok(f"Done. Result shape: {result.shape if hasattr(result, 'shape') else 'ok'}")
        except Exception as exc:
            _warn(f"Failed: {exc}")
            traceback.print_exc()

    def do_plot(self, _arg: str) -> None:
        """Quick plot of current signals"""
        if self._data is None and self._processor is None:
            _warn("Load data first")
            return
        import matplotlib.pyplot as plt

        data_src = self._processor if self._processor else self._data
        eloss = data_src.eloss
        signals = data_src.signals

        fig, ax = plt.subplots(figsize=(8, 5))
        if signals.ndim == 1:
            ax.plot(eloss, signals)
        else:
            for i in range(signals.shape[1]):
                ax.plot(eloss, signals[:, i], label=f"q-channel {i}")
            ax.legend(fontsize=8)
        ax.set_xlabel("Energy Loss (eV)")
        ax.set_ylabel("Intensity")
        ax.set_title("Quick Plot")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        plt.show(block=False)
        _ok("Plot window opened")

    def do_run(self, arg: str) -> None:
        """Run analysis script: run <script.py>"""
        filepath = arg.strip()
        if not filepath:
            _warn("Usage: run <script.py>")
            return
        script_path = Path(filepath).resolve()
        if not script_path.exists():
            script_path = ANALYSIS_DIR / "scripts" / filepath
        if not script_path.exists():
            _warn(f"Script not found: {script_path}")
            return
        import runpy
        _info(f"Running script: {script_path}")
        try:
            runpy.run_path(str(script_path))
            _ok("Script finished")
        except Exception as exc:
            _warn(f"Script failed: {exc}")
            traceback.print_exc()

    def do_script(self, arg: str) -> None:
        """Create a script template under analysis/scripts/: script <name>"""
        name = arg.strip()
        if not name:
            _warn("Usage: script <name>")
            return
        if not name.endswith(".py"):
            name += ".py"
        scripts_dir = ANALYSIS_DIR / "scripts"
        scripts_dir.mkdir(parents=True, exist_ok=True)
        out_path = scripts_dir / name
        if out_path.exists():
            _warn(f"File already exists: {out_path}")
            return
        content = textwrap.dedent(f"""\
            # XRSana analysis script: {name}
            import matplotlib.pyplot as plt
            from xrsana import xrs_read
            from xrsana import xrs_process

            DATA_PATH = "{ANALYSIS_DIR}/data/Ho"
            data = xrs_read.read_heps_id33(DATA_PATH, q_range=(9.2, 10))
            proc = xrs_process.XRSProcess(data, ["Ho"], [1], {{"Ho": ["N4"]}})

            proc.xrs_remove_elastic()
            proc.xrs_remove_stray_background()
            # proc.xrs_remove_poly_core(whichq=1, polyregion=(200, 600), coreregion=(150, 170))
            # proc.extractval(whichq=1)
            # result = proc.averageqs(whichq=[0, 2, 3])

            plt.plot(proc.eloss, proc.signals)
            plt.xlabel("Energy Loss (eV)")
            plt.ylabel("Intensity")
            plt.title("Quick Plot")
            plt.grid(True, alpha=0.3)
            plt.show()
        """)
        out_path.write_text(content)
        _ok(f"Template created: {out_path}")

    def do_save(self, arg: str) -> None:
        """Save results: save <output_path>"""
        if self._processor is None:
            _warn("Load data and run processing first")
            return
        out = arg.strip()
        if not out:
            if self._result_dir:
                out = str(self._result_dir / "result.dat")
            else:
                out = "result.dat"
        import numpy as np
        out_path = Path(out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            if hasattr(self._processor, "valence") and self._processor.valence is not None:
                np.savetxt(
                    str(out_path),
                    np.column_stack([self._processor.eloss, self._processor.valence]),
                    header="eloss_eV  valence",
                )
                _ok(f"Saved: {out_path}")
            else:
                _warn("No valence spectrum to save. Run 'extract' first.")
        except Exception as exc:
            _warn(f"Save failed: {exc}")

    def do_back(self, _arg: str) -> bool:
        """Return to main menu"""
        return True

    def do_quit(self, _arg: str) -> bool:
        return True

    do_q = do_back
    do_exit = do_back
    do_ls = do_list

    def emptyline(self) -> None:
        pass

    def default(self, line: str) -> None:
        if line.startswith("#"):
            return
        _warn(f"Unknown command: {line}. Type '?' for help.")


# ══════════════════════════════════════════════════════════════════════════════
#  Root Shell (Main Menu)
# ══════════════════════════════════════════════════════════════════════════════

class XRSanaShell(cmd.Cmd):
    intro = BANNER + textwrap.dedent("""\
           Workspace:
             PROJECT_ROOT  = %s
             DATA_ROOT     = %s
             RESOURCE_DIR  = %s

           Enter a number or a command:
             planning   (or 1)  enter Planning mode
             performing (or 2)  enter Performing mode
             analysis   (or 3)  enter Analysis mode
             run <f>            run a Python script (macro)
             help               show help
             quit               exit
    """ % (PROJECT_ROOT, DATA_ROOT, RESOURCE_DIR))

    prompt = "xrsana> "

    def __init__(self) -> None:
        super().__init__()

        _info("Loading modules...")
        self._loaded = _load_modules()

        loaded_info = textwrap.indent(
            "\n".join(self._loaded["lines"]),
            "  "
        )

        self.intro += (
            "\n\n"
            "Preloaded modules:\n"
            f"{loaded_info}"
        )

    # ── Main menu commands ──────────────────────────────────────────────────

    def do_planning(self, _arg: str) -> None:
        """Enter Planning mode (spectrum prediction)"""
        _print()
        PlanningShell().cmdloop()
        _print()

    def do_performing(self, _arg: str) -> None:
        """Enter Performing mode"""
        _print()
        PerformingShell().cmdloop()
        _print()

    def do_analysis(self, _arg: str) -> None:
        """Enter Analysis mode"""
        _print()
        AnalysisShell().cmdloop()
        _print()

    def do_run(self, arg: str) -> None:
        """Run a Python script: run <script.py> (or run --new <name> to create a template)"""
        filepath = arg.strip()
        if not filepath:
            _info("Usage: run <script.py>")
            _info("       run --new <name>  create a script template")
            return

        if filepath.startswith("--new"):
            name = filepath[5:].strip()
            if not name:
                _warn("Usage: run --new <script_name>")
                return
            if not name.endswith(".py"):
                name += ".py"
            scripts_dir = ANALYSIS_DIR / "scripts"
            scripts_dir.mkdir(parents=True, exist_ok=True)
            out_path = scripts_dir / name
            if out_path.exists():
                _warn(f"File already exists: {out_path}")
                return
            content = textwrap.dedent(f"""\
                # XRSana script: {name}
                import numpy as np
                import matplotlib.pyplot as plt
                from xrsana import xrs_read
                from xrsana import xrs_public
                from xrsana import xrs_process

                # Write your analysis code here
                print("Hello from {name}")
            """)
            out_path.write_text(content)
            _ok(f"Template created: {out_path}")
            return

        script_path = Path(filepath).resolve()
        if not script_path.exists():
            script_path = ANALYSIS_DIR / "scripts" / filepath
        if not script_path.exists():
            _warn(f"Script not found: {script_path}")
            return

        import runpy
        _info(f"Running script: {script_path}")
        try:
            runpy.run_path(str(script_path))
            _ok("Script finished")
        except Exception as exc:
            _warn(f"Script failed: {exc}")
            traceback.print_exc()

    def do_quit(self, _arg: str) -> bool:
        """Exit XRSana"""
        _print("\nBye!\n")
        return True

    def do_help(self, arg: str) -> None:
        if arg:
            super().do_help(arg)
        else:
            super().do_help("")

    # Numeric aliases: 1->planning, 2->performing, 3->analysis
    def do_1(self, _arg: str) -> None:
        self.do_planning("")

    def do_2(self, _arg: str) -> None:
        self.do_performing("")

    def do_3(self, _arg: str) -> None:
        self.do_analysis("")

    do_q = do_quit
    do_exit = do_quit

    def emptyline(self) -> None:
        pass

    def default(self, line: str) -> None:
        if line.startswith("#") or line.strip() == "":
            return
        _warn(f"Unknown command: {line}. Type 'help' to see commands.")


# ══════════════════════════════════════════════════════════════════════════════
#  Module loading
# ══════════════════════════════════════════════════════════════════════════════

_MODULE_TABLE = [
    ("xrs_read", "Data reader"),
    ("xrs_public", "Physics utilities"),
    ("xrs_prediction", "Experiment prediction"),
    ("xrs_process", "Analysis pipeline"),
    ("xrs_ComptonProfiles", "HF Compton profiles"),
    ("xrs_extraction", "Legacy extraction pipeline"),
    ("math_functions", "Fitting functions"),
    ("data_browser", "Web data browser"),
]


def _load_modules() -> dict:
    loaded = []
    for mod_name, desc in _MODULE_TABLE:
        try:
            __import__(f"xrsana.{mod_name}", fromlist=[mod_name])
            loaded.append(f"  [✓] {mod_name:<25s} {desc}")
        except ImportError:
            loaded.append(f"  [✗] {mod_name:<25s} {desc}  (unavailable)")
    return {"lines": loaded}


# ══════════════════════════════════════════════════════════════════════════════
#  Entry point
# ══════════════════════════════════════════════════════════════════════════════

def main(argv: list[str] | None = None) -> None:
    """CLI entry point."""

    import argparse
    import sys as _sys

    _sys.path.insert(0, str(PROJECT_ROOT))

    parser = argparse.ArgumentParser(
        prog="xrsana",
        description="XRSana -- Interactive XRS Analysis Toolkit",
    )
    parser.add_argument(
        "script",
        nargs="?",
        default=None,
        help="Run a Python script directly (skip interactive UI)",
    )
    parser.add_argument(
        "--no-interactive",
        action="store_true",
        help="Do not enter interactive mode after loading (use with script)",
    )
    parser.add_argument(
        "--browse",
        action="store_true",
        help="Launch data browser (same as xrsana-browse-data)",
    )
    args, _unknown = parser.parse_known_args(argv)

    # Configure matplotlib
    os.environ.setdefault("MPLCONFIGDIR", MPLCONFIGDIR)
    import matplotlib
    matplotlib.use(os.environ.get("MPLBACKEND", "TkAgg"))

    # Browser mode
    if args.browse:
        from xrsana.data_browser import main as browse_main
        browse_main()
        return

    # Script mode: run and exit
    if args.script:
        import runpy
        script_path = Path(args.script).resolve()
        if not script_path.exists():
            _warn(f"Script not found: {script_path}")
            _sys.exit(1)
        _info(f"Running script: {script_path}")
        runpy.run_path(str(script_path))
        return

    # Interactive mode
    try:
        XRSanaShell().cmdloop()
    except KeyboardInterrupt:
        _print("\n\nBye!\n")
    except EOFError:
        _print("\nBye!\n")


if __name__ == "__main__":
    main()