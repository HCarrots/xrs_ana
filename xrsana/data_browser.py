"""Small web browser for reduced XRS crystal spectra."""

from __future__ import annotations

import argparse
import io
import json
import math
import os
import threading
import webbrowser
from dataclasses import dataclass
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Dict, Iterable, List, Optional
from urllib.parse import parse_qs, unquote, urlparse

os.environ.setdefault("MPLCONFIGDIR", os.path.join("/tmp", "xrsana-matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_DATA_DIR = "/home/hushiqi/work/xrs_ana/ex_space/analysis/data/Ho"
ENERGY_COLUMN = "Energy Transfer (eV)"
HC_OVER_2PI = 0.5067731


@dataclass
class CrystalRecord:
    index: int
    crystal: str
    energy: np.ndarray
    intensity: np.ndarray
    q_ave: float
    scattering_angle: float
    center: float
    width: float
    r_square: float

    def to_dict(self) -> Dict[str, object]:
        return {
            "index": self.index,
            "crystal": self.crystal,
            "q_ave": _finite_or_none(self.q_ave),
            "scattering_angle": _finite_or_none(self.scattering_angle),
            "center": _finite_or_none(self.center),
            "width": _finite_or_none(self.width),
            "r_square": _finite_or_none(self.r_square),
            "points": int(len(self.energy)),
            "intensity_min": _finite_or_none(float(np.nanmin(self.intensity))),
            "intensity_max": _finite_or_none(float(np.nanmax(self.intensity))),
        }


class ReducedDataBrowser:
    """Load a reduced text export and render one labeled image per crystal."""

    def __init__(
        self,
        data_dir: str = DEFAULT_DATA_DIR,
        scan_name: Optional[str] = None,
        data_file: Optional[str] = None,
        roi_file: Optional[str] = None,
    ) -> None:
        self.data_dir = os.path.abspath(data_dir)
        self.scan_name = scan_name
        self.data_file = os.path.abspath(data_file) if data_file else None
        self.roi_file = os.path.abspath(roi_file) if roi_file else None
        self.records: List[CrystalRecord] = []
        self.records_by_crystal: Dict[str, CrystalRecord] = {}
        self.load()

    def load(self) -> None:
        files = self._resolve_files()
        dataframe = pd.read_csv(files["data_file"], sep="\t")
        rois = pd.read_csv(files["roi_file"], sep="\t")
        dataframe = _drop_unnamed_columns(dataframe)
        rois = _drop_unnamed_columns(rois)

        if ENERGY_COLUMN in dataframe.columns:
            energy_column = ENERGY_COLUMN
        else:
            energy_column = dataframe.columns[0]

        roi_by_crystal = rois.set_index("crystal", drop=False) if "crystal" in rois else pd.DataFrame()
        crystals = [
            column
            for column in dataframe.columns
            if column != energy_column and not str(column).startswith("Unnamed:")
        ]
        energy = dataframe[energy_column].to_numpy(dtype="float64")

        self.records = []
        for index, crystal in enumerate(crystals, start=1):
            roi = roi_by_crystal.loc[crystal] if crystal in roi_by_crystal.index else {}
            q_ave = _row_float(roi, "q_ave")
            center = _row_float(roi, "center")
            record = CrystalRecord(
                index=index,
                crystal=str(crystal),
                energy=energy,
                intensity=dataframe[crystal].to_numpy(dtype="float64"),
                q_ave=q_ave,
                scattering_angle=scattering_angle_degrees(q_ave, center),
                center=center,
                width=_row_float(roi, "width"),
                r_square=_row_float(roi, "r-square"),
            )
            self.records.append(record)
        self.records_by_crystal = {record.crystal: record for record in self.records}

    def _resolve_files(self) -> Dict[str, str]:
        if self.data_file:
            data_file = self.data_file
            scan_root = _scan_root(data_file)
            roi_file = self.roi_file or os.path.join(os.path.dirname(data_file), scan_root + "_rois.txt")
            return {"data_file": data_file, "roi_file": roi_file}

        if not os.path.isdir(self.data_dir):
            raise FileNotFoundError(f"Data directory does not exist: {self.data_dir}")

        candidates = []
        for filename in sorted(os.listdir(self.data_dir)):
            if filename.endswith("_all_data.txt") or filename.endswith("_data.txt"):
                if self.scan_name and not filename.startswith(self.scan_name):
                    continue
                data_file = os.path.join(self.data_dir, filename)
                scan_root = _scan_root(data_file)
                roi_file = os.path.join(self.data_dir, scan_root + "_rois.txt")
                if os.path.isfile(roi_file):
                    candidates.append({"data_file": data_file, "roi_file": roi_file})
        if not candidates:
            raise FileNotFoundError(
                f"No reduced data triplet found in {self.data_dir}; expected *_data.txt or *_all_data.txt with *_rois.txt"
            )
        return candidates[0]

    def filtered_records(
        self,
        query: str = "",
        q_min: Optional[float] = None,
        q_max: Optional[float] = None,
        sort: str = "index",
    ) -> List[CrystalRecord]:
        records: Iterable[CrystalRecord] = self.records
        if query:
            needle = query.lower()
            records = [record for record in records if needle in record.crystal.lower()]
        if q_min is not None:
            records = [record for record in records if math.isnan(record.q_ave) or record.q_ave >= q_min]
        if q_max is not None:
            records = [record for record in records if math.isnan(record.q_ave) or record.q_ave <= q_max]

        sorters = {
            "index": lambda record: record.index,
            "crystal": lambda record: record.crystal,
            "q": lambda record: (math.inf if math.isnan(record.q_ave) else record.q_ave),
            "angle": lambda record: (math.inf if math.isnan(record.scattering_angle) else record.scattering_angle),
        }
        return sorted(records, key=sorters.get(sort, sorters["index"]))

    def render_png(self, crystal: str) -> bytes:
        record = self.records_by_crystal[crystal]
        fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=150)
        ax.plot(record.energy, record.intensity, color="#165a72", linewidth=1.35)
        ax.set_xlabel("Energy transfer (eV)")
        ax.set_ylabel("Intensity")
        ax.grid(True, color="#d7dee3", linewidth=0.6, alpha=0.8)
        ax.set_title(
            (
                f"{record.crystal}  |  Crystal #{record.index}  |  "
                f"q={_format_number(record.q_ave, 3)} A^-1  |  "
                f"2theta={_format_number(record.scattering_angle, 2)} deg"
            ),
            fontsize=10,
        )
        ax.text(
            0.015,
            0.96,
            f"E0={_format_number(record.center, 2)} eV\nFWHM={_format_number(record.width, 3)} eV",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#c8d1d8", "alpha": 0.9},
        )
        fig.tight_layout()
        buffer = io.BytesIO()
        fig.savefig(buffer, format="png")
        plt.close(fig)
        return buffer.getvalue()


def scattering_angle_degrees(q_ave: float, e0_eV: float) -> float:
    """Return 2theta in degrees from q magnitude and incident energy."""

    if math.isnan(q_ave) or math.isnan(e0_eV) or e0_eV <= 0.0:
        return float("nan")
    k = HC_OVER_2PI * (e0_eV / 1000.0)
    cos_tth = 1.0 - (q_ave**2) / (2.0 * k**2)
    return float(np.degrees(np.arccos(np.clip(cos_tth, -1.0, 1.0))))


def run_browser(
    data_dir: str = DEFAULT_DATA_DIR,
    host: str = "127.0.0.1",
    port: int = 8765,
    open_browser: bool = False,
    scan_name: Optional[str] = None,
    data_file: Optional[str] = None,
    roi_file: Optional[str] = None,
) -> ThreadingHTTPServer:
    browser = ReducedDataBrowser(
        data_dir=data_dir,
        scan_name=scan_name,
        data_file=data_file,
        roi_file=roi_file,
    )
    handler = _make_handler(browser)
    server = ThreadingHTTPServer((host, port), handler)
    url = f"http://{host}:{port}/"
    print(f"xrsana data browser serving {len(browser.records)} crystals at {url}")
    if open_browser:
        threading.Timer(0.2, lambda: webbrowser.open(url)).start()
    server.serve_forever()
    return server


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description="Browse reduced XRS crystal spectra in a local web UI.")
    parser.add_argument("data_dir", nargs="?", default=DEFAULT_DATA_DIR, help="Directory containing *_data.txt and *_rois.txt")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--scan-name", default=None)
    parser.add_argument("--data-file", default=None)
    parser.add_argument("--roi-file", default=None)
    parser.add_argument("--open", action="store_true", help="Open the browser automatically")
    args = parser.parse_args(argv)
    run_browser(
        data_dir=args.data_dir,
        host=args.host,
        port=args.port,
        open_browser=args.open,
        scan_name=args.scan_name,
        data_file=args.data_file,
        roi_file=args.roi_file,
    )


def _make_handler(browser: ReducedDataBrowser):
    class DataBrowserHandler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            parsed = urlparse(self.path)
            if parsed.path == "/":
                self._send_text(INDEX_HTML, "text/html; charset=utf-8")
                return
            if parsed.path == "/api/records":
                params = parse_qs(parsed.query)
                records = browser.filtered_records(
                    query=_first(params, "q", ""),
                    q_min=_optional_float(_first(params, "q_min", "")),
                    q_max=_optional_float(_first(params, "q_max", "")),
                    sort=_first(params, "sort", "index"),
                )
                self._send_json(
                    {
                        "count": len(records),
                        "total": len(browser.records),
                        "records": [record.to_dict() for record in records],
                    }
                )
                return
            if parsed.path.startswith("/image/"):
                crystal = unquote(parsed.path[len("/image/") :])
                if crystal not in browser.records_by_crystal:
                    self.send_error(HTTPStatus.NOT_FOUND, "Unknown crystal")
                    return
                self._send_bytes(browser.render_png(crystal), "image/png")
                return
            self.send_error(HTTPStatus.NOT_FOUND)

        def log_message(self, format: str, *args: object) -> None:
            return

        def _send_json(self, payload: Dict[str, object]) -> None:
            self._send_text(json.dumps(payload), "application/json; charset=utf-8")

        def _send_text(self, payload: str, content_type: str) -> None:
            self._send_bytes(payload.encode("utf-8"), content_type)

        def _send_bytes(self, payload: bytes, content_type: str) -> None:
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(payload)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(payload)

    return DataBrowserHandler


def _scan_root(path: str) -> str:
    basename = os.path.basename(path)
    if basename.endswith("_all_data.txt"):
        return basename[: -len("_all_data.txt")]
    if basename.endswith("_data.txt"):
        return basename[: -len("_data.txt")]
    return os.path.splitext(basename)[0]


def _drop_unnamed_columns(dataframe: pd.DataFrame) -> pd.DataFrame:
    return dataframe.loc[:, [column for column in dataframe.columns if not str(column).startswith("Unnamed:")]]


def _row_float(row: object, key: str) -> float:
    try:
        value = row[key]  # type: ignore[index]
    except Exception:
        return float("nan")
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _finite_or_none(value: float) -> Optional[float]:
    return value if math.isfinite(value) else None


def _optional_float(value: str) -> Optional[float]:
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _first(params: Dict[str, List[str]], key: str, default: str) -> str:
    values = params.get(key)
    return values[0] if values else default


def _format_number(value: float, digits: int) -> str:
    if not math.isfinite(value):
        return "n/a"
    return f"{value:.{digits}f}"


INDEX_HTML = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>XRS Data Browser</title>
  <style>
    :root {
      color-scheme: light;
      --bg: #f6f8f9;
      --surface: #ffffff;
      --ink: #172126;
      --muted: #617078;
      --line: #d8e0e5;
      --accent: #0f6f78;
      --accent-soft: #e2f0f1;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: var(--bg);
      color: var(--ink);
    }
    header {
      position: sticky;
      top: 0;
      z-index: 5;
      display: grid;
      grid-template-columns: minmax(180px, 1fr) auto;
      gap: 16px;
      align-items: end;
      padding: 16px 22px;
      background: rgba(246, 248, 249, 0.96);
      border-bottom: 1px solid var(--line);
      backdrop-filter: blur(10px);
    }
    h1 {
      margin: 0;
      font-size: 20px;
      font-weight: 650;
      letter-spacing: 0;
    }
    .summary {
      margin-top: 4px;
      color: var(--muted);
      font-size: 13px;
    }
    .controls {
      display: grid;
      grid-template-columns: 170px 90px 90px 130px;
      gap: 8px;
      align-items: center;
    }
    input, select {
      width: 100%;
      min-height: 36px;
      border: 1px solid var(--line);
      background: var(--surface);
      color: var(--ink);
      border-radius: 6px;
      padding: 7px 9px;
      font: inherit;
      font-size: 13px;
    }
    main { padding: 18px 22px 28px; }
    .grid {
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(340px, 1fr));
      gap: 14px;
    }
    .card {
      overflow: hidden;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: var(--surface);
      box-shadow: 0 1px 2px rgba(18, 30, 35, 0.04);
    }
    .meta {
      display: grid;
      grid-template-columns: 1fr auto;
      gap: 8px;
      align-items: start;
      padding: 10px 12px 8px;
      border-bottom: 1px solid var(--line);
    }
    .crystal {
      font-size: 15px;
      font-weight: 650;
      min-width: 0;
    }
    .badge {
      justify-self: end;
      padding: 3px 7px;
      border-radius: 999px;
      background: var(--accent-soft);
      color: var(--accent);
      font-size: 12px;
      font-weight: 650;
    }
    .metrics {
      grid-column: 1 / -1;
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 6px;
      color: var(--muted);
      font-size: 12px;
    }
    .metric strong {
      display: block;
      color: var(--ink);
      font-size: 13px;
      font-weight: 600;
      margin-top: 2px;
    }
    img {
      display: block;
      width: 100%;
      aspect-ratio: 64 / 42;
      object-fit: contain;
      background: #fff;
    }
    .empty {
      padding: 26px;
      border: 1px dashed var(--line);
      border-radius: 8px;
      color: var(--muted);
      background: var(--surface);
    }
    @media (max-width: 760px) {
      header { grid-template-columns: 1fr; align-items: stretch; }
      .controls { grid-template-columns: 1fr 1fr; }
      main { padding: 14px; }
      .grid { grid-template-columns: 1fr; }
    }
  </style>
</head>
<body>
  <header>
    <div>
      <h1>XRS Data Browser</h1>
      <div class="summary" id="summary">Loading spectra...</div>
    </div>
    <div class="controls">
      <input id="query" type="search" placeholder="Crystal">
      <input id="qMin" type="number" step="0.1" placeholder="q min">
      <input id="qMax" type="number" step="0.1" placeholder="q max">
      <select id="sort">
        <option value="index">Crystal #</option>
        <option value="crystal">Crystal name</option>
        <option value="q">q magnitude</option>
        <option value="angle">Scattering angle</option>
      </select>
    </div>
  </header>
  <main>
    <div id="grid" class="grid"></div>
  </main>
  <script>
    const grid = document.getElementById('grid');
    const summary = document.getElementById('summary');
    const controls = ['query', 'qMin', 'qMax', 'sort'].map(id => document.getElementById(id));
    let timer = null;

    controls.forEach(control => control.addEventListener('input', () => {
      clearTimeout(timer);
      timer = setTimeout(loadRecords, 120);
    }));

    function value(id) {
      return document.getElementById(id).value.trim();
    }

    function fmt(value, digits) {
      return value === null || Number.isNaN(value) ? 'n/a' : Number(value).toFixed(digits);
    }

    async function loadRecords() {
      const params = new URLSearchParams();
      if (value('query')) params.set('q', value('query'));
      if (value('qMin')) params.set('q_min', value('qMin'));
      if (value('qMax')) params.set('q_max', value('qMax'));
      params.set('sort', value('sort') || 'index');
      const response = await fetch('/api/records?' + params.toString());
      const payload = await response.json();
      summary.textContent = `Showing ${payload.count} of ${payload.total} crystal spectra`;
      grid.innerHTML = '';
      if (!payload.records.length) {
        grid.innerHTML = '<div class="empty">No spectra match the current filters.</div>';
        return;
      }
      for (const record of payload.records) {
        const card = document.createElement('article');
        card.className = 'card';
        card.innerHTML = `
          <div class="meta">
            <div class="crystal">${record.crystal}</div>
            <div class="badge">#${record.index}</div>
            <div class="metrics">
              <div class="metric">q |A^-1|<strong>${fmt(record.q_ave, 3)}</strong></div>
              <div class="metric">2theta deg<strong>${fmt(record.scattering_angle, 2)}</strong></div>
              <div class="metric">R^2<strong>${fmt(record.r_square, 3)}</strong></div>
            </div>
          </div>
          <img loading="lazy" alt="${record.crystal} spectrum labeled with q, scattering angle, and crystal number" src="/image/${encodeURIComponent(record.crystal)}">
        `;
        grid.appendChild(card);
      }
    }

    loadRecords();
  </script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
