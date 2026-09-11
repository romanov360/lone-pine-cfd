"""
Assemble the standalone report page from its template.

The template in report/report.template.html holds the prose, the styling and
the chart code, with two placeholders:

  /*__DATA__*/     the JSON blob written by export_web.py
  __FIG_FIELDS__   the CFD temperature-field figure, as a data: URI

Keeping those out of the template means the page is regenerated from whatever
the latest run produced, and the template stays readable in a diff instead of
being a wall of base64.
"""

from __future__ import annotations

import base64
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
TEMPLATE = ROOT / "report" / "report.template.html"
DATA = ROOT / "results" / "data" / "web.json"
FIELDS = ROOT / "results" / "figures" / "09_fields.png"
OUT = ROOT / "report" / "report.html"


def data_uri(path: pathlib.Path) -> str:
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode()


def build(out: pathlib.Path = OUT) -> pathlib.Path:
    html = TEMPLATE.read_text()
    if not DATA.exists():
        raise SystemExit(f"missing {DATA} -- run src/export_web.py first")
    html = html.replace("/*__DATA__*/", DATA.read_text())
    if FIELDS.exists():
        html = html.replace("__FIG_FIELDS__", data_uri(FIELDS))
    else:
        print(f"warning: {FIELDS} missing, figure will be blank")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html)
    return out


if __name__ == "__main__":
    dest = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else OUT
    p = build(dest)
    print(f"wrote {p} ({p.stat().st_size/1024:.0f} kB)")
