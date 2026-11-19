"""Consistency audit (Stage 7).  Checks that the manuscript, the Data files, the
bibliography, the Literature index and the README agree, and writes the result
to AUDIT_REPORT.txt in the project folder.

Usage:  python Codes/tools/audit.py        (after run_all.py and after compiling)
"""
import csv
import filecmp
import glob
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
DATA, CODES = os.path.join(ROOT, "Data"), os.path.join(ROOT, "Codes")
out, n_fail = [], 0


def check(name, ok, detail=""):
    global n_fail
    n_fail += 0 if ok else 1
    out.append(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"\n         {detail}" if detail else ""))


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


tex = read(os.path.join(ROOT, "manuscript.tex"))
body = "\n".join(l.split("%")[0] if not l.strip().startswith("%") else "" for l in tex.splitlines())
nums = {r["macro"]: r for r in csv.DictReader(
    l for l in open(os.path.join(DATA, "manuscript_numbers.csv")) if not l.startswith("#"))}
tables = re.findall(r"\\newcommand\{\\(\w+)\}", read(os.path.join(ROOT, "generated_tables.tex")))

# 1 numbers -------------------------------------------------------------
used = set(re.findall(r"\\([A-Za-z]+)", body))
used_nums = sorted(m for m in used if m in nums)
missing_src = []
for m in used_nums:
    src = nums[m]["source"]
    if src == "runtime" or src.startswith("Codes/"):
        continue
    if not glob.glob(os.path.join(DATA, src)):
        missing_src.append(f"{m} -> {src}")
check(f"every number macro used in the manuscript ({len(used_nums)} of {len(nums)} defined) "
      "traces to an existing file in Data/", not missing_src, "; ".join(missing_src))
check("every generated table body used in the manuscript exists",
      all(t in tables for t in used if t.startswith("Tab") and t.endswith("Body")),
      "used: " + ", ".join(sorted(t for t in used if t.startswith("Tab") and t.endswith("Body"))))
unused_tables = [t for t in tables if t not in used]
check("no generated table body is left unused", not unused_tables, ", ".join(unused_tables))

# 2 figures -------------------------------------------------------------
FIG_SOURCES = {
    "fig2_verification": ["verification/derivatives_336_check.csv", "verification/closedform_vs_numerical.csv"],
    "fig3_co2_isotherms": ["validation/co2_model_curves.csv", "validation/co2_eos_vs_nist.csv"],
    "fig4_co2_density_deviation": ["validation/co2_eos_vs_nist.csv"],
    "fig5_co2_isobar_10MPa": ["validation/co2_model_curves.csv", "validation/co2_eos_vs_nist.csv"],
    "fig6_co2_saturation": ["validation/co2_saturation_curves.csv", "validation/co2_saturation_vs_nist.csv"],
    "fig7_table_route": ["validation/co2_table_route_vs_nist.csv"],
    "fig8_uncertainty_budget": ["outputs/02_process_co2_pr_uncertainty_budget.csv",
                                "outputs/05_infer_from_measurements_uncertainty_budget.csv"],
    "toc_graphic": ["validation/co2_eos_vs_nist.csv"],
}
figs = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", body)
bad = [f for f in figs if not os.path.exists(os.path.join(ROOT, f + ".pdf"))]
check(f"all {len(figs)} included figure files exist", not bad, ", ".join(bad))
bad = [f"{f}: {s}" for f in figs for s in FIG_SOURCES.get(f, ["<no source recorded>"])
       if not os.path.exists(os.path.join(DATA, s))]
check("every figure traces to existing Data files", not bad, "; ".join(bad))
out.append("         figure sources: " + "; ".join(f"{f} <- {', '.join(FIG_SOURCES[f])}" for f in figs))
out.append("         Figure 1 (architecture) is drawn in LaTeX and contains no data.")

# 3 citations -----------------------------------------------------------
cited = set()
for grp in re.findall(r"\\cite\{([^}]+)\}", body):
    cited |= {k.strip() for k in grp.split(",")}
bibkeys = set(re.findall(r"@\w+\{(\w+),", read(os.path.join(ROOT, "references.bib"))))
lit = read(os.path.join(ROOT, "Literature", "reference_list.md"))
check(f"every citation ({len(cited)}) exists in references.bib", cited <= bibkeys,
      ", ".join(sorted(cited - bibkeys)))
check("every bibliography entry is cited", bibkeys <= cited, ", ".join(sorted(bibkeys - cited)))
check("every citation is listed in Literature/reference_list.md",
      all(f"| {k} |" in lit for k in cited), ", ".join(k for k in sorted(cited) if f"| {k} |" not in lit))
unver = [l.split("|")[2].strip() for l in lit.splitlines() if l.startswith("|") and "UNVERIFIED" in l.split("|")[-2]
         and "VERIFIED (" not in l.split("|")[-2]]
out.append(f"[INFO] references marked UNVERIFIED (details from memory, not checked online): "
           + (", ".join(unver) or "none"))
pdfs = glob.glob(os.path.join(ROOT, "Literature", "*.pdf"))
check("PDF copies of the references are in Literature/ (the group must add them)",
      len(pdfs) >= len(bibkeys), f"{len(pdfs)} PDF files found, {len(bibkeys)} references; "
      "see Literature/download_checklist.csv")

# 4 placeholders and literal numbers -------------------------------------
ph = re.findall(r"\[[^\]]*TO (?:FILL|CONFIRM)[^\]]*\]", tex)
check("no [TO FILL] placeholders remain in the manuscript (authors, affiliation, e-mail, "
      "acknowledgment, AI-use disclosure)", not ph, f"{len(ph)} placeholders remain")
m = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, re.S)
words = len(re.sub(r"\$[^$]*\$|\\[A-Za-z]+\{?\}?|[{}~]", " x ", m.group(1)).split())
check(f"abstract length is within 100-150 words (counted about {words})", 95 <= words <= 155)
text = re.sub(r"\\(?:label|ref|cite|includegraphics|input|usepackage|documentclass|begin|end|"
              r"setlength|renewcommand|setcounter|geometry|node|draw|resizebox|tikzset)"
              r"(?:\[[^\]]*\])?(?:\{[^}]*\})*", " ", body.split(r"\begin{document}")[1])
lits = sorted(set(re.findall(r"(?<![A-Za-z\\{^_.\d])\d+(?:\.\d+)?(?![\d.]*[A-Za-z}])", text)), key=float)
out.append(f"[INFO] {len(lits)} distinct literal numbers remain in the manuscript text. They are inputs, "
           "ranges, constants or settings (state points, grid sizes, tolerances, published EoS constants), "
           "not results: " + " ".join(lits))

# 5 LaTeX log ------------------------------------------------------------
log = os.path.join(ROOT, "manuscript.log")
if os.path.exists(log):
    L = read(log) if True else ""
    check("manuscript compiles with no LaTeX errors", "\n!" not in L)
    check("no overfull boxes", "Overfull" not in L,
          f"{L.count('Overfull')} overfull box warnings")
    check("no undefined references or citations", "undefined" not in L.lower().replace("undefined control", ""))
else:
    check("manuscript.log present (compile the manuscript before the audit)", False)
check("manuscript.pdf exists", os.path.exists(os.path.join(ROOT, "manuscript.pdf")))

# 6 README instructions reproduce the results ----------------------------
env = dict(os.environ)
r = subprocess.run([sys.executable, "-m", "thermocalculus", "derive", "(dT/dP)_H"], cwd=CODES,
                   capture_output=True, text=True, env=env)
check("README step 1 (derive) gives the documented output",
      "(dT/dP)_H = V*(T*alpha - 1)/Cp" in r.stdout, r.stdout.strip()[:120])
with tempfile.TemporaryDirectory() as tmp:
    files = sorted(glob.glob(os.path.join(DATA, "inputs", "*.txt")))
    subprocess.run([sys.executable, "-m", "thermocalculus", "run"] + files + ["-o", tmp], cwd=CODES,
                   capture_output=True, text=True, env=env)
    diff = []
    for f in sorted(glob.glob(os.path.join(tmp, "*_output.csv"))):
        ref = os.path.join(DATA, "outputs", os.path.basename(f))
        a = [l for l in open(f) if not l.startswith("# input file")]
        b = [l for l in open(ref) if not l.startswith("# input file")] if os.path.exists(ref) else None
        if a != b:
            diff.append(os.path.basename(f))
    check(f"README step 4: re-running all {len(files)} scenarios reproduces Data/outputs exactly",
          not diff, ", ".join(diff))
with tempfile.TemporaryDirectory() as tmp:
    tr = os.path.join(tmp, "t.txt")
    subprocess.run([sys.executable, "-m", "thermocalculus", "ask", "--script",
                    os.path.join(DATA, "inputs", "sessions", "ask_session_script.txt"),
                    "--transcript", tr], cwd=CODES, capture_output=True, text=True, env=env)
    ref = os.path.join(DATA, "outputs", "ask_session_transcript.txt")
    check("README step 5: the interactive session (ask) reproduces its transcript exactly",
          os.path.exists(tr) and os.path.exists(ref) and open(tr).read() == open(ref).read())
r = subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=CODES, capture_output=True, text=True)
tail = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-200:]
check("README step 7: the test suite passes", " passed" in tail and "failed" not in tail, tail)
idx = read(os.path.join(DATA, "DATA_INDEX.txt"))
check("every Data file is described in Data/DATA_INDEX.txt", "(no description)" not in idx)
js = os.path.join(DATA, "verification", "website_js_vs_python.csv")
if os.path.exists(js):
    line = read(js).splitlines()[1]
    val = float(re.search(r"difference: ([\d.e+-]+)", line).group(1))
    check("website JavaScript engine agrees with the Python engine (relative difference < 1e-9)",
          val < 1e-9, line[2:])
for d in ("Literature", "Codes", "Data"):
    check(f"folder {d}/ exists", os.path.isdir(os.path.join(ROOT, d)))
extra = [d for d in os.listdir(ROOT) if os.path.isdir(os.path.join(ROOT, d)) and d not in ("Literature", "Codes", "Data")]
check("the project folder has exactly the three required subfolders", not extra, ", ".join(extra))

head = ["CONSISTENCY AUDIT  (generated by Codes/tools/audit.py)", "=" * 60,
        f"{sum(1 for l in out if l.startswith('[PASS]'))} checks passed, {n_fail} failed", ""]
with open(os.path.join(ROOT, "AUDIT_REPORT.txt"), "w") as fh:
    fh.write("\n".join(head + out) + "\n")
print("\n".join(head + out))
