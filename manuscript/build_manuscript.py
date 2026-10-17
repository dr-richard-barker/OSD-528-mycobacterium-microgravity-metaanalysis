#!/usr/bin/env python3
"""
build_manuscript.py
Build every manuscript output from one text.

The manuscript text lives only in chapters/*.tex. Two LaTeX wrappers lay it out:
  main.tex             article, two columns -> main.pdf, copied to
                       OSD528_Microbial_Microgravity_Manuscript.pdf
  npj_manuscript.tex   journal-style layout (wlpeerj.cls) -> npj_manuscript_latex.pdf
plus the supplementary information:
  npj_supplementary.tex -> npj_supplementary.pdf
pandoc converts main.tex (with references.bib) to:
  manuscript_formatted.html and OSD528_Microbial_Microgravity_Manuscript.docx

Each LaTeX document is built pdflatex -> bibtex -> pdflatex -> pdflatex, and the build
fails if the final log still reports an undefined citation or reference, so no PDF can
ship with "??" in it. Figures are referenced without an extension: pdflatex picks the
vector PDF, pandoc the PNG.
"""

import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ENV = dict(os.environ, PATH="/opt/homebrew/bin:/Library/TeX/texbin:" + os.environ.get("PATH", ""))
AUX = (".aux", ".bbl", ".blg", ".log", ".out", ".toc", ".fls", ".fdb_latexmk", ".synctex.gz")


def run(cmd):
    res = subprocess.run(cmd, cwd=HERE, env=ENV, capture_output=True, text=True)
    return res.returncode, res.stdout + res.stderr


def latex(stem, jobname=None):
    """Build stem.tex into jobname.pdf (default: stem.pdf).

    -jobname writes every intermediate and the PDF under the job name, so building
    npj_manuscript.tex as npj_manuscript_latex never creates or overwrites a file
    called npj_manuscript.pdf.
    """
    job = jobname or stem
    print(f"=== {stem}.tex -> {job}.pdf")
    tex = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"-jobname={job}", f"{stem}.tex"]
    for step in (tex, ["bibtex", job], tex, tex):
        if step[0] == "bibtex" and not re.search(r"\\bibdata", open(os.path.join(HERE, f"{job}.aux")).read()):
            continue
        code, out = run(step)
        if code != 0:
            sys.exit(f"{' '.join(step)} failed:\n{out[-2500:]}")
    log = open(os.path.join(HERE, f"{job}.log"), encoding="latin-1").read()
    problems = re.findall(r"(?:Citation|Reference) `[^']+' (?:on page \d+ )?undefined", log)
    if problems or "There were undefined references" in log:
        sys.exit(f"{job}: undefined citations/references remain: {sorted(set(problems))}")
    blg = os.path.join(HERE, f"{job}.blg")
    if os.path.exists(blg) and re.search(r"^Warning--I didn't find a database entry", open(blg).read(), re.M):
        sys.exit(f"{job}: bibtex could not find a cited entry; see {job}.blg")
    print(f"  ok: {job}.pdf")
    for ext in AUX:
        p = os.path.join(HERE, job + ext)
        if os.path.exists(p):
            os.remove(p)


def pandoc():
    common = ["pandoc", "main.tex", "--from=latex", "--citeproc", "--bibliography=references.bib",
              "--default-image-extension=png", "--resource-path=.:figures",
              "--metadata=link-citations:true"]
    for out, extra in (("manuscript_formatted.html", ["--standalone", "--mathml",
                                                      "--metadata=pagetitle:OSD-528 re-analysis"]),
                       ("OSD528_Microbial_Microgravity_Manuscript.docx", [])):
        print(f"=== pandoc -> {out}")
        code, msg = run(common + extra + ["-o", out])
        if code != 0:
            sys.exit(f"pandoc failed for {out}:\n{msg[-2500:]}")
        missing = re.findall(r"Citeproc: citation (\S+) not found", msg)
        if missing:
            sys.exit(f"pandoc: citations not found in references.bib: {missing}")
        print(f"  ok: {out}")


def main():
    latex("main")
    shutil.copyfile(os.path.join(HERE, "main.pdf"),
                    os.path.join(HERE, "OSD528_Microbial_Microgravity_Manuscript.pdf"))
    latex("npj_manuscript", "npj_manuscript_latex")
    latex("npj_supplementary")
    pandoc()


if __name__ == "__main__":
    main()
