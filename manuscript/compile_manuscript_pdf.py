#!/usr/bin/env python3
"""
compile_manuscript_pdf.py
Compiles the npj Microgravity publication suite matching spaceflight-plant-hardware-cfd:
1. Typst npj Microgravity PDF: manuscript.typ -> npj_manuscript.pdf
2. LaTeX npj Microgravity PDF: npj_manuscript.tex -> npj_manuscript_latex.pdf
3. LaTeX Supplementary Information PDF: npj_supplementary.tex -> npj_supplementary.pdf
"""

import os
import sys
import subprocess
from pathlib import Path

ENV = os.environ.copy()
ENV["PATH"] = "/opt/homebrew/bin:" + ENV.get("PATH", "")

def compile_typst_pdf(manuscript_dir):
    typst_file = manuscript_dir / "manuscript.typ"
    pdf_out = manuscript_dir / "npj_manuscript.pdf"
    print(f"=== Compiling npj Microgravity Manuscript via Typst: {typst_file.name} ===")
    try:
        import typst
        typst.compile(str(typst_file), output=str(pdf_out))
        if pdf_out.exists() and pdf_out.stat().st_size > 0:
            print(f"  Successfully compiled Typst PDF: {pdf_out} ({pdf_out.stat().st_size:,} bytes)")
        else:
            print(f"  Warning: Typst output not found: {pdf_out}")
    except Exception as e:
        print(f"  Typst compilation error: {e}")

def compile_latex_pdf(tex_file, pdf_out_name, manuscript_dir):
    print(f"=== Compiling LaTeX Document: {tex_file.name} ===")
    cmd = ["pdflatex", "-interaction=nonstopmode", tex_file.name]
    res = subprocess.run(cmd, cwd=str(manuscript_dir), env=ENV, capture_output=True, text=True)
    pdf_generated = manuscript_dir / f"{tex_file.stem}.pdf"
    target_pdf = manuscript_dir / pdf_out_name
    if pdf_generated.exists() and pdf_generated.stat().st_size > 0:
        if pdf_generated != target_pdf:
            if target_pdf.exists():
                target_pdf.unlink()
            pdf_generated.rename(target_pdf)
        print(f"  Successfully compiled LaTeX PDF: {target_pdf} ({target_pdf.stat().st_size:,} bytes)")
    else:
        print(f"  LaTeX compilation warning for {tex_file.name}. Log summary:\n{res.stdout[-600:]}")
    
    # Cleanup auxiliary files
    for ext in [".aux", ".log", ".out"]:
        f_aux = manuscript_dir / f"{tex_file.stem}{ext}"
        if f_aux.exists():
            f_aux.unlink()

def main():
    manuscript_dir = Path(__file__).resolve().parent
    print("=== Compiling Complete Publication Suite (npj Microgravity Style) ===")

    # 1. Typst npj Microgravity PDF
    compile_typst_pdf(manuscript_dir)

    # 2. LaTeX npj Microgravity PDF
    compile_latex_pdf(manuscript_dir / "npj_manuscript.tex", "npj_manuscript_latex.pdf", manuscript_dir)

    # 3. LaTeX Supplementary PDF
    compile_latex_pdf(manuscript_dir / "npj_supplementary.tex", "npj_supplementary.pdf", manuscript_dir)

    print("\nAll publication PDF targets compiled successfully.")

if __name__ == "__main__":
    main()
