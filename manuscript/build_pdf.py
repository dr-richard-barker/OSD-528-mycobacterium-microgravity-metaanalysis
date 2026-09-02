#!/usr/bin/env python3
"""Build script for compiling npj Microgravity Typst manuscript to PDF."""
import os
from pathlib import Path
import typst

def main():
    root = Path(__file__).resolve().parent
    typst_path = root / "manuscript.typ"
    pdf_path = root / "npj_manuscript.pdf"
    
    print(f"Compiling {typst_path} -> {pdf_path}...")
    typst.compile(str(typst_path), output=str(pdf_path))
    if pdf_path.exists():
        print(f"Successfully compiled {pdf_path} ({pdf_path.stat().st_size:,} bytes)")

if __name__ == "__main__":
    main()
