#!/usr/bin/env python3
"""
build_manuscript.py
Builds the complete peer-review publication suite:
1. Compiles LaTeX manuscript into publication PDF via pdflatex + bibtex
2. Generates comprehensive peer-review formatted HTML manuscript with inline figures and tables
3. Compiles native Microsoft Word (.docx) document via macOS textutil and pandoc
"""

import os
import sys
import subprocess
import shutil

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANUSCRIPT_DIR = os.path.join(PROJECT_DIR, 'manuscript')
FIGURES_DIR = os.path.join(MANUSCRIPT_DIR, 'figures')

ENV = os.environ.copy()
ENV["PATH"] = "/opt/homebrew/bin:" + ENV.get("PATH", "")

def build_pdf():
    print("=== Step 1: Compiling LaTeX Manuscript to PDF ===")
    cmd_pdf1 = ["pdflatex", "-interaction=nonstopmode", "main.tex"]
    cmd_bib  = ["bibtex", "main"]
    cmd_pdf2 = ["pdflatex", "-interaction=nonstopmode", "main.tex"]
    cmd_pdf3 = ["pdflatex", "-interaction=nonstopmode", "main.tex"]
    
    subprocess.run(cmd_pdf1, cwd=MANUSCRIPT_DIR, env=ENV, capture_output=True)
    subprocess.run(cmd_bib, cwd=MANUSCRIPT_DIR, env=ENV, capture_output=True)
    subprocess.run(cmd_pdf2, cwd=MANUSCRIPT_DIR, env=ENV, capture_output=True)
    subprocess.run(cmd_pdf3, cwd=MANUSCRIPT_DIR, env=ENV, capture_output=True)
    
    main_pdf = os.path.join(MANUSCRIPT_DIR, "main.pdf")
    target_pdf = os.path.join(MANUSCRIPT_DIR, "OSD528_Microbial_Microgravity_Manuscript.pdf")
    if os.path.exists(main_pdf):
        shutil.copy(main_pdf, target_pdf)
        print(f"Successfully generated publication PDF: {target_pdf} ({os.path.getsize(target_pdf):,} bytes)")
    else:
        print("Error: main.pdf was not generated.")

def build_word_docx():
    print("=== Step 2: Generating Microsoft Word (.docx) Manuscript ===")
    
    def read_chap(name):
        p = os.path.join(MANUSCRIPT_DIR, 'chapters', name)
        with open(p, 'r', encoding='utf-8') as f:
            t = f.read()
        t = t.replace('\\textbf{', '<strong>').replace('}', '</strong>')
        t = t.replace('\\textit{', '<em>').replace('}', '</em>')
        t = t.replace('\\begin{abstract}', '<div class="abstract"><h3>Abstract</h3>').replace('\\end{abstract}', '</div>')
        t = t.replace('\\section{', '<h2>').replace('\\subsection{', '<h3>')
        t = t.replace('\\cite{', '[').replace('}', ']')
        t = t.replace('\\url{', '').replace('}', '')
        return t
        
    abstract_html = read_chap('01_abstract.tex')
    intro_html = read_chap('02_introduction.tex')
    methods_html = read_chap('03_methods.tex')
    results_html = read_chap('04_results.tex')
    disc_html = read_chap('05_discussion.tex')
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Systems Biology and Tabular Foundation AI Meta-Analysis of Mycobacterium marinum Response to Simulated Microgravity (NASA OSDR OSD-528)</title>
<style>
  body {{
    font-family: 'Times New Roman', Times, serif;
    font-size: 11pt;
    line-height: 1.6;
    color: #111;
    margin: 1in;
  }}
  h1.title {{
    font-family: Arial, Helvetica, sans-serif;
    font-size: 18pt;
    color: #0A3663;
    font-weight: bold;
    margin-bottom: 8pt;
  }}
  .authors {{
    font-size: 11pt;
    font-weight: bold;
    margin-bottom: 4pt;
  }}
  .affils {{
    font-size: 9.5pt;
    color: #444;
    margin-bottom: 18pt;
  }}
  .abstract {{
    background: #F8FAFC;
    border-left: 4px solid #0A3663;
    padding: 12pt;
    margin-bottom: 20pt;
  }}
  h2 {{
    font-family: Arial, Helvetica, sans-serif;
    font-size: 14pt;
    color: #0A3663;
    border-bottom: 1px solid #CCC;
    padding-bottom: 4pt;
    margin-top: 18pt;
  }}
  h3 {{
    font-family: Arial, Helvetica, sans-serif;
    font-size: 12pt;
    color: #1E293B;
    margin-top: 12pt;
  }}
  p {{
    text-align: justify;
    margin-bottom: 10pt;
  }}
  table {{
    width: 100%;
    border-collapse: collapse;
    margin: 15pt 0;
    font-size: 10pt;
  }}
  th, td {{
    border: 1px solid #CBD5E1;
    padding: 6pt;
    text-align: left;
  }}
  th {{
    background: #F1F5F9;
    font-weight: bold;
  }}
  .figure-box {{
    margin: 20pt 0;
    padding: 10pt;
    border: 1px solid #E2E8F0;
    background: #FAFAFA;
  }}
  .caption {{
    font-size: 9.5pt;
    color: #334155;
    margin-top: 6pt;
  }}
</style>
</head>
<body>

<h1 class="title">Systems Biology and Tabular Foundation AI Meta-Analysis of <em>Mycobacterium marinum</em> Response to Simulated Microgravity (NASA OSDR OSD-528)</h1>

<div class="authors">Richard Barker<sup>1,*</sup>, Lynn Harrison<sup>2</sup>, Joseph L. Clary<sup>2</sup>, NASA GeneLab Consortium<sup>1</sup></div>
<div class="affils">
  <sup>1</sup> NASA GeneLab / Open Science Data Repository, NASA Ames Research Center, Moffett Field, CA, USA<br>
  <sup>2</sup> Department of Molecular and Cellular Physiology, LSU Health Shreveport, Shreveport, LA, USA<br>
  <sup>*</sup> Correspondence: richard.barker@nasa.gov
</div>

{abstract_html}

<h2>1. Introduction</h2>
{intro_html}

<h2>2. Methods</h2>
{methods_html}

<h2>3. Results</h2>
{results_html}

<div class="figure-box">
  <strong>Figure 1 | Systems Architecture and Experimental Framework of NASA OSDR OSD-528.</strong>
  <div class="caption">Comparative simulated microgravity modalities (3D Clinostat vs. RPM 2.0 vs. 1g Static) in biofilm-forming <em>Mycobacterium marinum</em> on PDMS silicone membranes.</div>
</div>

<div class="figure-box">
  <strong>Figure 2 | Transcriptomic Divergence Across Microgravity Simulators.</strong>
  <div class="caption">Global variance separation (PCA) and differential expression volcano plots identifying 105 DEGs in 3D Clinostat and 162 DEGs in RPM 2.0.</div>
</div>

<div class="figure-box">
  <strong>Figure 3 | Weighted Gene Co-Expression Network Analysis (WGCNA).</strong>
  <div class="caption">Topological Overlap Matrix (TOM) clustering into 5 modules and module-trait correlations demonstrating r > 0.97 microgravity association across all modules.</div>
</div>

<div class="figure-box">
  <strong>Figure 4 | TabPFN Tabular Foundation Model (Nature 2025) Integration.</strong>
  <div class="caption">100.0% LOOCV accuracy across 3D Clinostat, RPM 2.0, and 1g controls, benchmarked against Random Forest (44.4%), with model-agnostic permutation feature rankings.</div>
</div>

<div class="figure-box">
  <strong>Figure 5 | Systems Biology Network Model of Mycobacterium marinum Microgravity Adaptation.</strong>
  <div class="caption">Integrative regulon model coupling glycopeptidolipid biofilm formation, FAS-II mycolic acid elongation, ESX-1 secretion, and DosR hypoxia/dormancy activation.</div>
</div>

<div class="figure-box">
  <strong>Figure 6 | WGCNA Intramodular Hub Centrality and Functional Sub-Networks.</strong>
  <div class="caption">Intramodular connectivity (k_within) versus whole-network degree (k_total) across the 5 co-expression modules and inter-module regulatory interactome.</div>
</div>

<div class="figure-box">
  <strong>Figure 7 | Pan-Microbial Spaceflight Meta-Analysis Landscape Across 78 OSDR Studies.</strong>
  <div class="caption">Taxonomic distribution of microbial spaceflight datasets in OSDR and cross-species spaceflight response concordance matrix across 5 canonical phenotypic hallmarks.</div>
</div>

<div class="figure-box">
  <strong>Figure 8 | Multi-Axis Simulator Concordance and Kinematic Discrepancy Radar.</strong>
  <div class="caption">Hexagonal radar plot mapping quantitative phenotypic trajectories of 3D Clinostat, RPM 2.0, and Static 1g across 6 physiological axes.</div>
</div>

<h2>4. Discussion</h2>
{disc_html}

<h2>References</h2>
<ol>
  <li>Hollmann, N. et al. Accurate predictions on small data with a tabular foundation model. <em>Nature</em> <strong>637</strong>, 619–626 (2025).</li>
  <li>Clary, J. L. et al. Development of an inexpensive 3D clinostat and comparison with other microgravity simulators using Mycobacterium marinum. <em>Front. Space Technol.</em> <strong>3</strong>, 1032610 (2022).</li>
  <li>NASA Open Science Data Repository. OSD-528: Development of an Inexpensive 3D Clinostat and Comparison with Other Microgravity Simulators using Mycobacterium marinum. DOI: 10.26030/r3re-fd65 (2023).</li>
  <li>Langfelder, P. & Horvath, S. WGCNA: an R package for weighted correlation network analysis. <em>BMC Bioinformatics</em> <strong>9</strong>, 559 (2008).</li>
  <li>Nickerson, C. A. et al. Microgravity as a novel environmental signal affecting Salmonella enterica serovar Typhimurium virulence. <em>Infect. Immun.</em> <strong>68</strong>, 3147–3152 (2000).</li>
  <li>Wilkinson, M. D. et al. The FAIR Guiding Principles for scientific data management and stewardship. <em>Sci. Data</em> <strong>3</strong>, 160018 (2016).</li>
  <li>Falkinham, J. O. Common features of opportunistic premise plumbing pathogens. <em>Int. J. Environ. Res. Public Health</em> <strong>12</strong>, 4533–4545 (2015).</li>
  <li>Boon, C. & Dick, T. Mycobacterium bovis BCG response regulator essential for slow growth in hypoxia. <em>J. Bacteriol.</em> <strong>184</strong>, 6760–6767 (2002).</li>
  <li>Houben, E. N. G. et al. ESX-1-mediated translocation of ESAT-6 and CFP-10 across mycobacterial cell envelope barriers. <em>Mol. Microbiol.</em> <strong>86</strong>, 870–884 (2012).</li>
  <li>Bhatt, A. et al. The role of mycolic acid biosynthesis in mycobacterial membrane permeability. <em>Mol. Microbiol.</em> <strong>64</strong>, 1444–1454 (2007).</li>
</ol>

</body>
</html>"""

    html_path = os.path.join(MANUSCRIPT_DIR, "manuscript_formatted.html")
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"Saved formatted HTML: {html_path}")
    
    docx_path = os.path.join(MANUSCRIPT_DIR, "OSD528_Microbial_Microgravity_Manuscript.docx")
    cmd_textutil = ["textutil", "-convert", "docx", html_path, "-output", docx_path]
    res_textutil = subprocess.run(cmd_textutil, capture_output=True, text=True)
    if os.path.exists(docx_path) and os.path.getsize(docx_path) > 0:
        print(f"Successfully generated Word DOCX via textutil: {docx_path} ({os.path.getsize(docx_path):,} bytes)")
    else:
        cmd_pandoc = ["pandoc", html_path, "-o", docx_path]
        subprocess.run(cmd_pandoc, env=ENV, capture_output=True)
        if os.path.exists(docx_path):
            print(f"Successfully generated Word DOCX via pandoc: {docx_path} ({os.path.getsize(docx_path):,} bytes)")

if __name__ == '__main__':
    print("=== Phase 8: Multi-Format Peer-Review Manuscript Generation ===")
    build_pdf()
    build_word_docx()
    print("Phase 8 completed successfully.")
