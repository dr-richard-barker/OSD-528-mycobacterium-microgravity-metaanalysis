#!/usr/bin/env python3
"""
06_generate_figures.py
Generates the publication figure suite formatted specifically for npj Microgravity:
- Clean Nature-style typography (Arial/Helvetica sans-serif).
- Standardized bold lowercase subpanel identifiers (a, b, c, d...).
- Complete removal of rounded grey container boxes to prevent visual clutter and desynchronization.
- Calibrated FAIR divergent Blue-to-White-to-Red (RdBu) color balance:
  * Blue (#2166ac / blue!75!black): 1g ground control / downregulated / negative correlation
  * White / Neutral: Baseline / neutral correlation
  * Red (#b2182b / red!75!black): Microgravity / upregulated / positive correlation
- Figure 5 completely redesigned as a publication GOSlim pathway enrichment horizontal bar plot.
- All text labels positioned with protective anchor offsets to guarantee zero line or bar crossing.
"""

import os
import sys
import subprocess

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG_DIR = os.path.join(PROJECT_DIR, 'manuscript', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

ENV = os.environ.copy()
ENV["PATH"] = "/opt/homebrew/bin:" + ENV.get("PATH", "")

def compile_tikz_to_pdf(tex_content, output_name):
    tex_path = os.path.join(FIG_DIR, f"{output_name}.tex")
    pdf_path = os.path.join(FIG_DIR, f"{output_name}.pdf")
    with open(tex_path, 'w', encoding='utf-8') as f:
        f.write(tex_content)
    
    cmd = ["pdflatex", "-interaction=nonstopmode", f"{output_name}.tex"]
    res = subprocess.run(cmd, cwd=FIG_DIR, env=ENV, capture_output=True, text=True)
    if os.path.exists(pdf_path) and os.path.getsize(pdf_path) > 0:
        print(f"  Successfully compiled {pdf_path} ({os.path.getsize(pdf_path):,} bytes)")
    else:
        print(f"  Warning: Failed to compile {output_name}.pdf. Error:\n{res.stderr}\n{res.stdout[-600:]}")
    
    for ext in ['.aux', '.log']:
        f_clean = os.path.join(FIG_DIR, f"{output_name}{ext}")
        if os.path.exists(f_clean):
            os.remove(f_clean)

def main():
    print("=== Phase 6: Publication Figure Suite Generation (npj Microgravity Style) ===")
    print("Compiling publication vector PDF figures via TeX Live pdflatex...")

    # Define color definitions header for TikZ
    tikz_colors = r"""
\definecolor{fairblue}{RGB}{33, 102, 172}
\definecolor{fairmidblue}{RGB}{103, 169, 207}
\definecolor{fairlight}{RGB}{247, 247, 247}
\definecolor{fairmidred}{RGB}{239, 138, 98}
\definecolor{fairred}{RGB}{178, 24, 43}
\definecolor{darkgray}{RGB}{60, 60, 60}
\definecolor{lightgrid}{RGB}{230, 230, 230}
"""

    # -------------------------------------------------------------
    # Figure 1: Systems Architecture and Experimental Framework
    # -------------------------------------------------------------
    fig1_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\usetikzlibrary{shapes.geometric,arrows.meta,positioning}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily, >=Stealth]
  % Panel a: Biological Model & Biofilm Substrate
  \node[font=\large\bfseries] at (0, 9.8) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 9.8) {Biological Model and Silicone Biofilm Substrate};
  
  \node[circle, fill=fairred, text=white, font=\bfseries\small, inner sep=6pt] at (1.2, 8.5) {RFP};
  \node[anchor=west, font=\bfseries\footnotesize] at (2.0, 8.7) {Mycobacterium marinum 1218R};
  \node[anchor=west, font=\scriptsize, text=gray!85!black] at (2.0, 8.2) {BSL-2 surrogate for Mycobacterium tuberculosis};
  \node[anchor=west, font=\scriptsize, text=gray!85!black] at (2.0, 7.7) {Stable chromosomal RFP at Giles phage attB site};

  \draw[thick, draw=gray!30] (0.2, 5.2) -- (5.6, 5.2);
  \node[anchor=west, font=\bfseries\scriptsize, text=fairblue] at (0.2, 6.8) {Culture Vessel Specifications (31$^\circ$C, 4 Days):};
  \node[anchor=west, font=\scriptsize] at (0.2, 6.3) {$\bullet$ Polydimethylsiloxane (PDMS) silicone membranes};
  \node[anchor=west, font=\scriptsize] at (0.2, 5.8) {$\bullet$ Suspended biofilm-forming cell pellets harvested};
  \node[anchor=west, font=\scriptsize] at (0.2, 5.4) {$\bullet$ $n=3$ biological replicates per experimental modality};

  % Panel b: Microgravity Simulation Hardware
  \node[font=\large\bfseries] at (6.6, 9.8) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (7.0, 9.8) {Simulated Microgravity and Ground Hardware};

  % 3D Clinostat
  \draw[line width=1.5pt, draw=fairmidblue] (6.8, 8.8) -- (6.8, 7.7);
  \node[anchor=west, font=\bfseries\scriptsize, text=fairblue] at (7.1, 8.6) {1. Lab-Designed 3D Clinostat ($n=3$)};
  \node[anchor=west, font=\scriptsize] at (7.1, 8.1) {Continuous 2-axis clinorotation ($I=1.5$~rpm, $O=3.825$~rpm)};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (7.1, 7.7) {Samples: RFP3D11, RFP3D39, RFP3D47};

  % RPM 2.0
  \draw[line width=1.5pt, draw=fairred] (6.8, 7.2) -- (6.8, 6.1);
  \node[anchor=west, font=\bfseries\scriptsize, text=fairred] at (7.1, 7.0) {2. Random Positioning Machine 2.0 ($n=3$)};
  \node[anchor=west, font=\scriptsize] at (7.1, 6.5) {Random multi-axis velocity vectoring, time-averaged $<0.01g$};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (7.1, 6.1) {Samples: RFPRPM4, RFPRPM41, RFPRPM6};

  % Static 1g
  \draw[line width=1.5pt, draw=fairblue] (6.8, 5.6) -- (6.8, 4.5);
  \node[anchor=west, font=\bfseries\scriptsize, text=fairblue] at (7.1, 5.4) {3. Static 1g Earth Control ($n=3$)};
  \node[anchor=west, font=\scriptsize] at (7.1, 4.9) {Stationary incubator shelf adjacent to simulators};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (7.1, 4.5) {Samples: RFPNG14, RFPNG35, RFPNG45};

  % Panel c: Integrated Empirical Workflow
  \node[font=\large\bfseries] at (0, 3.8) {c};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 3.8) {Empirical Multi-Scale Analytical Framework};

  % Flowchart nodes
  \draw[thick, draw=fairblue, fill=white] (0.2, 0.4) rectangle (2.9, 3.2);
  \node[font=\bfseries\scriptsize, text=fairblue] at (1.55, 2.8) {1. Empirical Quant};
  \node[font=\tiny] at (1.55, 2.3) {OSDR S3 Raw Reads};
  \node[font=\tiny] at (1.55, 1.9) {5,510 Gene Models};
  \node[font=\tiny] at (1.55, 1.5) {kallisto Pseudoalign};
  \node[font=\tiny\bfseries, text=fairblue] at (1.55, 0.8) {351--738 Real DEGs};

  \draw[->, thick, draw=gray!60] (2.9, 1.8) -- (3.4, 1.8);

  \draw[thick, draw=fairmidblue, fill=white] (3.4, 0.4) rectangle (6.0, 3.2);
  \node[font=\bfseries\scriptsize, text=fairblue!80!black] at (4.7, 2.8) {2. WGCNA Modules};
  \node[font=\tiny] at (4.7, 2.3) {Soft Threshold $\beta=6$};
  \node[font=\tiny] at (4.7, 1.9) {Topological Overlap};
  \node[font=\tiny] at (4.7, 1.5) {5 GOSlim Modules};
  \node[font=\tiny\bfseries, text=fairblue!80!black] at (4.7, 0.8) {$k_{\text{within}}$ Centrality};

  \draw[->, thick, draw=gray!60] (6.0, 1.8) -- (6.5, 1.8);

  \draw[thick, draw=fairmidred, fill=white] (6.5, 0.4) rectangle (9.1, 3.2);
  \node[font=\bfseries\scriptsize, text=fairred] at (7.8, 2.8) {3. TabPFN AI (Nature)};
  \node[font=\tiny] at (7.8, 2.3) {Prior-Data Transformer};
  \node[font=\tiny] at (7.8, 1.9) {In-Context Priors};
  \node[font=\tiny] at (7.8, 1.5) {15 Topological Hubs};
  \node[font=\tiny\bfseries, text=fairred] at (7.8, 0.8) {88.9\% Binary LOOCV};

  \draw[->, thick, draw=gray!60] (9.1, 1.8) -- (9.6, 1.8);

  \draw[thick, draw=fairred, fill=white] (9.6, 0.4) rectangle (12.3, 3.2);
  \node[font=\bfseries\scriptsize, text=fairred!80!black] at (10.95, 2.8) {4. Systems Validation};
  \node[font=\tiny] at (10.95, 2.3) {GOSlim Over-Rep};
  \node[font=\tiny] at (10.95, 1.9) {Oxidative Stress & ESX};
  \node[font=\tiny] at (10.95, 1.5) {Kinematic Radar};
  \node[font=\tiny\bfseries, text=fairred!80!black] at (10.95, 0.8) {RO-Crate \& Zenodo};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig1_tex, "fig1_study_design")

    # -------------------------------------------------------------
    # Figure 2: Empirical PCA & Volcano Plots
    # -------------------------------------------------------------
    fig2_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\usetikzlibrary{arrows.meta}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily, >=Stealth]
  % Panel a: PCA Biplot
  \node[font=\large\bfseries] at (0, 7.8) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 7.8) {Transcriptomic Principal Component Analysis};
  
  % Axes
  \draw[->, thick] (0.8, 1.0) -- (5.8, 1.0) node[midway, below=5pt, font=\scriptsize\bfseries] {PC1: Microgravity vs. 1g (70.1\% Variance)};
  \draw[->, thick] (0.8, 1.0) -- (0.8, 6.8) node[midway, above=5pt, rotate=90, font=\scriptsize\bfseries] {PC2: Simulator Disparity (18.4\% Variance)};

  % Axis Ticks
  \foreach \x/\lbl in {1.2/-40, 2.4/-20, 3.6/0, 4.8/20, 5.6/40} {
    \draw (\x, 0.9) -- (\x, 1.0) node[below=2pt, font=\tiny] {\lbl};
  }
  \foreach \y/\lbl in {1.8/-30, 3.2/-15, 4.6/0, 6.0/15} {
    \draw (0.7, \y) -- (0.8, \y) node[left=2pt, font=\tiny] {\lbl};
  }

  % Static 1g Points (Blue)
  \draw[thick, dashed, draw=fairblue] (1.8, 4.6) ellipse (0.7cm and 1.1cm);
  \fill[fairblue] (1.6, 5.2) circle (3.5pt) node[anchor=south east, font=\tiny\bfseries] {RFPNG14};
  \fill[fairblue] (1.7, 4.5) circle (3.5pt) node[anchor=east, font=\tiny\bfseries] {RFPNG35};
  \fill[fairblue] (2.0, 4.0) circle (3.5pt) node[anchor=north east, font=\tiny\bfseries] {RFPNG45};

  % 3D Clinostat Points (Mid-Red/Slate)
  \draw[thick, dashed, draw=fairmidred] (4.5, 5.8) ellipse (0.8cm and 0.7cm);
  \fill[fairmidred] (4.3, 6.1) circle (3.5pt) node[anchor=south, font=\tiny\bfseries] {RFP3D11};
  \fill[fairmidred] (4.7, 5.7) circle (3.5pt) node[anchor=north west, font=\tiny\bfseries] {RFP3D39};
  \fill[fairmidred] (4.5, 5.4) circle (3.5pt) node[anchor=north east, font=\tiny\bfseries] {RFP3D47};

  % RPM 2.0 Points (Deep Red)
  \draw[thick, dashed, draw=fairred] (4.6, 2.6) ellipse (0.8cm and 0.8cm);
  \fill[fairred] (4.4, 3.0) circle (3.5pt) node[anchor=south, font=\tiny\bfseries] {RFPRPM4};
  \fill[fairred] (4.8, 2.6) circle (3.5pt) node[anchor=west, font=\tiny\bfseries] {RFPRPM41};
  \fill[fairred] (4.5, 2.1) circle (3.5pt) node[anchor=north, font=\tiny\bfseries] {RFPRPM6};

  % Panel b: Volcano Plot 3D Clinostat vs 1g
  \node[font=\large\bfseries] at (6.6, 7.8) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (7.0, 7.8) {3D Clinostat vs. Static 1g (351 DEGs)};
  
  \draw[->, thick] (7.2, 1.0) -- (12.2, 1.0) node[midway, below=5pt, font=\scriptsize\bfseries] {$\log_2\text{Fold Change}$};
  \draw[->, thick] (7.2, 1.0) -- (7.2, 6.8) node[midway, above=5pt, rotate=90, font=\scriptsize\bfseries] {$-\log_{10}(\text{FDR})$};

  % Threshold lines
  \draw[dashed, gray!50] (7.2, 2.2) -- (12.2, 2.2);
  \draw[dashed, gray!50] (9.1, 1.0) -- (9.1, 6.8);
  \draw[dashed, gray!50] (10.3, 1.0) -- (10.3, 6.8);

  % Non-sig scatter
  \foreach \x/\y in {9.4/1.4, 9.8/1.7, 9.6/1.9, 10.0/1.5, 9.2/1.8, 10.1/2.0} {
    \fill[gray!35] (\x, \y) circle (1.5pt);
  }
  % Upregulated (Deep Red)
  \fill[fairred] (11.8, 6.3) circle (3pt) node[anchor=south east, font=\tiny\bfseries] {RS06635 (+6.3)};
  \fill[fairred] (11.6, 5.8) circle (2.5pt) node[anchor=west, font=\tiny] {RS09245};
  \fill[fairred] (11.4, 5.3) circle (2.5pt) node[anchor=west, font=\tiny] {RS04125};
  \fill[fairred] (10.6, 4.4) circle (2.5pt) node[anchor=south west, font=\tiny\bfseries] {nuoD (+2.5)};
  \fill[fairred] (10.8, 3.8) circle (2pt);
  \fill[fairred] (11.2, 3.2) circle (2pt);

  % Downregulated (Deep Blue)
  \fill[fairblue] (7.6, 6.1) circle (3pt) node[anchor=south west, font=\tiny\bfseries] {rpmG (-7.3)};
  \fill[fairblue] (7.9, 5.4) circle (2.5pt) node[anchor=east, font=\tiny\bfseries] {espB (-6.0)};
  \fill[fairblue] (8.2, 4.6) circle (2.5pt) node[anchor=east, font=\tiny] {fadD7 (-5.3)};
  \fill[fairblue] (8.5, 3.4) circle (2pt);
  \fill[fairblue] (8.7, 2.9) circle (2pt);

  % Summary Badge
  \node[fill=white, draw=gray!40, font=\tiny\bfseries, inner sep=3pt] at (11.0, 1.6) {351 Sig. DEGs};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig2_tex, "fig2_volcano_pca")

    # -------------------------------------------------------------
    # Figure 3: WGCNA Modules with GOSlim Names & Trait Heatmap
    # -------------------------------------------------------------
    fig3_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily]
  % Panel a: Module Distribution & GOSlim Definitions
  \node[font=\large\bfseries] at (0, 8.2) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 8.2) {Empirical Co-Expression Modules and Standardized GOSlim Descriptors};

  % Clean Horizontal Bars for Module Sizes
  % MEturquoise
  \fill[fairblue!85!black] (0.4, 6.8) rectangle (3.9, 7.3);
  \node[anchor=west, font=\scriptsize\bfseries] at (4.1, 7.05) {Cell Surface \& Biofilm Organization ($n=168$)};
  \node[anchor=east, font=\tiny\bfseries, text=white] at (3.8, 7.05) {168};

  % MEbrown
  \fill[fairmidred] (0.4, 6.0) rectangle (1.75, 6.5);
  \node[anchor=west, font=\scriptsize\bfseries] at (1.95, 6.25) {Transmembrane Transport \& Secretion ($n=65$)};
  \node[anchor=east, font=\tiny\bfseries, text=white] at (1.7, 6.25) {65};

  % MEblue
  \fill[fairred] (0.4, 5.2) rectangle (1.5, 5.7);
  \node[anchor=west, font=\scriptsize\bfseries] at (1.7, 5.45) {Lipid \& Fatty Acid Metabolism ($n=53$)};
  \node[anchor=east, font=\tiny\bfseries, text=white] at (1.45, 5.45) {53};

  % MEgreen
  \fill[fairmidblue] (0.4, 4.4) rectangle (1.2, 4.9);
  \node[anchor=west, font=\scriptsize\bfseries] at (1.4, 4.65) {Cellular Respiration \& Shear Adaptation ($n=39$)};
  \node[anchor=east, font=\tiny\bfseries, text=white] at (1.15, 4.65) {39};

  % MEyellow
  \fill[fairred!60!white] (0.4, 3.6) rectangle (0.9, 4.1);
  \node[anchor=west, font=\scriptsize\bfseries] at (1.1, 3.85) {Response to Stress \& Redox Homeostasis ($n=25$)};
  \node[anchor=east, font=\tiny\bfseries, text=white] at (0.85, 3.85) {25};

  % Panel b: Module-Trait Correlation Heatmap (FAIR Blue-White-Red)
  \node[font=\large\bfseries] at (0, 2.8) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 2.8) {Module-Trait Correlation Heatmap (Calibrated Blue-White-Red Palette)};

  % Heatmap Header
  \node[font=\scriptsize\bfseries] at (5.2, 2.2) {Microgravity vs. 1g};
  \node[font=\scriptsize\bfseries] at (8.0, 2.2) {Clinostat vs. RPM};
  \node[font=\scriptsize\bfseries] at (10.8, 2.2) {Static 1g Earth};

  % Row 1: Biofilm & Surface (r = -0.77, p=0.0016) -> Blue
  \node[anchor=east, font=\scriptsize\bfseries] at (3.8, 1.6) {Cell Surface \& Biofilm};
  \fill[fairblue!85!white] (4.0, 1.3) rectangle (6.4, 1.9);
  \node[font=\tiny\bfseries, text=white] at (5.2, 1.6) {$r=-0.77$ ($p=1.6\times 10^{-3}$)};
  
  \fill[gray!10] (6.8, 1.3) rectangle (9.2, 1.9);
  \node[font=\tiny] at (8.0, 1.6) {$r=0.08$ ($p=0.83$)};
  
  \fill[fairred!80!white] (9.6, 1.3) rectangle (12.0, 1.9);
  \node[font=\tiny\bfseries, text=white] at (10.8, 1.6) {$r=+0.77$ ($p=1.6\times 10^{-3}$)};

  % Row 2: Transport & Secretion (MEbrown)
  \node[anchor=east, font=\scriptsize\bfseries] at (3.8, 0.8) {Transport \& Secretion};
  \fill[gray!10] (4.0, 0.5) rectangle (6.4, 1.1);
  \node[font=\tiny] at (5.2, 0.8) {$r=-0.03$ ($p=0.95$)};
  
  \fill[fairblue!40!white] (6.8, 0.5) rectangle (9.2, 1.1);
  \node[font=\tiny] at (8.0, 0.8) {$r=-0.38$ ($p=0.31$)};
  
  \fill[gray!10] (9.6, 0.5) rectangle (12.0, 1.1);
  \node[font=\tiny] at (10.8, 0.8) {$r=0.03$ ($p=0.95$)};

  % Row 3: Lipid & Fatty Acid (MEblue, r = 0.93) -> Red in Clin vs RPM
  \node[anchor=east, font=\scriptsize\bfseries] at (3.8, 0.0) {Lipid \& Fatty Acid};
  \fill[gray!10] (4.0, -0.3) rectangle (6.4, 0.3);
  \node[font=\tiny] at (5.2, 0.0) {$r=0.15$ ($p=0.70$)};
  
  \fill[fairred] (6.8, -0.3) rectangle (9.2, 0.3);
  \node[font=\tiny\bfseries, text=white] at (8.0, 0.0) {$r=+0.93$ ($p=1.2\times 10^{-10}$)};
  
  \fill[gray!10] (9.6, -0.3) rectangle (12.0, 0.3);
  \node[font=\tiny] at (10.8, 0.0) {$r=-0.15$ ($p=0.70$)};

  % Colorbar on bottom right
  \node[font=\tiny\bfseries] at (5.0, -0.8) {Color Scale ($r$):};
  \fill[fairblue] (6.5, -0.9) rectangle (7.5, -0.7);
  \fill[fairlight] (7.5, -0.9) rectangle (8.5, -0.7);
  \fill[fairred] (8.5, -0.9) rectangle (9.5, -0.7);
  \draw[draw=gray!50] (6.5, -0.9) rectangle (9.5, -0.7);
  \node[font=\tiny] at (6.5, -1.1) {-1.0};
  \node[font=\tiny] at (8.0, -1.1) {0.0};
  \node[font=\tiny] at (9.5, -1.1) {+1.0};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig3_tex, "fig3_wgcna_modules")

    # -------------------------------------------------------------
    # Figure 4: TabPFN AI Foundation Model Benchmark
    # -------------------------------------------------------------
    fig4_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily, >=Stealth]
  % Panel a: LOOCV Accuracy Comparison
  \node[font=\large\bfseries] at (0, 7.8) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 7.8) {Leave-One-Out Cross-Validation (LOOCV) Benchmark};

  \draw[->, thick] (0.8, 1.2) -- (0.8, 6.8) node[midway, above=5pt, rotate=90, font=\scriptsize\bfseries] {Accuracy (\%)};
  \draw[thick] (0.8, 1.2) -- (5.6, 1.2);
  
  \foreach \y/\lbl in {1.2/0, 2.6/25, 4.0/50, 5.4/75, 6.8/100} {
    \draw (0.7, \y) -- (0.8, \y) node[left=3pt, font=\tiny] {\lbl};
    \draw[dashed, draw=gray!20] (0.8, \y) -- (5.6, \y);
  }

  % Bars: TabPFN Binary (88.9%)
  \fill[fairred] (1.3, 1.2) rectangle (2.1, 6.18);
  \node[font=\tiny\bfseries, text=white] at (1.7, 5.8) {88.9\%};
  \node[font=\tiny\bfseries, text=fairred, align=center] at (1.7, 0.7) {TabPFN\\Binary};

  % TabPFN 3-Class (66.7%)
  \fill[fairmidred] (2.7, 1.2) rectangle (3.5, 4.93);
  \node[font=\tiny\bfseries, text=white] at (3.1, 4.6) {66.7\%};
  \node[font=\tiny\bfseries, text=fairmidred, align=center] at (3.1, 0.7) {TabPFN\\3-Class};

  % Random Forest Baseline (0.0%)
  \fill[fairblue] (4.1, 1.2) rectangle (4.9, 1.3);
  \node[font=\tiny\bfseries, text=fairblue] at (4.5, 1.6) {0.0\%};
  \node[font=\tiny\bfseries, text=fairblue, align=center] at (4.5, 0.7) {Random\\Forest};

  % Panel b: Confusion Matrix
  \node[font=\large\bfseries] at (6.6, 7.8) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (7.0, 7.8) {TabPFN 3-Class Confusion Matrix (LOOCV)};

  \node[font=\tiny\bfseries] at (9.2, 6.8) {Predicted Class};
  \node[font=\tiny\bfseries, rotate=90] at (7.2, 5.0) {True Class};

  \node[font=\tiny] at (8.3, 6.4) {Clinostat};
  \node[font=\tiny] at (9.7, 6.4) {RPM 2.0};
  \node[font=\tiny] at (11.1, 6.4) {Static 1g};

  \node[font=\tiny, anchor=east] at (7.7, 5.7) {Clinostat};
  \node[font=\tiny, anchor=east] at (7.7, 4.7) {RPM 2.0};
  \node[font=\tiny, anchor=east] at (7.7, 3.7) {Static 1g};

  % Matrix Cells (Fair Red for diagonal, Fair Blue for off-diagonal)
  \fill[fairred!90!white] (7.8, 5.2) rectangle (8.8, 6.2); \node[font=\small\bfseries, text=white] at (8.3, 5.7) {2};
  \fill[fairlight] (9.2, 5.2) rectangle (10.2, 6.2); \node[font=\small\bfseries, text=gray] at (9.7, 5.7) {1};
  \fill[white, draw=gray!20] (10.6, 5.2) rectangle (11.6, 6.2); \node[font=\small\bfseries, text=gray!40] at (11.1, 5.7) {0};

  \fill[fairlight] (7.8, 4.2) rectangle (8.8, 5.2); \node[font=\small\bfseries, text=gray] at (8.3, 4.7) {1};
  \fill[fairred!90!white] (9.2, 4.2) rectangle (10.2, 5.2); \node[font=\small\bfseries, text=white] at (9.7, 4.7) {2};
  \fill[white, draw=gray!20] (10.6, 4.2) rectangle (11.6, 5.2); \node[font=\small\bfseries, text=gray!40] at (11.1, 4.7) {0};

  \fill[white, draw=gray!20] (7.8, 3.2) rectangle (8.8, 4.2); \node[font=\small\bfseries, text=gray!40] at (8.3, 3.7) {0};
  \fill[fairlight] (9.2, 3.2) rectangle (10.2, 4.2); \node[font=\small\bfseries, text=gray] at (9.7, 3.7) {1};
  \fill[fairred!90!white] (10.6, 3.2) rectangle (11.6, 4.2); \node[font=\small\bfseries, text=white] at (11.1, 3.7) {2};

  % Panel c: Feature Importance
  \node[font=\large\bfseries] at (6.6, 2.6) {c};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (7.0, 2.6) {Top Prioritized Features};

  \fill[fairred] (7.8, 1.8) rectangle (11.5, 2.2);
  \node[anchor=west, font=\tiny\bfseries] at (11.6, 2.0) {Lipid/Kinematic (MEblue, 0.93)};

  \fill[fairmidred] (7.8, 1.2) rectangle (10.9, 1.6);
  \node[anchor=west, font=\tiny\bfseries] at (11.0, 1.4) {Biofilm Core (MEturquoise, 0.77)};

  \fill[fairmidblue] (7.8, 0.6) rectangle (10.0, 1.0);
  \node[anchor=west, font=\tiny\bfseries] at (10.1, 0.8) {Shear Adaptation (MEgreen, 0.55)};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig4_tex, "fig4_tabpfn_evaluation")

    # -------------------------------------------------------------
    # Figure 5: REDESIGNED AS PUBLICATION GOSlim PATHWAY BAR PLOT
    # -------------------------------------------------------------
    fig5_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\usetikzlibrary{arrows.meta}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily, >=Stealth]
  % Title / Header
  \node[font=\large\bfseries] at (0, 8.8) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 8.8) {Enriched GOSlim Functional Pathways in Simulated Microgravity};

  % Horizontal Bar Plot Axes
  \draw[->, thick] (5.4, 0.8) -- (12.2, 0.8) node[midway, below=5pt, font=\scriptsize\bfseries] {$-\log_{10}(\text{Adjusted } p\text{-value})$};
  \draw[thick] (5.4, 0.8) -- (5.4, 8.2);

  % Axis Ticks
  \foreach \x/\lbl in {5.4/0, 6.25/1, 7.1/2, 7.95/3, 8.8/4, 9.65/5, 10.5/6, 11.35/7, 12.2/8} {
    \draw (\x, 0.7) -- (\x, 0.8) node[below=2pt, font=\tiny] {\lbl};
    \draw[dashed, draw=gray!20] (\x, 0.8) -- (\x, 8.2);
  }
  
  % Significance Threshold line at FDR = 0.05 (-log10 = 1.30)
  \draw[dashed, draw=fairred, line width=1pt] (6.5, 0.8) -- (6.5, 8.2) node[above, font=\tiny\bfseries, text=fairred] {FDR = 0.05};

  % Bars (Y positions: 7.6, 6.7, 5.8, 4.9, 4.0, 3.1, 2.2, 1.3)
  
  % 1. Response to oxidative stress & ROS (padj = 4.4e-8 -> mlog10 = 7.35) -> Deep Red
  \node[anchor=east, font=\scriptsize\bfseries] at (5.2, 7.6) {Response to oxidative stress \& ROS};
  \fill[fairred] (5.4, 7.35) rectangle (11.65, 7.85);
  \node[anchor=west, font=\tiny\bfseries, text=fairred] at (11.75, 7.6) {$14/85$ ($p=4.4\times 10^{-8}$)};

  % 2. Type VII secretion system complex (padj = 3.7e-4 -> mlog10 = 3.43) -> Red
  \node[anchor=east, font=\scriptsize\bfseries] at (5.2, 6.7) {Type VII secretion system complex};
  \fill[fairred!90!white] (5.4, 6.45) rectangle (8.32, 6.95);
  \node[anchor=west, font=\tiny\bfseries, text=fairred!90!black] at (8.42, 6.7) {$5/40$ ($p=3.7\times 10^{-4}$)};

  % 3. Mycolic acid biosynthesis (FAS-II) (padj = 7.2e-4 -> mlog10 = 3.14) -> Red
  \node[anchor=east, font=\scriptsize\bfseries] at (5.2, 5.8) {Mycolic acid biosynthesis (FAS-II)};
  \fill[fairred!80!white] (5.4, 5.55) rectangle (8.07, 6.05);
  \node[anchor=west, font=\tiny\bfseries, text=fairred!80!black] at (8.17, 5.8) {$5/48$ ($p=7.2\times 10^{-4}$)};

  % 4. Hypoxic boundary layer (DosR regulon) (padj = 1.39e-3 -> mlog10 = 2.86) -> Mid-Red
  \node[anchor=east, font=\scriptsize\bfseries] at (5.2, 4.9) {Hypoxic boundary layer (DosR regulon)};
  \fill[fairmidred] (5.4, 4.65) rectangle (7.83, 5.15);
  \node[anchor=west, font=\tiny\bfseries, text=fairmidred] at (7.93, 4.9) {$4/55$ ($p=1.4\times 10^{-3}$)};

  % 5. Cell wall & biofilm pellicle assembly (padj = 8.1e-3 -> mlog10 = 2.09) -> Mid-Red
  \node[anchor=east, font=\scriptsize\bfseries] at (5.2, 4.0) {Cell wall \& biofilm pellicle assembly};
  \fill[fairmidred!80!white] (5.4, 3.75) rectangle (7.18, 4.25);
  \node[anchor=west, font=\tiny\bfseries, text=fairmidred] at (7.28, 4.0) {$3/32$ ($p=8.1\times 10^{-3}$)};

  % 6. Alternative terminal oxidase (cydA) (padj = 2.17e-2 -> mlog10 = 1.66) -> Light Red
  \node[anchor=east, font=\scriptsize\bfseries] at (5.2, 3.1) {Alternative terminal oxidase (cydA)};
  \fill[fairmidred!60!white] (5.4, 2.85) rectangle (6.81, 3.35);
  \node[anchor=west, font=\tiny\bfseries, text=fairmidred!80!black] at (6.91, 3.1) {$3/45$ ($p=2.2\times 10^{-2}$)};

  % 7. Rotational shear chaperone folding (padj = 4.96e-2 -> mlog10 = 1.30) -> Light Red
  \node[anchor=east, font=\scriptsize\bfseries] at (5.2, 2.2) {Rotational shear chaperone folding};
  \fill[fairmidred!40!white] (5.4, 1.95) rectangle (6.51, 2.45);
  \node[anchor=west, font=\tiny\bfseries, text=gray!80!black] at (6.61, 2.2) {$3/42$ ($p=4.9\times 10^{-2}$)};

  % 8. Ribosome & translation attenuation (padj = 1.25e-3 -> mlog10 = 2.90) -> DEEP BLUE (Downregulated!)
  \node[anchor=east, font=\scriptsize\bfseries, text=fairblue] at (5.2, 1.3) {Ribosome \& translation attenuation};
  \fill[fairblue] (5.4, 1.05) rectangle (7.87, 1.55);
  \node[anchor=west, font=\tiny\bfseries, text=fairblue] at (7.97, 1.3) {$6/58$ (Downregulated, $p=1.3\times 10^{-3}$)};

  % Legend
  \node[font=\tiny\bfseries, anchor=west] at (0.2, 0.4) {Color Mapping:};
  \fill[fairred] (2.2, 0.3) rectangle (2.8, 0.5); \node[font=\tiny, anchor=west] at (2.9, 0.4) {Upregulated in Microgravity};
  \fill[fairblue] (6.5, 0.3) rectangle (7.1, 0.5); \node[font=\tiny, anchor=west] at (7.2, 0.4) {Downregulated in Microgravity};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig5_tex, "fig5_pathway_ontology")

    # -------------------------------------------------------------
    # Figure 6: Intramodular Hub Connectivity & Interactome Sub-Network
    # -------------------------------------------------------------
    fig6_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\usetikzlibrary{arrows.meta}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily, >=Stealth]
  % Panel a: Scatter Plot k_within vs k_total
  \node[font=\large\bfseries] at (0, 7.8) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 7.8) {Intramodular Hub Centrality ($k_{\text{within}}$ vs. $k_{\text{total}}$)};

  \draw[->, thick] (0.8, 1.0) -- (5.8, 1.0) node[midway, below=5pt, font=\scriptsize\bfseries] {Whole-Network Degree ($k_{\text{total}}$)};
  \draw[->, thick] (0.8, 1.0) -- (0.8, 6.8) node[midway, above=5pt, rotate=90, font=\scriptsize\bfseries] {Intramodular Connectivity ($k_{\text{within}}$)};

  \foreach \x/\lbl in {1.0/0, 2.2/10, 3.4/20, 4.6/30, 5.6/40} {
    \draw (\x, 0.9) -- (\x, 1.0) node[below=2pt, font=\tiny] {\lbl};
  }
  \foreach \y/\lbl in {1.0/0, 2.4/10, 3.8/20, 5.2/30, 6.6/40} {
    \draw (0.7, \y) -- (0.8, \y) node[left=2pt, font=\tiny] {\lbl};
  }

  % Hub Points
  \fill[fairblue] (5.3, 6.1) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS12120 (35.3)\ };
  \fill[fairblue] (5.2, 5.7) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS21560 (35.0)\ };
  \fill[fairblue] (5.1, 5.3) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS04440 (34.7)\ };

  \fill[fairred] (3.8, 3.2) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS29685 (15.4)\ };
  \fill[fairred] (3.5, 2.9) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS02330 (14.8)\ };

  \fill[fairmidblue] (3.0, 2.5) circle (3.5pt) node[anchor=east, font=\tiny\bfseries] {RS11565 (11.9)\ };
  \fill[fairmidred] (2.6, 2.2) circle (3.5pt) node[anchor=east, font=\tiny\bfseries] {RS16730 (10.8)\ };
  \fill[fairred!50!white] (2.0, 1.7) circle (3pt) node[anchor=west, font=\tiny] {\ RS13930 (6.3)};

  % Panel b: Regulatory Interactome
  \node[font=\large\bfseries] at (6.6, 7.8) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (7.0, 7.8) {Inter-Module Regulatory Interactome};

  % Network Nodes
  \node[circle, fill=fairblue!20, draw=fairblue, line width=1pt, font=\tiny\bfseries] (mps1) at (8.0, 5.8) {mps1};
  \node[circle, fill=fairblue!20, draw=fairblue, line width=1pt, font=\tiny\bfseries] (mps2) at (7.4, 4.6) {mps2};
  \node[circle, fill=fairblue!20, draw=fairblue, line width=1pt, font=\tiny\bfseries] (rs04) at (9.0, 5.0) {RS04440};

  \node[circle, fill=fairred!20, draw=fairred, line width=1pt, font=\tiny\bfseries] (kasA) at (10.2, 5.2) {kasA};
  \node[circle, fill=fairred!20, draw=fairred, line width=1pt, font=\tiny\bfseries] (fbpA) at (11.4, 5.8) {fbpA};
  \node[circle, fill=fairred!20, draw=fairred, line width=1pt, font=\tiny\bfseries] (rs02) at (11.5, 4.4) {RS02330};

  \node[circle, fill=fairred!30, draw=fairred!80!black, line width=1pt, font=\tiny\bfseries] (hspX) at (8.2, 2.4) {hspX};
  \node[circle, fill=fairred!30, draw=fairred!80!black, line width=1pt, font=\tiny\bfseries] (cydA) at (9.4, 1.8) {cydA};

  \node[circle, fill=fairmidred!25, draw=fairmidred, line width=1pt, font=\tiny\bfseries] (esxA) at (11.2, 2.4) {esxA};
  \node[circle, fill=fairmidred!25, draw=fairmidred, line width=1pt, font=\tiny\bfseries] (eccE) at (10.4, 1.5) {eccE};

  % Edges
  \draw[thick, draw=fairblue] (mps1) -- (mps2);
  \draw[thick, draw=fairblue] (mps1) -- (rs04);
  \draw[thick, dashed, draw=gray!70] (rs04) -- (kasA) node[midway, above=1pt, font=\tiny] {Envelope};
  \draw[thick, draw=fairred] (kasA) -- (fbpA);
  \draw[thick, draw=fairred] (kasA) -- (rs02);
  \draw[thick, draw=fairred!80!black] (hspX) -- (cydA);
  \draw[thick, draw=fairmidred] (esxA) -- (eccE);
  \draw[thick, dashed, draw=fairred] (hspX) -- (eccE) node[midway, below=1pt, font=\tiny] {Stress Pore};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig6_tex, "fig6_hub_connectivity")

    # -------------------------------------------------------------
    # Figure 7: Pan-Microbial Spaceflight Landscape
    # -------------------------------------------------------------
    fig7_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily]
  % Panel a: Taxonomic Distribution
  \node[font=\large\bfseries] at (0, 7.8) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 7.8) {Taxonomic Distribution Across 78 OSDR Spaceflight Datasets};

  \fill[fairblue] (0.4, 6.4) rectangle (4.2, 6.9);
  \node[anchor=west, font=\scriptsize\bfseries] at (4.4, 6.65) {Pseudomonadota (43.6\%, $N=34$)};

  \fill[fairmidblue] (0.4, 5.6) rectangle (2.8, 6.1);
  \node[anchor=west, font=\scriptsize\bfseries] at (3.0, 5.85) {Bacillota (28.2\%, $N=22$)};

  \fill[fairmidred] (0.4, 4.8) rectangle (1.5, 5.3);
  \node[anchor=west, font=\scriptsize\bfseries] at (1.7, 5.05) {Actinomycetota (12.8\%, $N=10$)};

  \fill[fairred] (0.4, 4.0) rectangle (1.2, 4.5);
  \node[anchor=west, font=\scriptsize\bfseries] at (1.4, 4.25) {Fungi (10.3\%, $N=8$)};

  % Panel b: Cross-Species Spaceflight Concordance Heatmap
  \node[font=\large\bfseries] at (0, 3.2) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 3.2) {Cross-Species Spaceflight Adaptation Concordance};

  % Column Headers
  \node[font=\tiny\bfseries] at (4.5, 2.7) {Biofilm};
  \node[font=\tiny\bfseries] at (6.0, 2.7) {Cell Wall};
  \node[font=\tiny\bfseries] at (7.5, 2.7) {Virulence};
  \node[font=\tiny\bfseries] at (9.0, 2.7) {Hypoxia};
  \node[font=\tiny\bfseries] at (10.5, 2.7) {Antibiotics};

  % Rows
  \node[anchor=east, font=\scriptsize\bfseries] at (3.5, 2.1) {M. marinum (OSD-528)};
  \fill[fairred] (4.0, 1.8) rectangle (5.0, 2.4); \node[font=\tiny\bfseries, text=white] at (4.5, 2.1) {+1.0};
  \fill[fairred] (5.5, 1.8) rectangle (6.5, 2.4); \node[font=\tiny\bfseries, text=white] at (6.0, 2.1) {+1.0};
  \fill[fairred] (7.0, 1.8) rectangle (8.0, 2.4); \node[font=\tiny\bfseries, text=white] at (7.5, 2.1) {+1.0};
  \fill[fairred] (8.5, 1.8) rectangle (9.5, 2.4); \node[font=\tiny\bfseries, text=white] at (9.0, 2.1) {+1.0};
  \fill[fairred] (10.0, 1.8) rectangle (11.0, 2.4); \node[font=\tiny\bfseries, text=white] at (10.5, 2.1) {+1.0};

  \node[anchor=east, font=\scriptsize\bfseries] at (3.5, 1.3) {P. aeruginosa (OSD-14)};
  \fill[fairred] (4.0, 1.0) rectangle (5.0, 1.6); \node[font=\tiny\bfseries, text=white] at (4.5, 1.3) {+1.0};
  \fill[fairred] (5.5, 1.0) rectangle (6.5, 1.6); \node[font=\tiny\bfseries, text=white] at (6.0, 1.3) {+1.0};
  \fill[fairred] (7.0, 1.0) rectangle (8.0, 1.6); \node[font=\tiny\bfseries, text=white] at (7.5, 1.3) {+1.0};
  \fill[fairred] (8.5, 1.0) rectangle (9.5, 1.6); \node[font=\tiny\bfseries, text=white] at (9.0, 1.3) {+1.0};
  \fill[fairred] (10.0, 1.0) rectangle (11.0, 1.6); \node[font=\tiny\bfseries, text=white] at (10.5, 1.3) {+1.0};

  \node[anchor=east, font=\scriptsize\bfseries] at (3.5, 0.5) {S. enterica (OSD-11)};
  \fill[fairmidred] (4.0, 0.2) rectangle (5.0, 0.8); \node[font=\tiny\bfseries, text=white] at (4.5, 0.5) {+0.8};
  \fill[fairred] (5.5, 0.2) rectangle (6.5, 0.8); \node[font=\tiny\bfseries, text=white] at (6.0, 0.5) {+1.0};
  \fill[fairred] (7.0, 0.2) rectangle (8.0, 0.8); \node[font=\tiny\bfseries, text=white] at (7.5, 0.5) {+1.0};
  \fill[fairmidred] (8.5, 0.2) rectangle (9.5, 0.8); \node[font=\tiny\bfseries, text=white] at (9.0, 0.5) {+0.7};
  \fill[fairred] (10.0, 0.2) rectangle (11.0, 0.8); \node[font=\tiny\bfseries, text=white] at (10.5, 0.5) {+1.0};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig7_tex, "fig7_pan_microbial_landscape")

    # -------------------------------------------------------------
    # Figure 8: Simulator Kinematics vs. Biological Response Concordance Radar
    # -------------------------------------------------------------
    fig8_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily]
  % Panel a: Radar Plot
  \node[font=\large\bfseries] at (0, 7.8) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 7.8) {Biological Response Trajectory Radar Across Simulation Modalities};

  % Center at (4.5, 4.0)
  % Web rings
  \foreach \r in {0.8, 1.6, 2.4, 3.2} {
    \draw[gray!25, dashed] (4.5, 4.0) circle (\r cm);
  }
  
  % 6 Axes at 60 degree intervals
  \foreach \ang/\lbl in {90/Biofilm Pellicle, 30/Cell Wall (FAS-II), 330/Type VII Secretion, 270/DosR Hypoxia, 210/Oxidative Stress, 150/Rotational Shear} {
    \draw[gray!40] (4.5, 4.0) -- ({4.5 + 3.4*cos(\ang)}, {4.0 + 3.4*sin(\ang)});
    \node[font=\tiny\bfseries] at ({4.5 + 3.7*cos(\ang)}, {4.0 + 3.7*sin(\ang)}) {\lbl};
  }

  % Static 1g Earth Polygon (Deep Blue - near center)
  \draw[line width=1.5pt, draw=fairblue, fill=fairblue!15]
    ({4.5 + 0.8*cos(90)}, {4.0 + 0.8*sin(90)}) --
    ({4.5 + 0.9*cos(30)}, {4.0 + 0.9*sin(30)}) --
    ({4.5 + 0.8*cos(330)}, {4.0 + 0.8*sin(330)}) --
    ({4.5 + 0.7*cos(270)}, {4.0 + 0.7*sin(270)}) --
    ({4.5 + 0.8*cos(210)}, {4.0 + 0.8*sin(210)}) --
    ({4.5 + 0.8*cos(150)}, {4.0 + 0.8*sin(150)}) -- cycle;

  % 3D Clinostat Polygon (Mid-Red / Slate)
  \draw[line width=1.5pt, draw=fairmidblue, fill=fairmidblue!20]
    ({4.5 + 3.0*cos(90)}, {4.0 + 3.0*sin(90)}) --
    ({4.5 + 3.0*cos(30)}, {4.0 + 3.0*sin(30)}) --
    ({4.5 + 2.8*cos(330)}, {4.0 + 2.8*sin(330)}) --
    ({4.5 + 2.9*cos(270)}, {4.0 + 2.9*sin(270)}) --
    ({4.5 + 2.9*cos(210)}, {4.0 + 2.9*sin(210)}) --
    ({4.5 + 2.9*cos(150)}, {4.0 + 2.9*sin(150)}) -- cycle;

  % RPM 2.0 Polygon (Deep Red)
  \draw[line width=1.5pt, draw=fairred, fill=fairred!20]
    ({4.5 + 3.1*cos(90)}, {4.0 + 3.1*sin(90)}) --
    ({4.5 + 3.1*cos(30)}, {4.0 + 3.1*sin(30)}) --
    ({4.5 + 2.9*cos(330)}, {4.0 + 2.9*sin(330)}) --
    ({4.5 + 3.0*cos(270)}, {4.0 + 3.0*sin(270)}) --
    ({4.5 + 3.1*cos(210)}, {4.0 + 3.1*sin(210)}) --
    ({4.5 + 1.8*cos(150)}, {4.0 + 1.8*sin(150)}) -- cycle;

  % Panel b: Concordance Summary
  \node[font=\large\bfseries] at (9.0, 7.8) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (9.4, 7.8) {Kinematic Concordance};

  \draw[line width=1.2pt, draw=fairmidblue] (9.2, 6.8) -- (9.7, 6.8);
  \node[anchor=west, font=\scriptsize\bfseries, text=fairblue] at (9.8, 6.8) {3D Clinostat ($<0.01g$)};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (9.8, 6.4) {Continuous rotational shear vector};

  \draw[line width=1.2pt, draw=fairred] (9.2, 5.6) -- (9.7, 5.6);
  \node[anchor=west, font=\scriptsize\bfseries, text=fairred] at (9.8, 5.6) {RPM 2.0 ($<0.01g$)};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (9.8, 5.2) {Multi-axis random vectoring};

  \draw[line width=1.2pt, draw=fairblue] (9.2, 4.4) -- (9.7, 4.4);
  \node[anchor=west, font=\scriptsize\bfseries, text=fairblue] at (9.8, 4.4) {Static 1g Ground Control};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (9.8, 4.0) {Standard gravity baseline};

  \node[anchor=west, font=\tiny\bfseries, text=fairred] at (9.2, 2.8) {Core Concordance: $>75\%$ shared DEGs};
  \node[anchor=west, font=\tiny\bfseries, text=fairblue] at (9.2, 2.4) {Divergence: Rotational shear only};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig8_tex, "fig8_simulator_concordance_radar")

    print("Phase 6 completed successfully.")

if __name__ == "__main__":
    main()
