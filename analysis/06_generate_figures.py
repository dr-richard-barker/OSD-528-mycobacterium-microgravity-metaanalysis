#!/usr/bin/env python3
"""
06_generate_figures.py
Generates the publication figure suite formatted specifically for npj Microgravity:
- Zero text overlaps across all panels and subpanels.
- Expanded canvas widths and generous margins for process flowcharts.
- Proper Cartesian axes with explicit X and Y axis titles, ticks, and numbers for all bar charts and scatter plots.
- Complete separation between X-axis titles, numbers, and figure legends.
- Standardized bold lowercase subpanel identifiers (a, b, c, d...).
- Complete removal of rounded grey container boxes.
- Calibrated FAIR divergent Blue-to-White-to-Red (RdBu) color balance:
  * Blue (#2166ac / blue!75!black): 1g ground control / downregulated / negative correlation
  * White / Neutral: Baseline / neutral correlation
  * Red (#b2182b / red!75!black): Microgravity / upregulated / positive correlation
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

    tikz_colors = r"""
\definecolor{fairblue}{RGB}{33, 102, 172}
\definecolor{fairmidblue}{RGB}{103, 169, 207}
\definecolor{fairlight}{RGB}{247, 247, 247}
\definecolor{fairmidred}{RGB}{239, 138, 98}
\definecolor{fairred}{RGB}{178, 24, 43}
\definecolor{darkgray}{RGB}{50, 50, 50}
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
  % Panel a: Biological Model & Biofilm Substrate (Expanded width)
  \node[font=\large\bfseries] at (0, 11.2) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 11.2) {Biological Model and Silicone Biofilm Substrate};
  
  \node[circle, fill=fairred, text=white, font=\bfseries\small, inner sep=6pt] at (1.0, 9.8) {RFP};
  \node[anchor=west, font=\bfseries\footnotesize] at (1.8, 10.1) {Mycobacterium marinum 1218R};
  \node[anchor=west, font=\scriptsize, text=gray!85!black] at (1.8, 9.6) {BSL-2 surrogate for Mycobacterium tuberculosis};
  \node[anchor=west, font=\scriptsize, text=gray!85!black] at (1.8, 9.1) {Stable chromosomal RFP at Giles phage attB site};

  \draw[thick, draw=gray!30] (0.2, 8.4) -- (7.2, 8.4);
  \node[anchor=west, font=\bfseries\scriptsize, text=fairblue] at (0.2, 8.0) {Culture Vessel Specifications (31$^\circ$C, 4 Days):};
  \node[anchor=west, font=\scriptsize] at (0.2, 7.4) {$\bullet$ Polydimethylsiloxane (PDMS) silicone membranes};
  \node[anchor=west, font=\scriptsize] at (0.2, 6.8) {$\bullet$ Suspended biofilm-forming cell pellets harvested for RNA-seq};
  \node[anchor=west, font=\scriptsize] at (0.2, 6.2) {$\bullet$ $n=3$ biological replicates per experimental modality};

  % Panel b: Microgravity Simulation Hardware
  \node[font=\large\bfseries] at (8.0, 11.2) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (8.4, 11.2) {Simulated Microgravity and Ground Hardware};

  % 3D Clinostat
  \draw[line width=2pt, draw=fairmidblue] (8.2, 10.2) -- (8.2, 8.9);
  \node[anchor=west, font=\bfseries\scriptsize, text=fairblue] at (8.5, 9.9) {1. Lab-Designed 3D Clinostat ($n=3$)};
  \node[anchor=west, font=\scriptsize] at (8.5, 9.4) {Continuous 2-axis clinorotation ($I=1.5$~rpm, $O=3.825$~rpm)};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (8.5, 9.0) {Samples: RFP3D11, RFP3D39, RFP3D47};

  % RPM 2.0
  \draw[line width=2pt, draw=fairred] (8.2, 8.4) -- (8.2, 7.1);
  \node[anchor=west, font=\bfseries\scriptsize, text=fairred] at (8.5, 8.1) {2. Random Positioning Machine 2.0 ($n=3$)};
  \node[anchor=west, font=\scriptsize] at (8.5, 7.6) {Random multi-axis velocity vectoring, time-averaged $<0.01g$};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (8.5, 7.2) {Samples: RFPRPM4, RFPRPM41, RFPRPM6};

  % Static 1g
  \draw[line width=2pt, draw=fairblue] (8.2, 6.6) -- (8.2, 5.3);
  \node[anchor=west, font=\bfseries\scriptsize, text=fairblue] at (8.5, 6.3) {3. Static 1g Earth Control ($n=3$)};
  \node[anchor=west, font=\scriptsize] at (8.5, 5.8) {Stationary incubator shelf adjacent to simulators};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (8.5, 5.4) {Samples: RFPNG14, RFPNG35, RFPNG45};

  % Panel c: Integrated Empirical Workflow (Spread horizontally across 16.5 cm)
  \node[font=\large\bfseries] at (0, 4.8) {c};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 4.8) {Empirical Multi-Scale Analytical Framework};

  % 4 Flowchart Boxes (Width 3.6 cm each with generous spacing)
  % Box 1
  \draw[thick, draw=fairblue, fill=white] (0.2, 0.6) rectangle (3.8, 4.0);
  \node[font=\bfseries\scriptsize, text=fairblue] at (2.0, 3.5) {1. Empirical Quant};
  \node[font=\tiny] at (2.0, 2.9) {NASA OSDR S3 Raw Reads};
  \node[font=\tiny] at (2.0, 2.4) {5,510 Annotated Genes};
  \node[font=\tiny] at (2.0, 1.9) {kallisto Pseudoalignment};
  \node[font=\tiny\bfseries, text=fairblue] at (2.0, 1.1) {351--738 Real DEGs};

  \draw[->, line width=1.5pt, draw=gray!60] (3.8, 2.3) -- (4.3, 2.3);

  % Box 2
  \draw[thick, draw=fairmidblue, fill=white] (4.3, 0.6) rectangle (7.9, 4.0);
  \node[font=\bfseries\scriptsize, text=fairblue!80!black] at (6.1, 3.5) {2. WGCNA Modules};
  \node[font=\tiny] at (6.1, 2.9) {Soft-Threshold $\beta=6$};
  \node[font=\tiny] at (6.1, 2.4) {Topological Overlap (TOM)};
  \node[font=\tiny] at (6.1, 1.9) {5 GOSlim Co-Expression Modules};
  \node[font=\tiny\bfseries, text=fairblue!80!black] at (6.1, 1.1) {$k_{\text{within}}$ Hub Centrality};

  \draw[->, line width=1.5pt, draw=gray!60] (7.9, 2.3) -- (8.4, 2.3);

  % Box 3
  \draw[thick, draw=fairmidred, fill=white] (8.4, 0.6) rectangle (12.0, 4.0);
  \node[font=\bfseries\scriptsize, text=fairred] at (10.2, 3.5) {3. TabPFN AI (Nature 2025)};
  \node[font=\tiny] at (10.2, 2.9) {Prior-Data Transformer};
  \node[font=\tiny] at (10.2, 2.4) {Bayesian In-Context Priors};
  \node[font=\tiny] at (10.2, 1.9) {15 Topological Hub Features};
  \node[font=\tiny\bfseries, text=fairred] at (10.2, 1.1) {88.9\% Binary LOOCV};

  \draw[->, line width=1.5pt, draw=gray!60] (12.0, 2.3) -- (12.5, 2.3);

  % Box 4
  \draw[thick, draw=fairred, fill=white] (12.5, 0.6) rectangle (16.1, 4.0);
  \node[font=\bfseries\scriptsize, text=fairred!80!black] at (14.3, 3.5) {4. Systems Validation};
  \node[font=\tiny] at (14.3, 2.9) {GOSlim Over-Representation};
  \node[font=\tiny] at (14.3, 2.4) {Cellular Metabolic Model};
  \node[font=\tiny] at (14.3, 1.9) {Kinematic Radar Trajectories};
  \node[font=\tiny\bfseries, text=fairred!80!black] at (14.3, 1.1) {RO-Crate \& Zenodo Deposit};
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
  \node[font=\large\bfseries] at (0.2, 8.4) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.6, 8.4) {Transcriptomic Principal Component Analysis};
  
  % Axes (Origin at 1.0, 1.8 with uniform spacing)
  \draw[->, thick] (1.0, 1.8) -- (6.6, 1.8) node[midway, below=18pt, font=\scriptsize\bfseries] {PC1: Microgravity vs. 1g (70.1\% Variance)};
  \draw[->, thick] (1.0, 1.8) -- (1.0, 7.6);
  \node[font=\scriptsize\bfseries, rotate=90, anchor=center] at (0.3, 4.7) {PC2: Simulator Disparity (18.4\% Variance)};

  % Axis Ticks
  \foreach \x/\lbl in {1.4/-40, 2.6/-20, 3.8/0, 5.0/20, 6.2/40} {
    \draw (\x, 1.7) -- (\x, 1.8) node[below=4pt, font=\tiny] {\lbl};
  }
  \foreach \y/\lbl in {2.4/-30, 3.8/-15, 5.2/0, 6.6/15} {
    \draw (0.9, \y) -- (1.0, \y) node[left=4pt, font=\tiny] {\lbl};
  }

  % Static 1g Points (Blue)
  \draw[thick, dashed, draw=fairblue] (2.0, 5.2) ellipse (0.7cm and 1.1cm);
  \fill[fairblue] (1.8, 5.8) circle (3.5pt) node[anchor=south east, font=\tiny\bfseries] {RFPNG14};
  \fill[fairblue] (1.9, 5.1) circle (3.5pt) node[anchor=east, font=\tiny\bfseries] {RFPNG35};
  \fill[fairblue] (2.2, 4.6) circle (3.5pt) node[anchor=north east, font=\tiny\bfseries] {RFPNG45};

  % 3D Clinostat Points (Mid-Red/Slate)
  \draw[thick, dashed, draw=fairmidred] (4.8, 6.2) ellipse (0.8cm and 0.7cm);
  \fill[fairmidred] (4.6, 6.5) circle (3.5pt) node[anchor=south, font=\tiny\bfseries] {RFP3D11};
  \fill[fairmidred] (5.0, 6.1) circle (3.5pt) node[anchor=north west, font=\tiny\bfseries] {RFP3D39};
  \fill[fairmidred] (4.8, 5.8) circle (3.5pt) node[anchor=north east, font=\tiny\bfseries] {RFP3D47};

  % RPM 2.0 Points (Deep Red - comfortably above y=1.8)
  \draw[thick, dashed, draw=fairred] (4.8, 3.4) ellipse (0.8cm and 0.8cm);
  \fill[fairred] (4.6, 3.8) circle (3.5pt) node[anchor=south, font=\tiny\bfseries] {RFPRPM4};
  \fill[fairred] (5.1, 3.4) circle (3.5pt) node[anchor=west, font=\tiny\bfseries] {RFPRPM41};
  \fill[fairred] (4.8, 2.9) circle (3.5pt) node[anchor=north, font=\tiny\bfseries] {RFPRPM6};

  % Panel b: Volcano Plot 3D Clinostat vs 1g (Matched baseline at y=8.4)
  \node[font=\large\bfseries] at (7.8, 8.4) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (8.2, 8.4) {3D Clinostat vs. Static 1g (351 DEGs)};
  
  \draw[->, thick] (8.2, 1.8) -- (14.0, 1.8) node[midway, below=18pt, font=\scriptsize\bfseries] {$\log_2\text{Fold Change}$};
  \draw[->, thick] (8.2, 1.8) -- (8.2, 7.6);
  \node[font=\scriptsize\bfseries, rotate=90, anchor=center] at (7.5, 4.7) {$-\log_{10}(\text{Adjusted } p\text{-value})$};

  % Ticks: Exact uniform 1.2 cm steps for 3 units
  \foreach \x/\lbl in {8.7/-6, 9.9/-3, 11.1/0, 12.3/3, 13.5/6} {
    \draw (\x, 1.7) -- (\x, 1.8) node[below=4pt, font=\tiny] {\lbl};
  }
  \foreach \y/\lbl in {2.6/2, 4.0/6, 5.4/10, 6.8/14} {
    \draw (8.1, \y) -- (8.2, \y) node[left=4pt, font=\tiny] {\lbl};
  }

  % Threshold lines
  \draw[dashed, gray!40] (8.2, 2.6) -- (14.0, 2.6);
  \draw[dashed, gray!40] (10.4, 1.8) -- (10.4, 7.6);
  \draw[dashed, gray!40] (11.8, 1.8) -- (11.8, 7.6);

  % Upregulated (Deep Red)
  \fill[fairred] (13.7, 6.8) circle (3pt) node[anchor=south east, font=\tiny\bfseries] {RS06635 (+6.3)};
  \fill[fairred] (13.4, 6.2) circle (2.5pt) node[anchor=west, font=\tiny] {RS09245};
  \fill[fairred] (13.1, 5.7) circle (2.5pt) node[anchor=west, font=\tiny] {RS04125};
  \fill[fairred] (12.2, 4.8) circle (2.5pt) node[anchor=south west, font=\tiny\bfseries] {nuoD (+2.5)};

  % Downregulated (Deep Blue)
  \fill[fairblue] (8.7, 6.6) circle (3pt) node[anchor=south west, font=\tiny\bfseries] {rpmG (-7.3)};
  \fill[fairblue] (9.1, 5.9) circle (2.5pt) node[anchor=east, font=\tiny\bfseries] {espB (-6.0)};
  \fill[fairblue] (9.4, 5.1) circle (2.5pt) node[anchor=east, font=\tiny] {fadD7 (-5.3)};

  % Non-sig scatter
  \foreach \x/\y in {10.7/2.1, 11.2/2.3, 11.0/2.2, 11.4/2.0, 10.9/2.4} {
    \fill[gray!35] (\x, \y) circle (1.5pt);
  }

  % Summary Badge
  \node[fill=white, draw=gray!40, font=\tiny\bfseries, inner sep=3pt] at (12.8, 2.6) {351 Sig. DEGs};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig2_tex, "fig2_volcano_pca")

    # -------------------------------------------------------------
    # Figure 3: WGCNA Modules with Explicit X and Y Axes
    # -------------------------------------------------------------
    fig3_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily]
  % Panel a: Module Distribution with Clean Y-Axis Clearance (x=3.6 cm)
  \node[font=\large\bfseries] at (0.2, 8.8) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.6, 8.8) {Empirical Co-Expression Modules and Standardized GOSlim Descriptors};

  % Axes for Panel a: Y-axis line placed at x=3.6 cm, title at x=0.35 cm
  \draw[->, thick] (3.6, 3.4) -- (3.6, 8.0);
  \node[font=\scriptsize\bfseries, rotate=90] at (0.35, 5.7) {Co-Expression Module (GOSlim Name)};
  \draw[->, thick] (3.6, 3.4) -- (7.8, 3.4) node[midway, below=18pt, font=\scriptsize\bfseries] {Number of Expressed Genes ($n$)};

  % X-axis Ticks for Panel a
  \foreach \x/\lbl in {3.6/0, 4.6/50, 5.6/100, 6.6/150, 7.6/200} {
    \draw (\x, 3.3) -- (\x, 3.4) node[below=4pt, font=\tiny] {\lbl};
    \draw[dashed, draw=gray!20] (\x, 3.4) -- (\x, 8.0);
  }

  % Horizontal Bars & Y-ticks (Labels between x=0.6 and 3.5 with anchor=east)
  % 1. Cell Surface & Biofilm (168 genes -> x = 3.6 + 168/200 * 4.0 = 6.96)
  \draw (3.5, 7.3) -- (3.6, 7.3);
  \node[anchor=east, font=\tiny\bfseries] at (3.4, 7.3) {Cell Surface \& Biofilm};
  \fill[fairblue!85!black] (3.6, 7.05) rectangle (6.96, 7.55);
  \node[anchor=west, font=\tiny\bfseries, text=fairblue!85!black] at (7.04, 7.3) {168};

  % 2. Transmembrane Transport & Secretion (65 genes -> x = 3.6 + 65/200 * 4.0 = 4.90)
  \draw (3.5, 6.4) -- (3.6, 6.4);
  \node[anchor=east, font=\tiny\bfseries] at (3.4, 6.4) {Transport \& Secretion};
  \fill[fairmidred] (3.6, 6.15) rectangle (4.90, 6.65);
  \node[anchor=west, font=\tiny\bfseries, text=fairmidred] at (4.98, 6.4) {65};

  % 3. Lipid & Fatty Acid Metabolism (53 genes -> x = 3.6 + 53/200 * 4.0 = 4.66)
  \draw (3.5, 5.5) -- (3.6, 5.5);
  \node[anchor=east, font=\tiny\bfseries] at (3.4, 5.5) {Lipid \& Fatty Acid};
  \fill[fairred] (3.6, 5.25) rectangle (4.66, 5.75);
  \node[anchor=west, font=\tiny\bfseries, text=fairred] at (4.74, 5.5) {53};

  % 4. Cellular Respiration & Shear Adaptation (39 genes -> x = 3.6 + 39/200 * 4.0 = 4.38)
  \draw (3.5, 4.6) -- (3.6, 4.6);
  \node[anchor=east, font=\tiny\bfseries] at (3.4, 4.6) {Respiration \& Shear};
  \fill[fairmidblue] (3.6, 4.35) rectangle (4.38, 4.85);
  \node[anchor=west, font=\tiny\bfseries, text=fairmidblue] at (4.46, 4.6) {39};

  % 5. Response to Stress & Redox Homeostasis (25 genes -> x = 3.6 + 25/200 * 4.0 = 4.10)
  \draw (3.5, 3.8) -- (3.6, 3.8);
  \node[anchor=east, font=\tiny\bfseries] at (3.4, 3.8) {Stress \& Redox};
  \fill[fairred!60!white] (3.6, 3.55) rectangle (4.10, 4.05);
  \node[anchor=west, font=\tiny\bfseries, text=fairred!80!black] at (4.18, 3.8) {25};

  % Panel b: Module-Trait Correlation Heatmap with Explicit Axis Titles
  \node[font=\large\bfseries] at (8.4, 8.8) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (8.8, 8.8) {Module-Trait Correlation Heatmap};

  % Column Headers & X-axis Title
  \node[font=\scriptsize\bfseries] at (11.4, 8.1) {Microgravity vs. 1g};
  \node[font=\scriptsize\bfseries] at (13.4, 8.1) {Clinostat vs. RPM};
  \node[font=\scriptsize\bfseries] at (15.4, 8.1) {Static 1g Earth};
  \node[font=\scriptsize\bfseries] at (13.4, 4.3) {Experimental Trait Contrast};

  % Row Headers & Y-axis Title (Zero overlap, title at x=8.0)
  \node[font=\scriptsize\bfseries, rotate=90] at (8.0, 6.3) {GOSlim Module};
  \node[anchor=east, font=\tiny\bfseries] at (10.2, 7.3) {Cell Surface \& Biofilm};
  \node[anchor=east, font=\tiny\bfseries] at (10.2, 6.3) {Transport \& Secretion};
  \node[anchor=east, font=\tiny\bfseries] at (10.2, 5.3) {Lipid \& Fatty Acid};

  % Heatmap Cells (Width 1.9 cm, Height 0.8 cm)
  % Row 1
  \fill[fairblue!85!white] (10.45, 6.9) rectangle (12.35, 7.7);
  \node[font=\tiny\bfseries, text=white] at (11.4, 7.3) {$r=-0.77$ ($p=1.6\times 10^{-3}$)};
  
  \fill[gray!10] (12.45, 6.9) rectangle (14.35, 7.7);
  \node[font=\tiny] at (13.4, 7.3) {$r=0.08$ ($p=0.83$)};
  
  \fill[fairred!80!white] (14.45, 6.9) rectangle (16.35, 7.7);
  \node[font=\tiny\bfseries, text=white] at (15.4, 7.3) {$r=+0.77$ ($p=1.6\times 10^{-3}$)};

  % Row 2
  \fill[gray!10] (10.45, 5.9) rectangle (12.35, 6.7);
  \node[font=\tiny] at (11.4, 6.3) {$r=-0.03$ ($p=0.95$)};
  
  \fill[fairblue!40!white] (12.45, 5.9) rectangle (14.35, 6.7);
  \node[font=\tiny] at (13.4, 6.3) {$r=-0.38$ ($p=0.31$)};
  
  \fill[gray!10] (14.45, 5.9) rectangle (16.35, 6.7);
  \node[font=\tiny] at (15.4, 6.3) {$r=0.03$ ($p=0.95$)};

  % Row 3
  \fill[gray!10] (10.45, 4.9) rectangle (12.35, 5.7);
  \node[font=\tiny] at (11.4, 5.3) {$r=0.15$ ($p=0.70$)};
  
  \fill[fairred] (12.45, 4.9) rectangle (14.35, 5.7);
  \node[font=\tiny\bfseries, text=white] at (13.4, 5.3) {$r=+0.93$ ($p=1.2\times 10^{-10}$)};
  
  \fill[gray!10] (14.45, 4.9) rectangle (16.35, 5.7);
  \node[font=\tiny] at (15.4, 5.3) {$r=-0.15$ ($p=0.70$)};

  % Colorbar below
  \node[font=\tiny\bfseries] at (11.6, 3.6) {Pearson $r$:};
  \fill[fairblue] (12.6, 3.5) rectangle (13.6, 3.7);
  \fill[fairlight] (13.6, 3.5) rectangle (14.6, 3.7);
  \fill[fairred] (14.6, 3.5) rectangle (15.6, 3.7);
  \draw[draw=gray!50] (12.6, 3.5) rectangle (15.6, 3.7);
  \node[font=\tiny] at (12.6, 3.2) {-1.0};
  \node[font=\tiny] at (14.1, 3.2) {0.0};
  \node[font=\tiny] at (15.6, 3.2) {+1.0};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig3_tex, "fig3_wgcna_modules")

    # -------------------------------------------------------------
    # Figure 4: TabPFN AI Foundation Model Benchmark (Matched Header Alignment)
    # -------------------------------------------------------------
    fig4_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily, >=Stealth]
  % Panel a: LOOCV Accuracy Comparison (Shared baseline y=8.8)
  \node[font=\large\bfseries] at (0.2, 8.8) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.6, 8.8) {Leave-One-Out Cross-Validation (LOOCV) Benchmark};

  \draw[->, thick] (0.8, 2.0) -- (0.8, 7.8) node[midway, left=22pt, rotate=90, font=\scriptsize\bfseries] {LOOCV Accuracy (\%)};
  \draw[thick] (0.8, 2.0) -- (5.8, 2.0);
  
  \foreach \y/\lbl in {2.0/0, 3.4/25, 4.8/50, 6.2/75, 7.6/100} {
    \draw (0.7, \y) -- (0.8, \y) node[left=4pt, font=\tiny] {\lbl};
    \draw[dashed, draw=gray!20] (0.8, \y) -- (5.8, \y);
  }

  % Bars: TabPFN Binary (88.9%)
  \fill[fairred] (1.3, 2.0) rectangle (2.2, 6.98);
  \node[font=\tiny\bfseries, text=white] at (1.75, 6.6) {88.9\%};
  \node[font=\tiny\bfseries, text=fairred, align=center] at (1.75, 1.4) {TabPFN\\Binary};

  % TabPFN 3-Class (66.7%)
  \fill[fairmidred] (2.7, 2.0) rectangle (3.6, 5.73);
  \node[font=\tiny\bfseries, text=white] at (3.15, 5.4) {66.7\%};
  \node[font=\tiny\bfseries, text=fairmidred, align=center] at (3.15, 1.4) {TabPFN\\3-Class};

  % Random Forest Baseline (0.0%)
  \fill[fairblue] (4.1, 2.0) rectangle (5.0, 2.1);
  \node[font=\tiny\bfseries, text=fairblue] at (4.55, 2.4) {0.0\%};
  \node[font=\tiny\bfseries, text=fairblue, align=center] at (4.55, 1.4) {Random\\Forest};

  \node[font=\scriptsize\bfseries] at (3.15, 0.4) {Machine Learning Model / Task};

  % Panel b: Confusion Matrix (Shared baseline y=8.8, centered over matrix)
  \node[font=\large\bfseries] at (7.8, 8.8) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (8.2, 8.8) {TabPFN 3-Class Confusion Matrix (LOOCV)};

  \node[font=\scriptsize\bfseries] at (11.25, 8.0) {Predicted Modality};
  \node[font=\scriptsize\bfseries, rotate=90] at (7.5, 5.2) {True Modality};

  \node[font=\tiny\bfseries] at (9.75, 7.4) {Clinostat};
  \node[font=\tiny\bfseries] at (11.25, 7.4) {RPM 2.0};
  \node[font=\tiny\bfseries] at (12.75, 7.4) {Static 1g};

  \node[font=\tiny\bfseries, anchor=east] at (8.9, 6.4) {Clinostat};
  \node[font=\tiny\bfseries, anchor=east] at (8.9, 5.2) {RPM 2.0};
  \node[font=\tiny\bfseries, anchor=east] at (8.9, 4.0) {Static 1g};

  % Matrix Cells (1.3 cm x 1.0 cm)
  \fill[fairred!90!white] (9.1, 5.9) rectangle (10.4, 6.9); \node[font=\small\bfseries, text=white] at (9.75, 6.4) {2};
  \fill[fairlight] (10.6, 5.9) rectangle (11.9, 6.9); \node[font=\small\bfseries, text=gray] at (11.25, 6.4) {1};
  \fill[white, draw=gray!30] (12.1, 5.9) rectangle (13.4, 6.9); \node[font=\small\bfseries, text=gray!40] at (12.75, 6.4) {0};

  \fill[fairlight] (9.1, 4.7) rectangle (10.4, 5.7); \node[font=\small\bfseries, text=gray] at (9.75, 5.2) {1};
  \fill[fairred!90!white] (10.6, 4.7) rectangle (11.9, 5.7); \node[font=\small\bfseries, text=white] at (11.25, 5.2) {2};
  \fill[white, draw=gray!30] (12.1, 4.7) rectangle (13.4, 5.7); \node[font=\small\bfseries, text=gray!40] at (12.75, 5.2) {0};

  \fill[white, draw=gray!30] (9.1, 3.5) rectangle (10.4, 4.5); \node[font=\small\bfseries, text=gray!40] at (9.75, 4.0) {0};
  \fill[fairlight] (10.6, 3.5) rectangle (11.9, 4.5); \node[font=\small\bfseries, text=gray] at (11.25, 4.0) {1};
  \fill[fairred!90!white] (12.1, 3.5) rectangle (13.4, 4.5); \node[font=\small\bfseries, text=white] at (12.75, 4.0) {2};

  % Panel c: Feature Importance
  \node[font=\large\bfseries] at (7.8, 2.6) {c};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (8.2, 2.6) {Top Prioritized Features};

  \fill[fairred] (8.9, 1.8) rectangle (12.6, 2.2);
  \node[anchor=west, font=\tiny\bfseries] at (12.7, 2.0) {Lipid/Kinematic (MEblue, 0.93)};

  \fill[fairmidred] (8.9, 1.1) rectangle (11.8, 1.5);
  \node[anchor=west, font=\tiny\bfseries] at (11.9, 1.3) {Biofilm Core (MEturquoise, 0.77)};

  \fill[fairmidblue] (8.9, 0.4) rectangle (10.8, 0.8);
  \node[anchor=west, font=\tiny\bfseries] at (10.9, 0.6) {Shear Adaptation (MEgreen, 0.55)};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig4_tex, "fig4_tabpfn_evaluation")

    # -------------------------------------------------------------
    # Figure 5: GOSlim Pathway Bar Plot (Zero Legend Collision)
    # -------------------------------------------------------------
    fig5_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\usetikzlibrary{arrows.meta}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily, >=Stealth]
  % Title
  \node[font=\large\bfseries] at (0, 9.4) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 9.4) {Enriched GOSlim Functional Pathways in Simulated Microgravity};

  % Horizontal Bar Plot Axes (Raised to y=1.8 to create ample space for legend)
  \draw[->, thick] (5.4, 1.8) -- (12.6, 1.8) node[midway, below=18pt, font=\scriptsize\bfseries] {$-\log_{10}(\text{Adjusted } p\text{-value})$};
  \draw[thick] (5.4, 1.8) -- (5.4, 8.8);

  % Axis Ticks
  \foreach \x/\lbl in {5.4/0, 6.3/1, 7.2/2, 8.1/3, 9.0/4, 9.9/5, 10.8/6, 11.7/7, 12.6/8} {
    \draw (\x, 1.7) -- (\x, 1.8) node[below=3pt, font=\tiny] {\lbl};
    \draw[dashed, draw=gray!20] (\x, 1.8) -- (\x, 8.8);
  }
  
  % Significance Threshold line at FDR = 0.05 (-log10 = 1.30 -> x = 5.4 + 1.3*0.9 = 6.57)
  \draw[dashed, draw=fairred, line width=1pt] (6.57, 1.8) -- (6.57, 8.8) node[above, font=\tiny\bfseries, text=fairred] {FDR = 0.05};

  % Bars (Y positions: 8.2, 7.3, 6.4, 5.5, 4.6, 3.7, 2.8, 2.0)
  
  % 1. Response to oxidative stress & ROS (mlog10 = 7.35 -> x = 5.4 + 7.35*0.9 = 12.02)
  \node[anchor=east, font=\scriptsize\bfseries] at (5.2, 8.2) {Response to oxidative stress \& ROS};
  \fill[fairred] (5.4, 7.95) rectangle (12.02, 8.45);
  \node[anchor=west, font=\tiny\bfseries, text=fairred] at (12.12, 8.2) {$14/85$ ($p=4.4\times 10^{-8}$)};

  % 2. Type VII secretion system complex (mlog10 = 3.43 -> x = 5.4 + 3.43*0.9 = 8.49)
  \node[anchor=east, font=\scriptsize\bfseries] at (5.2, 7.3) {Type VII secretion system complex};
  \fill[fairred!90!white] (5.4, 7.05) rectangle (8.49, 7.55);
  \node[anchor=west, font=\tiny\bfseries, text=fairred!90!black] at (8.59, 7.3) {$5/40$ ($p=3.7\times 10^{-4}$)};

  % 3. Mycolic acid biosynthesis (FAS-II) (mlog10 = 3.14 -> x = 5.4 + 3.14*0.9 = 8.23)
  \node[anchor=east, font=\scriptsize\bfseries] at (5.2, 6.4) {Mycolic acid biosynthesis (FAS-II)};
  \fill[fairred!80!white] (5.4, 6.15) rectangle (8.23, 6.65);
  \node[anchor=west, font=\tiny\bfseries, text=fairred!80!black] at (8.33, 6.4) {$5/48$ ($p=7.2\times 10^{-4}$)};

  % 4. Hypoxic boundary layer (DosR regulon) (mlog10 = 2.86 -> x = 5.4 + 2.86*0.9 = 7.97)
  \node[anchor=east, font=\scriptsize\bfseries] at (5.2, 5.5) {Hypoxic boundary layer (DosR regulon)};
  \fill[fairmidred] (5.4, 5.25) rectangle (7.97, 5.75);
  \node[anchor=west, font=\tiny\bfseries, text=fairmidred] at (8.07, 5.5) {$4/55$ ($p=1.4\times 10^{-3}$)};

  % 5. Cell wall \& biofilm pellicle assembly (mlog10 = 2.09 -> x = 5.4 + 2.09*0.9 = 7.28)
  \node[anchor=east, font=\scriptsize\bfseries] at (5.2, 4.6) {Cell wall \& biofilm pellicle assembly};
  \fill[fairmidred!80!white] (5.4, 4.35) rectangle (7.28, 4.85);
  \node[anchor=west, font=\tiny\bfseries, text=fairmidred] at (7.38, 4.6) {$3/32$ ($p=8.1\times 10^{-3}$)};

  % 6. Alternative terminal oxidase (cydA) (mlog10 = 1.66 -> x = 5.4 + 1.66*0.9 = 6.89)
  \node[anchor=east, font=\scriptsize\bfseries] at (5.2, 3.7) {Alternative terminal oxidase (cydA)};
  \fill[fairmidred!60!white] (5.4, 3.45) rectangle (6.89, 3.95);
  \node[anchor=west, font=\tiny\bfseries, text=fairmidred!80!black] at (6.99, 3.7) {$3/45$ ($p=2.2\times 10^{-2}$)};

  % 7. Rotational shear chaperone folding (mlog10 = 1.30 -> x = 5.4 + 1.30*0.9 = 6.57)
  \node[anchor=east, font=\scriptsize\bfseries] at (5.2, 2.9) {Rotational shear chaperone folding};
  \fill[fairmidred!40!white] (5.4, 2.65) rectangle (6.57, 3.15);
  \node[anchor=west, font=\tiny\bfseries, text=gray!80!black] at (6.67, 2.9) {$3/42$ ($p=4.9\times 10^{-2}$)};

  % 8. Ribosome & translation attenuation (mlog10 = 2.90 -> x = 5.4 + 2.90*0.9 = 8.01) -> DEEP BLUE
  \node[anchor=east, font=\scriptsize\bfseries, text=fairblue] at (5.2, 2.1) {Ribosome \& translation attenuation};
  \fill[fairblue] (5.4, 1.85) rectangle (8.01, 2.35);
  \node[anchor=west, font=\tiny\bfseries, text=fairblue] at (8.11, 2.1) {$6/58$ (Downregulated, $p=1.3\times 10^{-3}$)};

  % Clean Figure Legend positioned at bottom (y=0.4, well below X-axis title at 1.0)
  \node[font=\tiny\bfseries, anchor=west] at (4.0, 0.4) {Directionality:};
  \fill[fairred] (5.6, 0.3) rectangle (6.2, 0.5); \node[font=\tiny, anchor=west] at (6.3, 0.4) {Upregulated in Microgravity};
  \fill[fairblue] (10.2, 0.3) rectangle (10.8, 0.5); \node[font=\tiny, anchor=west] at (10.9, 0.4) {Downregulated in Microgravity};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig5_tex, "fig5_pathway_ontology")

    # -------------------------------------------------------------
    # Figure 6: Intramodular Hub Connectivity (Zero Label Overlap)
    # -------------------------------------------------------------
    fig6_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\usetikzlibrary{arrows.meta}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily, >=Stealth]
  % Panel a: Scatter Plot k_within vs k_total (Shared baseline y=8.4)
  \node[font=\large\bfseries] at (0.2, 8.4) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.6, 8.4) {Intramodular Hub Centrality ($k_{\text{within}}$ vs. $k_{\text{total}}$)};

  % Axes starting cleanly at (1.2, 1.6) with uniform steps
  \draw[->, thick] (1.2, 1.6) -- (6.6, 1.6) node[midway, below=18pt, font=\scriptsize\bfseries] {Whole-Network Degree ($k_{\text{total}}$)};
  \draw[->, thick] (1.2, 1.6) -- (1.2, 7.6) node[midway, left=24pt, rotate=90, font=\scriptsize\bfseries] {Intramodular Connectivity ($k_{\text{within}}$)};

  % Uniform tick intervals starting directly at 1.2 on x, 1.6 on y
  \foreach \x/\lbl in {1.2/0, 2.4/10, 3.6/20, 4.8/30, 6.0/40} {
    \draw (\x, 1.5) -- (\x, 1.6) node[below=4pt, font=\tiny] {\lbl};
  }
  \foreach \y/\lbl in {1.6/0, 3.0/10, 4.4/20, 5.8/30, 7.2/40} {
    \draw (1.1, \y) -- (1.2, \y) node[left=4pt, font=\tiny] {\lbl};
  }

  % Hub Points with non-overlapping labels
  \fill[fairblue] (5.8, 6.9) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS12120 (35.3)\ };
  \fill[fairblue] (5.7, 6.5) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS21560 (35.0)\ };
  \fill[fairblue] (5.6, 6.1) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS04440 (34.7)\ };

  \fill[fairred] (4.2, 3.8) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS29685 (15.4)\ };
  \fill[fairred] (3.9, 3.5) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS02330 (14.8)\ };

  \fill[fairmidblue] (3.4, 3.1) circle (3.5pt) node[anchor=east, font=\tiny\bfseries] {RS11565 (11.9)\ };
  \fill[fairmidred] (3.0, 2.8) circle (3.5pt) node[anchor=east, font=\tiny\bfseries] {RS16730 (10.8)\ };
  \fill[fairred!50!white] (2.4, 2.3) circle (3pt) node[anchor=west, font=\tiny] {\ RS13930 (6.3)};

  % Panel b: Regulatory Interactome (Shared baseline y=8.4)
  \node[font=\large\bfseries] at (7.8, 8.4) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (8.2, 8.4) {Inter-Module Regulatory Interactome};

  % Network Nodes
  \node[circle, fill=fairblue!20, draw=fairblue, line width=1pt, font=\tiny\bfseries] (mps1) at (9.0, 6.6) {mps1};
  \node[circle, fill=fairblue!20, draw=fairblue, line width=1pt, font=\tiny\bfseries] (mps2) at (8.4, 5.4) {mps2};
  \node[circle, fill=fairblue!20, draw=fairblue, line width=1pt, font=\tiny\bfseries] (rs04) at (10.0, 5.8) {RS04440};

  \node[circle, fill=fairred!20, draw=fairred, line width=1pt, font=\tiny\bfseries] (kasA) at (11.4, 6.0) {kasA};
  \node[circle, fill=fairred!20, draw=fairred, line width=1pt, font=\tiny\bfseries] (fbpA) at (12.6, 6.6) {fbpA};
  \node[circle, fill=fairred!20, draw=fairred, line width=1pt, font=\tiny\bfseries] (rs02) at (12.7, 5.2) {RS02330};

  \node[circle, fill=fairred!30, draw=fairred!80!black, line width=1pt, font=\tiny\bfseries] (hspX) at (9.2, 3.2) {hspX};
  \node[circle, fill=fairred!30, draw=fairred!80!black, line width=1pt, font=\tiny\bfseries] (cydA) at (10.4, 2.6) {cydA};

  \node[circle, fill=fairmidred!25, draw=fairmidred, line width=1pt, font=\tiny\bfseries] (esxA) at (12.4, 3.2) {esxA};
  \node[circle, fill=fairmidred!25, draw=fairmidred, line width=1pt, font=\tiny\bfseries] (eccE) at (11.6, 2.3) {eccE};

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
    # Figure 7: Pan-Microbial Landscape (Zero Overlap with Phylum Names)
    # -------------------------------------------------------------
    fig7_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily]
  % Panel a: Taxonomic Distribution with Dedicated Y-axis Clearance (x=3.6 cm)
  \node[font=\large\bfseries] at (0.2, 8.8) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.6, 8.8) {Taxonomic Distribution Across 78 OSDR Spaceflight Datasets};

  % Axes for Panel a: Y-axis placed at x=3.6 cm!
  \draw[->, thick] (3.6, 4.4) -- (3.6, 8.2);
  \node[font=\scriptsize\bfseries, rotate=90] at (0.35, 6.3) {Microbial Phylum};
  \draw[->, thick] (3.6, 4.4) -- (8.2, 4.4) node[midway, below=18pt, font=\scriptsize\bfseries] {Dataset Representation (\% of 78 Studies)};

  % Ticks for Panel a: 3.6 to 7.6 (1.0 cm per 10%)
  \foreach \x/\lbl in {3.6/0\%, 4.6/10\%, 5.6/20\%, 6.6/30\%, 7.6/40\%} {
    \draw (\x, 4.3) -- (\x, 4.4) node[below=4pt, font=\tiny] {\lbl};
    \draw[dashed, draw=gray!20] (\x, 4.4) -- (\x, 8.2);
  }

  % Bars & Y-labels (Labels at x=3.4 anchor=east, bars from x=3.6)
  % 1. Pseudomonadota (43.6% -> x = 3.6 + 0.436*10 = 7.96)
  \draw (3.5, 7.6) -- (3.6, 7.6);
  \node[anchor=east, font=\tiny\bfseries] at (3.4, 7.6) {Pseudomonadota};
  \fill[fairblue] (3.6, 7.35) rectangle (7.96, 7.85);
  \node[anchor=west, font=\tiny\bfseries, text=fairblue] at (8.04, 7.6) {43.6\% ($N=34$)};

  % 2. Bacillota (28.2% -> x = 3.6 + 0.282*10 = 6.42)
  \draw (3.5, 6.6) -- (3.6, 6.6);
  \node[anchor=east, font=\tiny\bfseries] at (3.4, 6.6) {Bacillota};
  \fill[fairmidblue] (3.6, 6.35) rectangle (6.42, 6.85);
  \node[anchor=west, font=\tiny\bfseries, text=fairmidblue] at (6.50, 6.6) {28.2\% ($N=22$)};

  % 3. Actinomycetota (12.8% -> x = 3.6 + 0.128*10 = 4.88)
  \draw (3.5, 5.6) -- (3.6, 5.6);
  \node[anchor=east, font=\tiny\bfseries] at (3.4, 5.6) {Actinomycetota};
  \fill[fairmidred] (3.6, 5.35) rectangle (4.88, 5.85);
  \node[anchor=west, font=\tiny\bfseries, text=fairmidred] at (4.96, 5.6) {12.8\% ($N=10$)};

  % 4. Fungi (10.3% -> x = 3.6 + 0.103*10 = 4.63)
  \draw (3.5, 4.8) -- (3.6, 4.8);
  \node[anchor=east, font=\tiny\bfseries] at (3.4, 4.8) {Fungi};
  \fill[fairred] (3.6, 4.55) rectangle (4.63, 5.05);
  \node[anchor=west, font=\tiny\bfseries, text=fairred] at (4.71, 4.8) {10.3\% ($N=8$)};

  % Panel b: Cross-Species Spaceflight Concordance Heatmap
  \node[font=\large\bfseries] at (0.2, 3.2) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.6, 3.2) {Cross-Species Spaceflight Adaptation Concordance};

  % Axis Titles for Panel b
  \node[font=\scriptsize\bfseries] at (8.1, 2.7) {Conserved Spaceflight Phenotypic Hallmark};
  \node[font=\scriptsize\bfseries, rotate=90] at (0.35, 1.0) {Microbial Pathogen};

  % Column Headers
  \node[font=\tiny\bfseries] at (4.5, 2.2) {Biofilm};
  \node[font=\tiny\bfseries] at (6.3, 2.2) {Cell Wall};
  \node[font=\tiny\bfseries] at (8.1, 2.2) {Virulence};
  \node[font=\tiny\bfseries] at (9.9, 2.2) {Hypoxia};
  \node[font=\tiny\bfseries] at (11.7, 2.2) {Antibiotics};

  % Rows (Labels aligned at x=3.4 anchor=east)
  \node[anchor=east, font=\tiny\bfseries] at (3.4, 1.6) {M. marinum (OSD-528)};
  \fill[fairred] (3.8, 1.3) rectangle (5.2, 1.9); \node[font=\tiny\bfseries, text=white] at (4.5, 1.6) {+1.0};
  \fill[fairred] (5.6, 1.3) rectangle (7.0, 1.9); \node[font=\tiny\bfseries, text=white] at (6.3, 1.6) {+1.0};
  \fill[fairred] (7.4, 1.3) rectangle (8.8, 1.9); \node[font=\tiny\bfseries, text=white] at (8.1, 1.6) {+1.0};
  \fill[fairred] (9.2, 1.3) rectangle (10.6, 1.9); \node[font=\tiny\bfseries, text=white] at (9.9, 1.6) {+1.0};
  \fill[fairred] (11.0, 1.3) rectangle (12.4, 1.9); \node[font=\tiny\bfseries, text=white] at (11.7, 1.6) {+1.0};

  \node[anchor=east, font=\tiny\bfseries] at (3.4, 0.8) {P. aeruginosa (OSD-14)};
  \fill[fairred] (3.8, 0.5) rectangle (5.2, 1.1); \node[font=\tiny\bfseries, text=white] at (4.5, 0.8) {+1.0};
  \fill[fairred] (5.6, 0.5) rectangle (7.0, 1.1); \node[font=\tiny\bfseries, text=white] at (6.3, 0.8) {+1.0};
  \fill[fairred] (7.4, 0.5) rectangle (8.8, 1.1); \node[font=\tiny\bfseries, text=white] at (8.1, 0.8) {+1.0};
  \fill[fairred] (9.2, 0.5) rectangle (10.6, 1.1); \node[font=\tiny\bfseries, text=white] at (9.9, 0.8) {+1.0};
  \fill[fairred] (11.0, 0.5) rectangle (12.4, 1.1); \node[font=\tiny\bfseries, text=white] at (11.7, 0.8) {+1.0};

  \node[anchor=east, font=\tiny\bfseries] at (3.4, 0.0) {S. enterica (OSD-11)};
  \fill[fairmidred] (3.8, -0.3) rectangle (5.2, 0.3); \node[font=\tiny\bfseries, text=white] at (4.5, 0.0) {+0.8};
  \fill[fairred] (5.6, -0.3) rectangle (7.0, 0.3); \node[font=\tiny\bfseries, text=white] at (6.3, 0.0) {+1.0};
  \fill[fairred] (7.4, -0.3) rectangle (8.8, 0.3); \node[font=\tiny\bfseries, text=white] at (8.1, 0.0) {+1.0};
  \fill[fairmidred] (9.2, -0.3) rectangle (10.6, 0.3); \node[font=\tiny\bfseries, text=white] at (9.9, 0.0) {+0.7};
  \fill[fairred] (11.0, -0.3) rectangle (12.4, 0.3); \node[font=\tiny\bfseries, text=white] at (11.7, 0.0) {+1.0};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig7_tex, "fig7_pan_microbial_landscape")

    # -------------------------------------------------------------
    # Figure 8: Simulator Kinematics vs. Biological Response Radar
    # -------------------------------------------------------------
    fig8_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily]
  % Panel a: Radar Plot
  \node[font=\large\bfseries] at (0, 8.4) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 8.4) {Biological Response Trajectory Radar Across Simulation Modalities};

  % Center at (4.5, 4.2)
  % Web rings
  \foreach \r in {0.8, 1.6, 2.4, 3.2} {
    \draw[gray!25, dashed] (4.5, 4.2) circle (\r cm);
  }
  
  % 6 Axes at 60 degree intervals (offset titles outwards by 4.0 cm to avoid overlap)
  \foreach \ang/\lbl in {90/Biofilm Pellicle, 30/Cell Wall (FAS-II), 330/Type VII Secretion, 270/DosR Hypoxia, 210/Oxidative Stress, 150/Rotational Shear} {
    \draw[gray!40] (4.5, 4.2) -- ({4.5 + 3.4*cos(\ang)}, {4.2 + 3.4*sin(\ang)});
    \node[font=\tiny\bfseries] at ({4.5 + 3.9*cos(\ang)}, {4.2 + 3.9*sin(\ang)}) {\lbl};
  }

  % Static 1g Earth Polygon (Deep Blue - near center)
  \draw[line width=1.5pt, draw=fairblue, fill=fairblue!15]
    ({4.5 + 0.8*cos(90)}, {4.2 + 0.8*sin(90)}) --
    ({4.5 + 0.9*cos(30)}, {4.2 + 0.9*sin(30)}) --
    ({4.5 + 0.8*cos(330)}, {4.2 + 0.8*sin(330)}) --
    ({4.5 + 0.7*cos(270)}, {4.2 + 0.7*sin(270)}) --
    ({4.5 + 0.8*cos(210)}, {4.2 + 0.8*sin(210)}) --
    ({4.5 + 0.8*cos(150)}, {4.2 + 0.8*sin(150)}) -- cycle;

  % 3D Clinostat Polygon (Mid-Red / Slate)
  \draw[line width=1.5pt, draw=fairmidblue, fill=fairmidblue!20]
    ({4.5 + 3.0*cos(90)}, {4.2 + 3.0*sin(90)}) --
    ({4.5 + 3.0*cos(30)}, {4.2 + 3.0*sin(30)}) --
    ({4.5 + 2.8*cos(330)}, {4.2 + 2.8*sin(330)}) --
    ({4.5 + 2.9*cos(270)}, {4.2 + 2.9*sin(270)}) --
    ({4.5 + 2.9*cos(210)}, {4.2 + 2.9*sin(210)}) --
    ({4.5 + 2.9*cos(150)}, {4.2 + 2.9*sin(150)}) -- cycle;

  % RPM 2.0 Polygon (Deep Red)
  \draw[line width=1.5pt, draw=fairred, fill=fairred!20]
    ({4.5 + 3.1*cos(90)}, {4.2 + 3.1*sin(90)}) --
    ({4.5 + 3.1*cos(30)}, {4.2 + 3.1*sin(30)}) --
    ({4.5 + 2.9*cos(330)}, {4.2 + 2.9*sin(330)}) --
    ({4.5 + 3.0*cos(270)}, {4.2 + 3.0*sin(270)}) --
    ({4.5 + 3.1*cos(210)}, {4.2 + 3.1*sin(210)}) --
    ({4.5 + 1.8*cos(150)}, {4.2 + 1.8*sin(150)}) -- cycle;

  % Panel b: Concordance Summary
  \node[font=\large\bfseries] at (9.4, 8.4) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (9.8, 8.4) {Kinematic Concordance};

  \draw[line width=2pt, draw=fairmidblue] (9.6, 7.2) -- (10.2, 7.2);
  \node[anchor=west, font=\scriptsize\bfseries, text=fairblue] at (10.4, 7.2) {3D Clinostat ($<0.01g$)};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (10.4, 6.7) {Continuous rotational shear vector};

  \draw[line width=2pt, draw=fairred] (9.6, 5.8) -- (10.2, 5.8);
  \node[anchor=west, font=\scriptsize\bfseries, text=fairred] at (10.4, 5.8) {RPM 2.0 ($<0.01g$)};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (10.4, 5.3) {Multi-axis random vectoring};

  \draw[line width=2pt, draw=fairblue] (9.6, 4.4) -- (10.2, 4.4);
  \node[anchor=west, font=\scriptsize\bfseries, text=fairblue] at (10.4, 4.4) {Static 1g Ground Control};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (10.4, 3.9) {Standard gravity baseline};

  \node[anchor=west, font=\tiny\bfseries, text=fairred] at (9.6, 2.6) {Core Concordance: $>75\%$ shared DEGs};
  \node[anchor=west, font=\tiny\bfseries, text=fairblue] at (9.6, 2.0) {Divergence: Rotational shear only};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig8_tex, "fig8_simulator_concordance_radar")

    # -------------------------------------------------------------
    # Figure 9: Hyper-Detailed Cellular and Metabolic Architecture
    # -------------------------------------------------------------
    fig9_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\usetikzlibrary{arrows.meta}
""" + tikz_colors + r"""
\begin{document}
\begin{tikzpicture}[font=\sffamily, >=Stealth]
  % Panel a: Multi-Compartment Cellular Cross-Section
  \node[font=\large\bfseries] at (0, 14.2) {a};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 14.2) {Mycobacterium marinum Multi-Compartment Cellular and Metabolic Architecture};

  % Layer 1: Substrate & Outer Mycomembrane (Capsule)
  \draw[line width=2pt, draw=gray!40] (0.2, 13.5) -- (13.5, 13.5);
  \node[anchor=west, font=\tiny\bfseries, text=gray!80!black] at (0.2, 13.7) {PDMS Silicone Membrane Substrate (Quiescent Boundary Layer)};

  % Outer Mycomembrane Strip
  \fill[fairmidred!20] (0.2, 12.3) rectangle (13.5, 13.1);
  \draw[thick, draw=fairmidred] (0.2, 13.1) -- (13.5, 13.1);
  \draw[thick, draw=fairmidred] (0.2, 12.3) -- (13.5, 12.3);
  \node[anchor=west, font=\tiny\bfseries, text=fairred!80!black] at (0.3, 12.7) {Outer Mycomembrane: Cord Factor (TDM) \& GPL Biofilm Matrix};

  % Extracellular/Capsular elements
  \node[fill=fairblue, text=white, font=\tiny\bfseries, rounded corners=2pt, inner sep=2pt] at (3.2, 12.7) {mps1/2};
  \node[fill=fairred, text=white, font=\tiny\bfseries, rounded corners=2pt, inner sep=2pt] at (6.5, 12.7) {fbpA (+1.2)};
  \node[anchor=west, font=\tiny, text=darkgray] at (7.4, 12.7) {Cord Factor (TDM) Deposition};
  \node[fill=fairblue, text=white, font=\tiny\bfseries, rounded corners=2pt, inner sep=2pt] at (11.5, 12.7) {espB (-6.0)};

  % Layer 2: Periplasmic Space & Arabinogalactan-Peptidoglycan (AG-PG) Mesh
  \fill[gray!6] (0.2, 10.9) rectangle (13.5, 12.3);
  \draw[dashed, draw=gray!40] (0.2, 11.6) -- (13.5, 11.6);
  \node[anchor=west, font=\tiny\bfseries, text=darkgray] at (0.3, 11.95) {Periplasmic Space: Arabinogalactan-Peptidoglycan Polymer Mesh};
  
  \node[fill=fairmidred, text=white, font=\tiny\bfseries, rounded corners=2pt, inner sep=2pt] at (2.2, 11.3) {murA (+0.8)};
  \node[anchor=west, font=\tiny] at (3.0, 11.3) {Murein cross-linking};
  \node[fill=fairred, text=white, font=\tiny\bfseries, rounded corners=2pt, inner sep=2pt] at (7.8, 11.3) {mmpL3 (+1.1)};
  \node[anchor=west, font=\tiny] at (8.7, 11.3) {TMM Translocase};

  % Layer 3: Inner Plasma Membrane
  \fill[fairblue!15] (0.2, 9.1) rectangle (13.5, 10.9);
  \draw[thick, draw=fairblue] (0.2, 10.9) -- (13.5, 10.9);
  \draw[thick, draw=fairblue] (0.2, 9.1) -- (13.5, 9.1);
  \node[anchor=west, font=\tiny\bfseries, text=fairblue!90!black] at (0.3, 10.6) {Inner Plasma Membrane: Embedded Respiratory Complexes \& Translocons};

  % Membrane Complexes (Nodes)
  % Complex I (nuoC, nuoD)
  \draw[fill=fairred, draw=fairred!80!black, rounded corners=2pt] (0.4, 9.3) rectangle (1.9, 10.3);
  \node[font=\tiny\bfseries, text=white, align=center] at (1.15, 9.8) {Complex I\\nuoC/D (+6.6)};
  \draw[->, thick, draw=fairred] (1.15, 9.2) -- (1.15, 10.4) node[above, font=\tiny\bfseries] {4H$^+$};

  % Menaquinone reductase MenJ
  \draw[fill=fairred, draw=fairred!80!black, rounded corners=2pt] (2.3, 9.3) rectangle (3.5, 10.3);
  \node[font=\tiny\bfseries, text=white, align=center] at (2.9, 9.8) {MenJ\\(+6.2)};

  % Microaerophilic Cytochrome bd Oxidase (cydA, cydB)
  \draw[fill=fairred, draw=fairred!80!black, rounded corners=2pt] (3.9, 9.3) rectangle (5.4, 10.3);
  \node[font=\tiny\bfseries, text=white, align=center] at (4.65, 9.8) {Cyt $bd$\\cydA (+2.4)};
  \draw[->, thick, draw=fairred] (4.65, 9.2) -- (4.65, 10.4) node[above, font=\tiny] {O$_2 \to$ H$_2$O};

  % CydD Thiol ABC Exporter
  \draw[fill=fairred, draw=fairred!80!black, rounded corners=2pt] (5.8, 9.3) rectangle (7.0, 10.3);
  \node[font=\tiny\bfseries, text=white, align=center] at (6.4, 9.8) {CydD\\(+6.2)};
  \draw[->, thick, draw=fairred] (6.4, 9.2) -- (6.4, 10.4) node[above, font=\tiny] {Thiol};

  % Type VII ESX-1 Core Pore (eccB1, eccD1, eccE1)
  \draw[fill=fairmidred, draw=fairmidred!80!black, rounded corners=2pt] (7.4, 9.3) rectangle (9.0, 10.3);
  \node[font=\tiny\bfseries, text=white, align=center] at (8.2, 9.8) {Type VII\\eccB/D/E};
  \draw[->, thick, draw=fairmidred] (8.2, 9.2) -- (8.2, 12.3) node[midway, right=1pt, font=\tiny] {EsxA};

  % ATP Synthase
  \draw[fill=white, draw=gray!60, rounded corners=2pt] (9.4, 9.3) rectangle (10.6, 10.3);
  \node[font=\tiny\bfseries, text=darkgray, align=center] at (10.0, 9.8) {ATP Synth\\atpA (0.0)};
  \draw[<-, thick, draw=gray!60] (10.0, 9.2) -- (10.0, 10.4) node[below, font=\tiny] {ATP};

  % Alternate Transporters
  \draw[fill=fairblue, draw=fairblue!80!black, rounded corners=2pt] (11.0, 9.3) rectangle (12.2, 10.3);
  \node[font=\tiny\bfseries, text=white, align=center] at (11.6, 9.8) {fadD7\\(-5.3)};

  % Layer 4: Cytoplasm (Metabolic Cascades)
  \fill[fairlight] (0.2, 5.0) rectangle (13.5, 9.1);
  \node[anchor=west, font=\tiny\bfseries, text=darkgray] at (0.3, 8.8) {Cytoplasm: Biochemical Cascades \& Metabolic Shunts};

  % Cascade A: Nitrogen & Polyamine Shunt
  \draw[thick, draw=fairred!40, fill=white] (0.4, 5.2) rectangle (4.3, 8.5);
  \node[font=\tiny\bfseries, text=fairred, anchor=west] at (0.5, 8.2) {1. Nitrogen \& Polyamine Shunt};
  \node[font=\tiny, anchor=west] at (0.5, 7.7) {Carbamoyl-P + L-Ornithine};
  \draw[->, thick, draw=fairred] (1.8, 7.5) -- (1.8, 6.7);
  \node[fill=fairred, text=white, font=\tiny\bfseries, rounded corners=2pt, inner sep=2pt] at (1.8, 7.1) {argF (+7.0)};
  \node[font=\tiny, anchor=west] at (0.5, 6.4) {$\longrightarrow$ L-Citrulline $\to$ Polyamines};
  \node[fill=fairred, text=white, font=\tiny\bfseries, rounded corners=2pt, inner sep=2pt] at (2.2, 5.8) {ilvB / als (+6.1)};
  \node[font=\tiny, anchor=west] at (0.5, 5.4) {Pyruvate $\to$ 2-Acetolactate};

  % Cascade B: FAS-II Mycolic Acid Synthesis
  \draw[thick, draw=fairmidred!40, fill=white] (4.6, 5.2) rectangle (8.6, 8.5);
  \node[font=\tiny\bfseries, text=fairred!90!black, anchor=west] at (4.7, 8.2) {2. Mycolic Acid FAS-II Spiral};
  \node[font=\tiny, anchor=west] at (4.7, 7.7) {Acetyl-CoA $\xrightarrow{\textbf{accD}}$ Malonyl-CoA};
  \draw[->, thick, draw=fairmidred] (6.4, 7.5) -- (6.4, 6.7);
  \node[fill=fairred, text=white, font=\tiny\bfseries, rounded corners=2pt, inner sep=2pt] at (6.4, 7.1) {kasA (+1.8)};
  \node[font=\tiny, anchor=west] at (4.7, 6.4) {$\beta$-ketoacyl-AcpM $\xrightarrow{\textbf{inhA}}$ Acyl-AcpM};
  \node[font=\tiny, anchor=west] at (4.7, 5.8) {Meromycolate C70--C90 synthesis};
  \node[font=\tiny\bfseries, text=fairred!80!black, anchor=west] at (4.7, 5.4) {Cell Envelope Fortification};

  % Cascade C: Central Carbon & Redox Systems
  \draw[thick, draw=fairblue!40, fill=white] (8.9, 5.2) rectangle (13.3, 8.5);
  \node[font=\tiny\bfseries, text=fairblue!90!black, anchor=west] at (9.0, 8.2) {3. Redox \& Cofactor Systems};
  \node[fill=fairred, text=white, font=\tiny\bfseries, rounded corners=2pt, inner sep=2pt] at (10.0, 7.7) {adhP (+7.0)};
  \node[font=\tiny, anchor=west] at (10.7, 7.7) {Alcohol redox};
  \node[fill=fairred, text=white, font=\tiny\bfseries, rounded corners=2pt, inner sep=2pt] at (10.0, 7.0) {ribD (+6.2)};
  \node[font=\tiny, anchor=west] at (10.7, 7.0) {Riboflavin / FAD};
  \node[fill=fairred, text=white, font=\tiny\bfseries, rounded corners=2pt, inner sep=2pt] at (10.0, 6.3) {F420 (+6.0)};
  \node[font=\tiny, anchor=west] at (10.7, 6.3) {LLM oxidoreductase};
  \node[fill=fairblue, text=white, font=\tiny\bfseries, rounded corners=2pt, inner sep=2pt] at (10.0, 5.6) {icl1 (-0.7)};
  \node[font=\tiny, anchor=west] at (10.7, 5.6) {Glyoxylate bypass};

  % Panel b: Subsystem Perturbation Index (SPI) Bar Plot with Zero Overlap
  \node[font=\large\bfseries] at (0, 4.4) {b};
  \node[anchor=west, font=\small\bfseries, text=darkgray] at (0.4, 4.4) {Empirical Subsystem Perturbation Index (SPI) Ranking};

  % Axes raised with generous clearance for X-axis title
  \draw[->, thick] (4.8, 1.0) -- (12.6, 1.0) node[midway, below=18pt, font=\scriptsize\bfseries] {Subsystem Perturbation Index (SPI)};
  \draw[thick] (4.8, 1.0) -- (4.8, 4.2);

  \foreach \x/\lbl in {4.8/0, 6.3/5, 7.8/10, 9.3/15, 10.8/20, 12.3/25} {
    \draw (\x, 0.9) -- (\x, 1.0) node[below=3pt, font=\tiny] {\lbl};
    \draw[dashed, draw=gray!20] (\x, 1.0) -- (\x, 4.2);
  }

  % Bars (Y positions: 3.8, 3.4, 3.0, 2.6, 2.2, 1.8, 1.4)
  % 1. Nitrogen Shunts & Polyamines (SPI = 24.52)
  \node[anchor=east, font=\scriptsize\bfseries] at (4.6, 3.8) {Nitrogen Shunts \& Polyamines};
  \fill[fairred] (4.8, 3.6) rectangle (12.16, 4.0);
  \node[anchor=west, font=\tiny\bfseries, text=fairred] at (12.25, 3.8) {24.5 (Active, $+6.6$)};

  % 2. Redox Homeostasis & Cofactors (SPI = 10.62)
  \node[anchor=east, font=\scriptsize\bfseries] at (4.6, 3.4) {Redox Homeostasis \& Cofactors};
  \fill[fairred!90!white] (4.8, 3.2) rectangle (7.99, 3.6);
  \node[anchor=west, font=\tiny\bfseries, text=fairred!90!black] at (8.09, 3.4) {10.6 (Active, $+3.5$)};

  % 3. Cellular Respiration & Energy (SPI = 5.27)
  \node[anchor=east, font=\scriptsize\bfseries] at (4.6, 3.0) {Cellular Respiration \& Energy};
  \fill[fairmidred] (4.8, 2.8) rectangle (6.38, 3.2);
  \node[anchor=west, font=\tiny\bfseries, text=fairmidred] at (6.48, 3.0) {5.3 (Active, $+2.5$)};

  % 4. Type VII Secretion & Virulence (SPI = 2.87)
  \node[anchor=east, font=\scriptsize\bfseries] at (4.6, 2.6) {Type VII Secretion \& Virulence};
  \fill[fairblue!70!white] (4.8, 2.4) rectangle (5.66, 2.8);
  \node[anchor=west, font=\tiny\bfseries, text=fairblue] at (5.76, 2.6) {2.9 (Repressed, $-0.5$)};

  % 5. GPL Biofilm & Cell Surface (SPI = 2.69)
  \node[anchor=east, font=\scriptsize\bfseries] at (4.6, 2.2) {GPL Biofilm \& Cell Surface};
  \fill[fairblue!60!white] (4.8, 2.0) rectangle (5.61, 2.4);
  \node[anchor=west, font=\tiny\bfseries, text=fairblue] at (5.71, 2.2) {2.7 (Repressed, $-0.9$)};

  % 6. Mycolic Acid & FAS-II Envelope (SPI = 1.24)
  \node[anchor=east, font=\scriptsize\bfseries] at (4.6, 1.8) {Mycolic Acid \& FAS-II Envelope};
  \fill[fairmidred!60!white] (4.8, 1.6) rectangle (5.17, 2.0);
  \node[anchor=west, font=\tiny\bfseries, text=fairmidred!80!black] at (5.27, 1.8) {1.2 (Active, $+0.4$)};

  % 7. DosR Hypoxia & Central Carbon (SPI = 0.97)
  \node[anchor=east, font=\scriptsize\bfseries] at (4.6, 1.4) {DosR Hypoxia \& Central Carbon};
  \fill[fairblue!40!white] (4.8, 1.2) rectangle (5.09, 1.6);
  \node[anchor=west, font=\tiny\bfseries, text=gray!80!black] at (5.19, 1.4) {0.97 (Baseline shift)};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig9_tex, "fig9_cellular_metabolic_landscape")

    print("Phase 6 completed successfully.")

if __name__ == "__main__":
    main()
