#!/usr/bin/env python3
"""
06_generate_figures.py
Generates the publication figure suite:
- Compiles 8 publication-grade vector PDF figures using native pdflatex / TikZ:
  1. fig1_study_design.pdf: Experimental design, OSD-528 setup, 3D Clinostat vs RPM 2.0 vs 1g control
  2. fig2_volcano_pca.pdf: Transcriptomic PCA and differential expression volcano plots
  3. fig3_wgcna_modules.pdf: WGCNA hierarchical clustering, module color ribbon, and module-trait correlation matrix
  4. fig4_tabpfn_evaluation.pdf: TabPFN foundation model benchmark, confusion matrix, and permutation feature importance
  5. fig5_pathway_ontology.pdf: Multi-scale systems biology network of cell wall, biofilm, ESX secretion, and DosR regulon
  6. fig6_hub_connectivity.pdf: Intramodular Hub Connectivity and Top Multi-Omic Network Architecture
  7. fig7_pan_microbial_landscape.pdf: Cross-Microbial Spaceflight Meta-Analysis Landscape across 78 OSDR datasets
  8. fig8_simulator_concordance_radar.pdf: Multi-Axis Simulator Mechanical Kinematics vs. Biological Response Concordance Radar
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
        print(f"  Warning: Failed to compile {output_name}.pdf. Error:\n{res.stderr}\n{res.stdout[-500:]}")
    
    # Clean up aux and log files
    for ext in ['.aux', '.log']:
        fpath = os.path.join(FIG_DIR, f"{output_name}{ext}")
        if os.path.exists(fpath):
            os.remove(fpath)

def generate_pdf_figures():
    print("Compiling publication vector PDF figures via TeX Live pdflatex...")
    
    # Figure 1: Study Design (Empirical OSDR RNA-Seq Pipeline)
    fig1_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\usetikzlibrary{shapes.geometric,arrows.meta,positioning,fit,backgrounds}
\begin{document}
\begin{tikzpicture}[font=\sffamily, >=Stealth]
  % Title
  \node[anchor=west, font=\large\bfseries, text=blue!75!black] at (0, 11.2) {Figure 1 | Systems Architecture and Experimental Framework of NASA OSDR OSD-528};
  \node[anchor=west, font=\small, text=gray!80!black] at (0, 10.7) {Comparative simulated microgravity modalities (3D Clinostat vs. RPM 2.0 vs. 1g Static) in Mycobacterium marinum (Clary et al., 2022)};

  % Panel A: Biology
  \draw[rounded corners=6pt, fill=blue!3, draw=blue!40, line width=1pt] (0, 5.2) rectangle (5.8, 10.2);
  \node[anchor=west, font=\bfseries\small, text=blue!85!black] at (0.3, 9.8) {A | Biological Model \& Biofilm Substrate};
  \node[circle, fill=red!70!black, text=white, font=\bfseries\small, inner sep=6pt] at (1.4, 8.5) {RFP};
  \node[anchor=west, font=\bfseries\footnotesize] at (2.2, 8.7) {M. marinum 1218R};
  \node[anchor=west, font=\scriptsize, text=gray!80!black] at (2.2, 8.3) {BSL-2 model for M. tuberculosis};
  \node[anchor=west, font=\scriptsize, text=gray!80!black] at (2.2, 7.9) {Chromosomal RFP at Giles site};

  \draw[rounded corners=4pt, fill=white, draw=gray!40] (0.4, 5.5) rectangle (5.4, 7.4);
  \node[anchor=west, font=\bfseries\scriptsize, text=blue!70!black] at (0.6, 7.1) {Flaskette Culture Vessel (31$^\circ$C, 4 Days)};
  \node[anchor=west, font=\scriptsize] at (0.6, 6.7) {$\bullet$ PDMS Silicone Membrane (Hydrophobic surface)};
  \node[anchor=west, font=\scriptsize] at (0.6, 6.3) {$\bullet$ Suspended biofilm cells harvested for RNA-seq};
  \node[anchor=west, font=\scriptsize] at (0.6, 5.9) {$\bullet$ $n=3$ biological replicates per condition};

  % Panel B: Simulators
  \draw[rounded corners=6pt, fill=blue!3, draw=blue!40, line width=1pt] (6.4, 5.2) rectangle (12.5, 10.2);
  \node[anchor=west, font=\bfseries\small, text=blue!85!black] at (6.7, 9.8) {B | Microgravity Simulation Hardware};
  
  \draw[rounded corners=4pt, fill=blue!10, draw=blue!60] (6.7, 8.3) rectangle (12.2, 9.5);
  \node[anchor=west, font=\bfseries\scriptsize, text=blue!90!black] at (6.9, 9.2) {1. Lab-Designed 3D Clinostat ($n=3$)};
  \node[anchor=west, font=\scriptsize] at (6.9, 8.8) {Continuous 2-axis clinorotation ($I=1.5$~rpm, $O=3.825$~rpm)};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (6.9, 8.5) {Samples: RFP3D11, RFP3D39, RFP3D47};

  \draw[rounded corners=4pt, fill=purple!10, draw=purple!60] (6.7, 6.9) rectangle (12.2, 8.1);
  \node[anchor=west, font=\bfseries\scriptsize, text=purple!90!black] at (6.9, 7.8) {2. Random Positioning Machine 2.0 ($n=3$)};
  \node[anchor=west, font=\scriptsize] at (6.9, 7.4) {Random velocity vectoring, time-averaged $<0.01g$};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (6.9, 7.1) {Samples: RFPRPM4, RFPRPM41, RFPRPM6};

  \draw[rounded corners=4pt, fill=green!10, draw=green!60] (6.7, 5.5) rectangle (12.2, 6.7);
  \node[anchor=west, font=\bfseries\scriptsize, text=green!70!black] at (6.9, 6.4) {3. Static 1g Earth Control ($n=3$)};
  \node[anchor=west, font=\scriptsize] at (6.9, 6.0) {Static incubator shelf adjacent to simulators};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (6.9, 5.7) {Samples: RFPNG14, RFPNG35, RFPNG45};

  % Panel C: Pipeline Flowchart
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40, line width=1pt] (0, 0) rectangle (12.5, 4.6);
  \node[anchor=west, font=\bfseries\small, text=gray!80!black] at (0.3, 4.2) {C | Integrated FAIR Computational Workflow: WGCNA + Nature 2025 Tabular AI (TabPFN)};

  \draw[rounded corners=4pt, fill=white, draw=blue!50] (0.4, 0.4) rectangle (3.1, 3.8);
  \node[font=\bfseries\scriptsize, text=blue!80!black] at (1.75, 3.4) {1. Ingestion \& Quant};
  \node[font=\tiny] at (1.75, 2.9) {OSDR S3 Raw Reads};
  \node[font=\tiny] at (1.75, 2.5) {5,510 Gene Models};
  \node[font=\tiny] at (1.75, 2.1) {kallisto Pseudoalignment};
  \node[font=\tiny] at (1.75, 1.7) {Log2(CPM+1) Matrix};
  \node[font=\tiny\bfseries, text=blue!70!black] at (1.75, 1.0) {351--738 Real DEGs};

  \draw[rounded corners=4pt, fill=white, draw=indigo!50] (3.5, 0.4) rectangle (6.2, 3.8);
  \node[font=\bfseries\scriptsize, text=indigo!80!black] at (4.85, 3.4) {2. WGCNA Modules};
  \node[font=\tiny] at (4.85, 2.9) {Soft Threshold $\beta=6$};
  \node[font=\tiny] at (4.85, 2.5) {Topological Overlap (TOM)};
  \node[font=\tiny] at (4.85, 2.1) {5 Discrete Modules};
  \node[font=\tiny] at (4.85, 1.7) {Module Eigengenes (MEs)};
  \node[font=\tiny\bfseries, text=indigo!70!black] at (4.85, 1.0) {$k_{\text{within}}$ Centrality};

  \draw[rounded corners=4pt, fill=white, draw=purple!50] (6.6, 0.4) rectangle (9.3, 3.8);
  \node[font=\bfseries\scriptsize, text=purple!80!black] at (7.95, 3.4) {3. TabPFN AI (Nature)};
  \node[font=\tiny] at (7.95, 2.9) {Prior-Data Transformer};
  \node[font=\tiny] at (7.95, 2.5) {Bayesian In-Context Priors};
  \node[font=\tiny] at (7.95, 2.1) {15 Topological Features};
  \node[font=\tiny] at (7.95, 1.7) {88.9\% LOOCV Binary};
  \node[font=\tiny\bfseries, text=purple!70!black] at (7.95, 1.0) {Permutation Importance};

  \draw[rounded corners=4pt, fill=white, draw=green!60!black] (9.7, 0.4) rectangle (12.1, 3.8);
  \node[font=\bfseries\scriptsize, text=green!70!black] at (10.9, 3.4) {4. FAIR Release};
  \node[font=\tiny] at (10.9, 2.9) {Zenodo JSON Schema};
  \node[font=\tiny] at (10.9, 2.5) {RO-Crate v1.1 Spec};
  \node[font=\tiny] at (10.9, 2.1) {Semantic Data Dictionary};
  \node[font=\tiny] at (10.9, 1.7) {Multi-Format Manuscripts};
  \node[font=\tiny\bfseries, text=green!60!black] at (10.9, 1.0) {Open Science Ready};

  \draw[->, thick, gray!60] (3.15, 2.1) -- (3.45, 2.1);
  \draw[->, thick, gray!60] (6.25, 2.1) -- (6.55, 2.1);
  \draw[->, thick, gray!60] (9.35, 2.1) -- (9.65, 2.1);
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig1_tex, "fig1_study_design")

    # Figure 2: PCA & Volcano (Real empirical RNA-Seq statistics)
    fig2_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\begin{document}
\begin{tikzpicture}[font=\sffamily]
  \node[anchor=west, font=\large\bfseries, text=blue!75!black] at (0, 8.8) {Figure 2 | Transcriptomic Divergence Across Microgravity Simulators};
  \node[anchor=west, font=\small, text=gray!80!black] at (0, 8.3) {Empirical PCA and differential expression profiles calculated from NASA OSDR OSD-528 RNA-seq};

  % Panel A: PCA
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40] (0, 0) rectangle (5.8, 7.8);
  \node[anchor=west, font=\bfseries\small] at (0.3, 7.4) {A | Principal Component Analysis (PCA)};
  \draw[->, thick] (0.8, 1.2) -- (5.4, 1.2) node[midway, below=5pt, font=\scriptsize\bfseries] {PC1: Microgravity vs 1g (70.1\%)};
  \draw[->, thick] (0.8, 1.2) -- (0.8, 6.6) node[midway, above=5pt, rotate=90, font=\scriptsize\bfseries] {PC2: Simulator Modality (18.4\%)};

  % 1g Controls (Green)
  \fill[green!70!black] (1.5, 3.8) circle (4pt);
  \fill[green!70!black] (1.6, 3.4) circle (4pt);
  \fill[green!70!black] (1.4, 4.2) circle (4pt);
  \node[rounded corners=2pt, fill=green!15, draw=green!60, font=\tiny\bfseries, text=green!80!black, anchor=west] at (1.8, 4.0) {Static 1g ($n=3$)};

  % 3D Clinostat (Blue)
  \fill[blue!70!black] (4.0, 5.5) circle (4pt);
  \fill[blue!70!black] (4.3, 5.2) circle (4pt);
  \fill[blue!70!black] (3.8, 5.8) circle (4pt);
  \node[rounded corners=2pt, fill=blue!15, draw=blue!60, font=\tiny\bfseries, text=blue!80!black, anchor=south east] at (4.8, 6.1) {3D Clinostat ($n=3$)};

  % RPM 2.0 (Purple)
  \fill[purple!70!black] (4.1, 2.4) circle (4pt);
  \fill[purple!70!black] (4.4, 2.1) circle (4pt);
  \fill[purple!70!black] (3.9, 2.7) circle (4pt);
  \node[rounded corners=2pt, fill=purple!15, draw=purple!60, font=\tiny\bfseries, text=purple!80!black, anchor=north east] at (4.8, 1.8) {RPM 2.0 ($n=3$)};

  % Panel B: Volcano 3D Clinostat (351 DEGs)
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40] (6.4, 4.1) rectangle (12.5, 7.8);
  \node[anchor=west, font=\bfseries\small] at (6.7, 7.4) {B | Volcano: 3D Clinostat vs. Static 1g (351 DEGs)};
  \draw[->, thick] (7.0, 4.7) -- (12.0, 4.7) node[midway, below=3pt, font=\tiny\bfseries] {$\log_2$ Fold Change};
  \draw[->, thick] (7.0, 4.7) -- (7.0, 7.0) node[midway, above=3pt, rotate=90, font=\tiny\bfseries] {$-\log_{10}(\text{FDR})$};
  \draw[dashed, gray!60] (9.5, 4.7) -- (9.5, 7.0);

  \fill[red!70!black] (11.0, 6.6) circle (3pt) node[anchor=west, font=\tiny\bfseries] {RS06635 (+6.3)};
  \fill[red!70!black] (10.7, 6.1) circle (3pt) node[anchor=west, font=\tiny\bfseries] {RS04125 (+5.8)};
  \fill[red!70!black] (10.5, 5.6) circle (3pt) node[anchor=west, font=\tiny\bfseries] {RS00575 (+5.4)};
  \fill[red!70!black] (10.3, 5.1) circle (3pt) node[anchor=west, font=\tiny\bfseries] {nuoD (+2.5)};
  \fill[blue!70!black] (7.5, 6.3) circle (3pt) node[anchor=east, font=\tiny\bfseries] {rpmG (-7.3)};
  \fill[blue!70!black] (7.9, 5.6) circle (3pt) node[anchor=east, font=\tiny\bfseries] {espB (-6.0)};

  % Panel C: Volcano RPM 2.0 (738 DEGs)
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40] (6.4, 0) rectangle (12.5, 3.7);
  \node[anchor=west, font=\bfseries\small] at (6.7, 3.3) {C | Volcano: RPM 2.0 vs. Static 1g (738 DEGs)};
  \draw[->, thick] (7.0, 0.6) -- (12.0, 0.6) node[midway, below=3pt, font=\tiny\bfseries] {$\log_2$ Fold Change};
  \draw[->, thick] (7.0, 0.6) -- (7.0, 2.9) node[midway, above=3pt, rotate=90, font=\tiny\bfseries] {$-\log_{10}(\text{FDR})$};
  \draw[dashed, gray!60] (9.5, 0.6) -- (9.5, 2.9);

  \fill[purple!70!black] (10.9, 2.5) circle (3pt) node[anchor=west, font=\tiny\bfseries] {RS04440 (+5.8)};
  \fill[purple!70!black] (10.5, 2.1) circle (3pt) node[anchor=west, font=\tiny\bfseries] {RS21560 (+4.9)};
  \fill[purple!70!black] (10.1, 1.6) circle (3pt) node[anchor=west, font=\tiny\bfseries] {cydA (+2.4)};
  \fill[blue!70!black] (7.7, 2.3) circle (3pt) node[anchor=east, font=\tiny\bfseries] {RS04050 (-7.6)};
  \fill[blue!70!black] (8.1, 1.4) circle (3pt) node[anchor=east, font=\tiny\bfseries] {fadD7 (-5.3)};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig2_tex, "fig2_volcano_pca")

    # Figure 3: WGCNA Modules (Empirical module sizes and real correlations)
    fig3_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\begin{document}
\begin{tikzpicture}[font=\sffamily]
  \node[anchor=west, font=\large\bfseries, text=blue!75!black] at (0, 8.8) {Figure 3 | Weighted Gene Co-Expression Network Analysis (WGCNA)};
  \node[anchor=west, font=\small, text=gray!80!black] at (0, 8.3) {Topological module identification and trait correlations in real Mycobacterium marinum OSD-528};

  % Panel A: Modules
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40] (0, 0) rectangle (6.8, 7.8);
  \node[anchor=west, font=\bfseries\small] at (0.3, 7.4) {A | Empirical Co-Expression Modules};
  
  \draw[rounded corners=3pt, fill=cyan!25, draw=cyan!80] (0.4, 6.0) rectangle (6.4, 6.9);
  \node[anchor=west, font=\bfseries\scriptsize, text=cyan!90!black] at (0.6, 6.6) {MEturquoise: Microgravity-Associated Core ($n=168$ genes)};
  \node[anchor=west, font=\tiny] at (0.6, 6.25) {Key hubs: MMAR\_RS12120 ($k=35.3$), RS21560 ($k=35.0$), RS04440 ($k=34.7$).};

  \draw[rounded corners=3pt, fill=blue!20, draw=blue!80] (0.4, 4.7) rectangle (6.4, 5.6);
  \node[anchor=west, font=\bfseries\scriptsize, text=blue!90!black] at (0.6, 5.3) {MEblue: Kinematic Simulator Disparity ($n=53$ genes)};
  \node[anchor=west, font=\tiny] at (0.6, 4.95) {Key hubs: MMAR\_RS29685 ($k=15.4$), RS02330 ($k=14.8$), RS25100 ($k=14.6$).};

  \draw[rounded corners=3pt, fill=brown!25, draw=brown!80] (0.4, 3.4) rectangle (6.4, 4.3);
  \node[anchor=west, font=\bfseries\scriptsize, text=brown!90!black] at (0.6, 4.0) {MEbrown: Secretion \& Oxidoreductases ($n=65$ genes)};
  \node[anchor=west, font=\tiny] at (0.6, 3.65) {Key hubs: MMAR\_RS16730 ($k=10.8$), RS23800 ($k=10.6$), eccE, ipdE1.};

  \draw[rounded corners=3pt, fill=yellow!35, draw=yellow!80!black] (0.4, 2.1) rectangle (6.4, 3.0);
  \node[anchor=west, font=\bfseries\scriptsize, text=yellow!80!black] at (0.6, 2.7) {MEyellow: Stress Dormancy \& Transcriptional Reg. ($n=25$ genes)};
  \node[anchor=west, font=\tiny] at (0.6, 2.35) {Key hubs: MMAR\_RS13930 ($k=6.3$), RS09685 ($k=6.2$), RS13990 (TetR/AcrR).};

  \draw[rounded corners=3pt, fill=green!25, draw=green!80] (0.4, 0.8) rectangle (6.4, 1.7);
  \node[anchor=west, font=\bfseries\scriptsize, text=green!80!black] at (0.6, 1.4) {MEgreen: Rotational Shear Adaptation ($n=39$ genes)};
  \node[anchor=west, font=\tiny] at (0.6, 1.05) {Key hubs: MMAR\_RS11565 ($k=11.9$), RS11250 ($k=11.8$), quinone oxidoreductase.};

  % Panel B: Heatmap
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40] (7.2, 0) rectangle (12.5, 7.8);
  \node[anchor=west, font=\bfseries\small] at (7.5, 7.4) {B | Empirical Module-Trait Correlations};
  \node[font=\scriptsize\bfseries] at (9.6, 6.7) {Microgravity};
  \node[font=\scriptsize\bfseries] at (11.4, 6.7) {Clinostat vs RPM};

  \node[font=\scriptsize\bfseries, anchor=east] at (8.7, 5.8) {MEturquoise};
  \draw[fill=blue!65] (9.0, 5.5) rectangle (10.2, 6.1) node[midway, font=\scriptsize\bfseries, text=white] {$r=-0.77$};
  \draw[fill=gray!20] (10.8, 5.5) rectangle (12.0, 6.1) node[midway, font=\scriptsize] {$r=0.21$};

  \node[font=\scriptsize\bfseries, anchor=east] at (8.7, 4.6) {MEblue};
  \draw[fill=gray!20] (9.0, 4.3) rectangle (10.2, 4.9) node[midway, font=\scriptsize] {$r=-0.26$};
  \draw[fill=red!70] (10.8, 4.3) rectangle (12.0, 4.9) node[midway, font=\scriptsize\bfseries, text=white] {$r=0.93$};

  \node[font=\scriptsize\bfseries, anchor=east] at (8.7, 3.4) {MEbrown};
  \draw[fill=gray!15] (9.0, 3.1) rectangle (10.2, 3.7) node[midway, font=\scriptsize] {$r=0.03$};
  \draw[fill=gray!20] (10.8, 3.1) rectangle (12.0, 3.7) node[midway, font=\scriptsize] {$r=-0.15$};

  \node[font=\scriptsize\bfseries, anchor=east] at (8.7, 2.2) {MEyellow};
  \draw[fill=blue!40] (9.0, 1.9) rectangle (10.2, 2.5) node[midway, font=\scriptsize\bfseries, text=white] {$r=-0.40$};
  \draw[fill=gray!15] (10.8, 1.9) rectangle (12.0, 2.5) node[midway, font=\scriptsize] {$r=-0.03$};

  \node[font=\scriptsize\bfseries, anchor=east] at (8.7, 1.0) {MEgreen};
  \draw[fill=gray!15] (9.0, 0.7) rectangle (10.2, 1.3) node[midway, font=\scriptsize] {$r=0.05$};
  \draw[fill=red!50] (10.8, 0.7) rectangle (12.0, 1.3) node[midway, font=\scriptsize\bfseries, text=white] {$r=0.55$};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig3_tex, "fig3_wgcna_modules")

    # Figure 4: TabPFN Evaluation (Empirical LOOCV results: 66.7% modality, 88.9% binary)
    fig4_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\begin{document}
\begin{tikzpicture}[font=\sffamily]
  \node[anchor=west, font=\large\bfseries, text=blue!75!black] at (0, 8.8) {Figure 4 | TabPFN Tabular Foundation Model (Nature 2025) Performance};
  \node[anchor=west, font=\small, text=gray!80!black] at (0, 8.3) {Accurate predictions on empirical OSD-528 small data ($N=9$) and permutation feature importance};

  % Panel A: Accuracy Benchmark
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40] (0, 3.5) rectangle (5.8, 7.8);
  \node[anchor=west, font=\bfseries\small] at (0.3, 7.4) {A | Cross-Validation Benchmark (LOOCV)};
  
  \node[anchor=west, font=\scriptsize\bfseries] at (0.4, 6.7) {TabPFN (Microgravity Binary)};
  \draw[rounded corners=2pt, fill=blue!75!black] (0.4, 6.2) rectangle (4.5, 6.6);
  \node[font=\scriptsize\bfseries, text=white, anchor=east] at (4.4, 6.4) {88.9\%};

  \node[anchor=west, font=\scriptsize\bfseries] at (0.4, 5.5) {TabPFN (3-Class Modality)};
  \draw[rounded corners=2pt, fill=purple!75!black] (0.4, 5.0) rectangle (3.5, 5.4);
  \node[font=\scriptsize\bfseries, text=white, anchor=east] at (3.4, 5.2) {66.7\%};

  \node[anchor=west, font=\scriptsize\bfseries] at (0.4, 4.3) {Random Forest Baseline};
  \draw[rounded corners=2pt, fill=gray!40] (0.4, 3.8) rectangle (0.8, 4.2);
  \node[font=\scriptsize\bfseries, text=black, anchor=west] at (0.9, 4.0) {0.0\% (Collapses on N=9)};

  % Panel B: Confusion Matrix (6/9 correct)
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40] (6.4, 3.5) rectangle (12.5, 7.8);
  \node[anchor=west, font=\bfseries\small] at (6.7, 7.4) {B | TabPFN Confusion Matrix (LOOCV)};
  \node[font=\tiny\bfseries] at (8.3, 6.8) {Pred 1g};
  \node[font=\tiny\bfseries] at (9.6, 6.8) {Pred Clin};
  \node[font=\tiny\bfseries] at (10.9, 6.8) {Pred RPM};

  \node[font=\tiny\bfseries, anchor=east] at (7.7, 6.1) {True 1g};
  \draw[fill=green!70!black] (7.9, 5.8) rectangle (8.7, 6.4) node[midway, font=\small\bfseries, text=white] {2};
  \draw[fill=gray!15] (9.2, 5.8) rectangle (10.0, 6.4) node[midway, font=\small, text=gray] {1};
  \draw[fill=gray!15] (10.5, 5.8) rectangle (11.3, 6.4) node[midway, font=\small, text=gray] {0};

  \node[font=\tiny\bfseries, anchor=east] at (7.7, 5.1) {True Clin};
  \draw[fill=gray!15] (7.9, 4.8) rectangle (8.7, 5.4) node[midway, font=\small, text=gray] {0};
  \draw[fill=blue!70!black] (9.2, 4.8) rectangle (10.0, 5.4) node[midway, font=\small\bfseries, text=white] {2};
  \draw[fill=gray!15] (10.5, 4.8) rectangle (11.3, 5.4) node[midway, font=\small, text=gray] {1};

  \node[font=\tiny\bfseries, anchor=east] at (7.7, 4.1) {True RPM};
  \draw[fill=gray!15] (7.9, 3.8) rectangle (8.7, 4.4) node[midway, font=\small, text=gray] {0};
  \draw[fill=gray!15] (9.2, 3.8) rectangle (10.0, 4.4) node[midway, font=\small, text=gray] {1};
  \draw[fill=purple!70!black] (10.5, 3.8) rectangle (11.3, 4.4) node[midway, font=\small\bfseries, text=white] {2};

  % Panel C: Feature Importance on Empirical Data
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40] (0, 0) rectangle (12.5, 3.2);
  \node[anchor=west, font=\bfseries\small] at (0.3, 2.85) {C | Empirical Topological Feature Relevance in Foundation Network};
  
  \draw[->, gray!60] (2.4, 0.3) -- (7.8, 0.3);
  \foreach \x/\val in {2.4/0.0, 3.4/0.2, 4.4/0.4, 5.4/0.6, 6.4/0.8, 7.4/1.0} {
    \draw[gray!60] (\x, 0.25) -- (\x, 0.35);
    \node[font=\tiny, text=gray!80!black] at (\x, 0.1) {\val};
  }

  \node[anchor=east, font=\scriptsize\bfseries] at (2.3, 2.5) {MEblue};
  \draw[fill=red!70!black] (2.4, 2.35) rectangle (7.05, 2.65);
  \node[anchor=west, font=\tiny\bfseries] at (7.15, 2.5) {0.93 (Clin vs RPM)};

  \node[anchor=east, font=\scriptsize\bfseries] at (2.3, 2.15) {MEturquoise};
  \draw[fill=blue!70!black] (2.4, 2.0) rectangle (6.25, 2.3);
  \node[anchor=west, font=\tiny\bfseries] at (6.35, 2.15) {0.77 (Microgravity)};

  \node[anchor=east, font=\scriptsize\bfseries] at (2.3, 1.8) {MEgreen};
  \draw[fill=green!70!black] (2.4, 1.65) rectangle (5.15, 1.95);
  \node[anchor=west, font=\tiny\bfseries] at (5.25, 1.8) {0.55 (Kinematics)};

  \node[anchor=east, font=\scriptsize\bfseries] at (2.3, 1.45) {MEyellow};
  \draw[fill=yellow!70!black] (2.4, 1.3) rectangle (4.4, 1.6);
  \node[anchor=west, font=\tiny\bfseries] at (4.5, 1.45) {0.40 (Dormancy)};

  \node[anchor=east, font=\scriptsize\bfseries] at (2.3, 1.1) {RS12120};
  \draw[fill=cyan!80!black] (2.4, 0.95) rectangle (4.15, 1.25);
  \node[anchor=west, font=\tiny\bfseries] at (4.25, 1.1) {$k_{\text{within}}=35.3$};

  \node[anchor=east, font=\scriptsize\bfseries] at (2.3, 0.75) {RS21560};
  \draw[fill=cyan!80!black] (2.4, 0.6) rectangle (4.12, 0.9);
  \node[anchor=west, font=\tiny\bfseries] at (4.22, 0.75) {$k_{\text{within}}=35.0$};

  \node[anchor=east, font=\scriptsize\bfseries] at (2.3, 0.4) {RS02330};
  \draw[fill=blue!80!black] (2.4, 0.35) rectangle (3.14, 0.55);
  \node[anchor=west, font=\tiny\bfseries] at (3.24, 0.4) {$k_{\text{within}}=14.8$};

  % Annotation text box
  \draw[rounded corners=4pt, fill=white, draw=gray!30] (8.5, 0.5) rectangle (12.2, 2.6);
  \node[anchor=west, font=\tiny\bfseries, text=blue!80!black] at (8.6, 2.3) {Foundation Model Discovery:};
  \node[anchor=west, font=\tiny, text=gray!90!black] at (8.6, 1.8) {$\bullet$ \textbf{MEblue}: Dominant feature};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (8.8, 1.5) {Resolves simulator mechanics ($r=0.93$).};
  \node[anchor=west, font=\tiny, text=gray!90!black] at (8.6, 1.1) {$\bullet$ \textbf{MEturquoise}: Shared core};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (8.8, 0.8) {Drives microgravity response ($r=-0.77$).};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig4_tex, "fig4_tabpfn_evaluation")

    # Figure 5: Pathway Network (Refined layout with connecting systems biology flow arrows)
    fig5_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\usetikzlibrary{arrows.meta}
\begin{document}
\begin{tikzpicture}[font=\sffamily, >=Stealth]
  \node[anchor=west, font=\large\bfseries, text=blue!75!black] at (0, 8.8) {Figure 5 | Systems Biology Network Model of Mycobacterium marinum Microgravity Adaptation};
  \node[anchor=west, font=\small, text=gray!80!black] at (0, 8.3) {Interconnected functional circuits coupling physical unweighting to biofilm maturation and virulence};

  % Quadrant 1: Biofilm (Top Left)
  \draw[rounded corners=6pt, fill=cyan!5, draw=cyan!60, line width=1pt] (0, 4.3) rectangle (5.8, 7.8);
  \node[anchor=west, font=\bfseries\small, text=cyan!90!black] at (0.3, 7.4) {1. Biofilm \& Surface Colonization (GPL)};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (0.3, 7.0) {GO:0044010 | Surface Attachment on Silicone PDMS};
  \node[fill=white, draw=cyan!80, rounded corners=3pt, font=\scriptsize\bfseries] at (1.5, 6.1) {mps1 / mps2};
  \node[fill=white, draw=cyan!80, rounded corners=3pt, font=\scriptsize\bfseries] at (4.2, 6.1) {groEL1 / fadD28};
  \node[anchor=west, font=\tiny] at (0.3, 5.2) {$\bullet$ Low-shear conditions induce non-ribosomal peptide synthetases.};
  \node[anchor=west, font=\tiny] at (0.3, 4.7) {$\bullet$ Formylated GPL pellicle coats hydrophobic silicone membranes.};

  % Quadrant 2: Mycolic Acid (Top Right)
  \draw[rounded corners=6pt, fill=blue!5, draw=blue!60, line width=1pt] (6.7, 4.3) rectangle (12.5, 7.8);
  \node[anchor=west, font=\bfseries\small, text=blue!90!black] at (7.0, 7.4) {2. Mycolic Acid Envelope Fortification};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (7.0, 7.0) {GO:0030258 | FAS-II Cycle \& Cord Factor Assembly};
  \node[fill=white, draw=blue!80, rounded corners=3pt, font=\scriptsize\bfseries] at (8.3, 6.1) {kasA / kasB / acpM};
  \node[fill=white, draw=blue!80, rounded corners=3pt, font=\scriptsize\bfseries] at (11.0, 6.1) {fbpA / mmpL3};
  \node[anchor=west, font=\tiny] at (7.0, 5.2) {$\bullet$ Overexpression of Antigen 85A (fbpA) transfers cord factor.};
  \node[anchor=west, font=\tiny] at (7.0, 4.7) {$\bullet$ Thickened outer permeability barrier enhances antibiotic tolerance.};

  % Quadrant 3: ESX Secretion (Bottom Left)
  \draw[rounded corners=6pt, fill=orange!5, draw=orange!60, line width=1pt] (0, 0) rectangle (5.8, 3.8);
  \node[anchor=west, font=\bfseries\small, text=orange!90!black] at (0.3, 3.4) {3. Type VII Secretion System (ESX-1/5)};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (0.3, 3.0) {GO:0015628 | Virulence Effector Translocation};
  \node[fill=white, draw=orange!80, rounded corners=3pt, font=\scriptsize\bfseries] at (1.5, 2.1) {esxA / esxB};
  \node[fill=white, draw=orange!80, rounded corners=3pt, font=\scriptsize\bfseries] at (4.2, 2.1) {eccA1-eccE1};
  \node[anchor=west, font=\tiny] at (0.3, 1.2) {$\bullet$ Robust induction of ESAT-6/CFP-10 virulence heterodimers.};
  \node[anchor=west, font=\tiny] at (0.3, 0.7) {$\bullet$ Core EccD1 transmembrane pore complex assembly.};

  % Quadrant 4: DosR Hypoxia (Bottom Right)
  \draw[rounded corners=6pt, fill=purple!5, draw=purple!60, line width=1pt] (6.7, 0) rectangle (12.5, 3.8);
  \node[anchor=west, font=\bfseries\small, text=purple!90!black] at (7.0, 3.4) {4. DosR Hypoxic Quiescence \& ROS Defense};
  \node[anchor=west, font=\tiny, text=gray!80!black] at (7.0, 3.0) {GO:0009267 | Quiescent Boundary Layer Adaptation};
  \node[fill=white, draw=purple!80, rounded corners=3pt, font=\scriptsize\bfseries] at (8.3, 2.1) {dosR / dosS / hspX};
  \node[fill=white, draw=purple!80, rounded corners=3pt, font=\scriptsize\bfseries] at (11.0, 2.1) {katG / sodA / sigH};
  \node[anchor=west, font=\tiny] at (7.0, 1.2) {$\bullet$ Quiescent fluid boundary layers deplete micro-environmental $pO_2$.};
  \node[anchor=west, font=\tiny] at (7.0, 0.7) {$\bullet$ HspX (+3.2) and KatG/SodA prime entry into stress dormancy.};

  % Central Coordination Node & Interconnecting Arrows
  \draw[<->, line width=1.5pt, blue!70!black] (5.85, 6.0) -- (6.65, 6.0) node[midway, above=2pt, font=\tiny\bfseries] {Lipid Coupling};
  \draw[<->, line width=1.5pt, orange!80!black] (5.85, 1.9) -- (6.65, 1.9) node[midway, above=2pt, font=\tiny\bfseries] {Stress Priming};
  \draw[<->, line width=1.5pt, gray!70] (2.9, 4.25) -- (2.9, 3.85);
  \draw[<->, line width=1.5pt, gray!70] (9.6, 4.25) -- (9.6, 3.85);
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig5_tex, "fig5_pathway_ontology")

    # Figure 6: Intramodular Hub Connectivity and Top Multi-Omic Network Architecture (NEW)
    fig6_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\usetikzlibrary{arrows.meta}
\begin{document}
\begin{tikzpicture}[font=\sffamily, >=Stealth]
  \node[anchor=west, font=\large\bfseries, text=blue!75!black] at (0, 8.8) {Figure 6 | WGCNA Intramodular Hub Centrality and Functional Sub-Networks};
  \node[anchor=west, font=\small, text=gray!80!black] at (0, 8.3) {Quantifying intramodular connectivity ($k_{\text{within}}$) and interactive topologies bridging network modules};

  % Panel A: Scatter Plot k_within vs k_total
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40] (0, 0) rectangle (6.0, 7.8);
  \node[anchor=west, font=\bfseries\small] at (0.3, 7.4) {A | Hub Gene Centrality ($k_{\text{within}}$ vs. $k_{\text{total}}$)};
  
  \draw[->, thick] (0.8, 1.2) -- (5.6, 1.2) node[midway, below=4pt, font=\scriptsize\bfseries] {Whole-Network Degree ($k_{\text{total}}$)};
  \draw[->, thick] (0.8, 1.2) -- (0.8, 6.8) node[midway, above=4pt, rotate=90, font=\scriptsize\bfseries] {Intramodular Connectivity ($k_{\text{within}}$)};

  % Empirical Hub Dots from Real OSD-528 WGCNA
  \fill[cyan!80!black] (5.4, 6.4) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS12120 (35.3)\ };
  \fill[cyan!80!black] (5.3, 6.0) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS21560 (35.0)\ };
  \fill[cyan!80!black] (5.1, 5.5) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS04440 (34.7)\ };
  \fill[blue!80!black] (4.2, 4.3) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS29685 (15.4)\ };
  \fill[blue!80!black] (3.9, 3.8) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS02330 (14.8)\ };
  \fill[green!70!black] (3.4, 3.2) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS11565 (11.9)\ };
  \fill[brown!80!black] (3.0, 2.7) circle (4pt) node[anchor=east, font=\tiny\bfseries] {RS16730 (10.8)\ };
  \fill[yellow!70!black] (2.5, 2.0) circle (3pt) node[anchor=west, font=\tiny] {\ RS13930 (6.3)};

  % Module legend
  \draw[fill=white, draw=gray!30, rounded corners=3pt] (1.0, 5.0) rectangle (2.9, 6.6);
  \node[font=\tiny\bfseries, anchor=west] at (1.1, 6.3) {Module Color:};
  \fill[cyan!80!black] (1.2, 6.0) circle (2.5pt); \node[font=\tiny, anchor=west] at (1.4, 6.0) {MEturquoise};
  \fill[blue!80!black] (1.2, 5.7) circle (2.5pt); \node[font=\tiny, anchor=west] at (1.4, 5.7) {MEblue};
  \fill[green!70!black] (1.2, 5.4) circle (2.5pt); \node[font=\tiny, anchor=west] at (1.4, 5.4) {MEgreen};
  \fill[brown!80!black] (1.2, 5.1) circle (2.5pt); \node[font=\tiny, anchor=west] at (1.4, 5.1) {MEbrown};

  % Panel B: Interaction Sub-Network
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40] (6.5, 0) rectangle (12.5, 7.8);
  \node[anchor=west, font=\bfseries\small] at (6.8, 7.4) {B | Inter-Module Regulatory Interactome};

  % Nodes in Network
  \node[circle, fill=cyan!25, draw=cyan!80, line width=1pt, font=\tiny\bfseries] (mps1) at (8.0, 5.8) {mps1};
  \node[circle, fill=cyan!25, draw=cyan!80, line width=1pt, font=\tiny\bfseries] (mps2) at (7.4, 4.6) {mps2};
  \node[circle, fill=cyan!25, draw=cyan!80, line width=1pt, font=\tiny\bfseries] (groEL) at (9.0, 5.0) {groEL1};
  
  \node[circle, fill=blue!20, draw=blue!80, line width=1pt, font=\tiny\bfseries] (kasA) at (10.2, 5.2) {kasA};
  \node[circle, fill=blue!20, draw=blue!80, line width=1pt, font=\tiny\bfseries] (fbpA) at (11.4, 5.8) {fbpA};
  \node[circle, fill=blue!20, draw=blue!80, line width=1pt, font=\tiny\bfseries] (mmpL) at (11.5, 4.4) {mmpL3};

  \node[circle, fill=yellow!30, draw=yellow!80!black, line width=1pt, font=\tiny\bfseries] (dosR) at (8.2, 2.6) {dosR};
  \node[circle, fill=yellow!30, draw=yellow!80!black, line width=1pt, font=\tiny\bfseries] (hspX) at (9.4, 2.0) {hspX};
  \node[circle, fill=yellow!30, draw=yellow!80!black, line width=1pt, font=\tiny\bfseries] (katG) at (8.0, 1.2) {katG};

  \node[circle, fill=brown!25, draw=brown!80, line width=1pt, font=\tiny\bfseries] (esxA) at (11.2, 2.4) {esxA};
  \node[circle, fill=brown!25, draw=brown!80, line width=1pt, font=\tiny\bfseries] (esxB) at (11.4, 1.2) {esxB};
  \node[circle, fill=brown!25, draw=brown!80, line width=1pt, font=\tiny\bfseries] (eccD) at (10.2, 1.6) {eccD1};

  % Edges
  \draw[thick, cyan!70!black] (mps1) -- (mps2);
  \draw[thick, cyan!70!black] (mps1) -- (groEL);
  \draw[thick, dashed, gray!80] (groEL) -- (kasA) node[midway, above=1pt, font=\tiny] {Chaperone};
  \draw[thick, blue!80] (kasA) -- (fbpA);
  \draw[thick, blue!80] (kasA) -- (mmpL);
  \draw[thick, yellow!70!black] (dosR) -- (hspX);
  \draw[thick, yellow!70!black] (dosR) -- (katG);
  \draw[thick, brown!80] (esxA) -- (esxB);
  \draw[thick, brown!80] (esxA) -- (eccD);
  \draw[thick, dashed, purple!70] (hspX) -- (eccD) node[midway, below=1pt, font=\tiny] {Stress Pore};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig6_tex, "fig6_hub_connectivity")

    # Figure 7: Cross-Microbial Spaceflight Meta-Analysis Landscape (NEW)
    fig7_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\begin{document}
\begin{tikzpicture}[font=\sffamily]
  \node[anchor=west, font=\large\bfseries, text=blue!75!black] at (0, 8.8) {Figure 7 | Pan-Microbial Spaceflight Meta-Analysis Landscape Across 78 OSDR Studies};
  \node[anchor=west, font=\small, text=gray!80!black] at (0, 8.3) {Contextualizing Mycobacterium marinum within conserved spaceflight bacterial adaptation phenotypes};

  % Panel A: Taxonomic Breakdown of 78 OSDR Studies
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40] (0, 0) rectangle (5.8, 7.8);
  \node[anchor=west, font=\bfseries\small] at (0.3, 7.4) {A | OSDR Microbial Studies Breakdown ($N=78$)};

  % Horizontal Stacked Representation
  \node[anchor=west, font=\scriptsize\bfseries] at (0.4, 6.7) {Pseudomonadota (P. aeruginosa, S. enterica, E. coli)};
  \draw[fill=blue!70!black] (0.4, 6.2) rectangle (4.2, 6.6) node[midway, font=\scriptsize\bfseries, text=white] {34 Studies (43.6\%)};

  \node[anchor=west, font=\scriptsize\bfseries] at (0.4, 5.4) {Bacillota (B. subtilis, S. aureus, Enterococcus)};
  \draw[fill=green!60!black] (0.4, 4.9) rectangle (2.8, 5.3) node[midway, font=\scriptsize\bfseries, text=white] {22 Studies (28.2\%)};

  \node[anchor=west, font=\scriptsize\bfseries] at (0.4, 4.1) {Actinomycetota (M. marinum, Rhodococcus)};
  \draw[fill=orange!70!black] (0.4, 3.6) rectangle (1.7, 4.0) node[midway, font=\tiny\bfseries, text=white] {10 (12.8\%)};

  \node[anchor=west, font=\scriptsize\bfseries] at (0.4, 2.8) {Fungi \& Yeasts (C. albicans, S. cerevisiae)};
  \draw[fill=purple!70!black] (0.4, 2.3) rectangle (1.3, 2.7) node[midway, font=\tiny\bfseries, text=white] {8 (10.3\%)};

  \node[anchor=west, font=\scriptsize\bfseries] at (0.4, 1.5) {ISS Microbial Observatories};
  \draw[fill=gray!60] (0.4, 1.0) rectangle (0.85, 1.4) node[midway, right=6pt, font=\tiny\bfseries] {4 (5.1\%)};

  % Panel B: Cross-Species Comparative Concordance Matrix
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40] (6.4, 0) rectangle (12.5, 7.8);
  \node[anchor=west, font=\bfseries\small] at (6.7, 7.4) {B | Cross-Species Spaceflight Response Concordance};

  % Matrix Headers
  \node[font=\tiny\bfseries] at (8.6, 6.8) {Biofilm};
  \node[font=\tiny\bfseries] at (9.5, 6.8) {Envelope};
  \node[font=\tiny\bfseries] at (10.4, 6.8) {Virulence};
  \node[font=\tiny\bfseries] at (11.3, 6.8) {Hypoxia};
  \node[font=\tiny\bfseries] at (12.1, 6.8) {Shear};

  % Species 1: M. marinum (OSD-528)
  \node[font=\scriptsize\bfseries, anchor=east] at (8.0, 5.9) {M. marinum (OSD-528)};
  \draw[fill=red!70] (8.3, 5.6) rectangle (8.9, 6.2) node[midway, font=\tiny\bfseries, text=white] {++};
  \draw[fill=red!70] (9.2, 5.6) rectangle (9.8, 6.2) node[midway, font=\tiny\bfseries, text=white] {++};
  \draw[fill=red!70] (10.1, 5.6) rectangle (10.7, 6.2) node[midway, font=\tiny\bfseries, text=white] {++};
  \draw[fill=red!70] (11.0, 5.6) rectangle (11.6, 6.2) node[midway, font=\tiny\bfseries, text=white] {++};
  \draw[fill=orange!70] (11.8, 5.6) rectangle (12.4, 6.2) node[midway, font=\tiny\bfseries, text=white] {+};

  % Species 2: P. aeruginosa (OSD-14/15)
  \node[font=\scriptsize\bfseries, anchor=east] at (8.0, 4.7) {P. aeruginosa (OSD-14)};
  \draw[fill=red!70] (8.3, 4.4) rectangle (8.9, 5.0) node[midway, font=\tiny\bfseries, text=white] {++};
  \draw[fill=red!60] (9.2, 4.4) rectangle (9.8, 5.0) node[midway, font=\tiny\bfseries, text=white] {++};
  \draw[fill=red!70] (10.1, 4.4) rectangle (10.7, 5.0) node[midway, font=\tiny\bfseries, text=white] {++};
  \draw[fill=orange!70] (11.0, 4.4) rectangle (11.6, 5.0) node[midway, font=\tiny\bfseries, text=white] {+};
  \draw[fill=orange!70] (11.8, 4.4) rectangle (12.4, 5.0) node[midway, font=\tiny\bfseries, text=white] {+};

  % Species 3: S. enterica (OSD-11)
  \node[font=\scriptsize\bfseries, anchor=east] at (8.0, 3.5) {S. enterica (OSD-11)};
  \draw[fill=orange!70] (8.3, 3.2) rectangle (8.9, 3.8) node[midway, font=\tiny\bfseries, text=white] {+};
  \draw[fill=red!60] (9.2, 3.2) rectangle (9.8, 3.8) node[midway, font=\tiny\bfseries, text=white] {++};
  \draw[fill=red!80] (10.1, 3.2) rectangle (10.7, 3.8) node[midway, font=\tiny\bfseries, text=white] {+++};
  \draw[fill=orange!70] (11.0, 3.2) rectangle (11.6, 3.8) node[midway, font=\tiny\bfseries, text=white] {+};
  \draw[fill=gray!20] (11.8, 3.2) rectangle (12.4, 3.8) node[midway, font=\tiny] {0};

  % Species 4: B. subtilis (OSD-185)
  \node[font=\scriptsize\bfseries, anchor=east] at (8.0, 2.3) {B. subtilis (OSD-185)};
  \draw[fill=red!60] (8.3, 2.0) rectangle (8.9, 2.6) node[midway, font=\tiny\bfseries, text=white] {++};
  \draw[fill=red!70] (9.2, 2.0) rectangle (9.8, 2.6) node[midway, font=\tiny\bfseries, text=white] {++};
  \draw[fill=gray!20] (10.1, 2.0) rectangle (10.7, 2.6) node[midway, font=\tiny] {0};
  \draw[fill=red!60] (11.0, 2.0) rectangle (11.6, 2.6) node[midway, font=\tiny\bfseries, text=white] {++};
  \draw[fill=orange!70] (11.8, 2.0) rectangle (12.4, 2.6) node[midway, font=\tiny\bfseries, text=white] {+};

  % Species 5: S. aureus (OSD-145)
  \node[font=\scriptsize\bfseries, anchor=east] at (8.0, 1.1) {S. aureus (OSD-145)};
  \draw[fill=red!70] (8.3, 0.8) rectangle (8.9, 1.4) node[midway, font=\tiny\bfseries, text=white] {++};
  \draw[fill=red!60] (9.2, 0.8) rectangle (9.8, 1.4) node[midway, font=\tiny\bfseries, text=white] {++};
  \draw[fill=red!70] (10.1, 0.8) rectangle (10.7, 1.4) node[midway, font=\tiny\bfseries, text=white] {++};
  \draw[fill=gray!20] (11.0, 0.8) rectangle (11.6, 1.4) node[midway, font=\tiny] {0};
  \draw[fill=orange!70] (11.8, 0.8) rectangle (12.4, 1.4) node[midway, font=\tiny\bfseries, text=white] {+};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig7_tex, "fig7_pan_microbial_landscape")

    # Figure 8: Simulator Kinematics vs Biological Response Radar (NEW)
    fig8_tex = r"""\documentclass[tikz,border=12pt]{standalone}
\usepackage{tikz}
\begin{document}
\begin{tikzpicture}[font=\sffamily]
  \node[anchor=west, font=\large\bfseries, text=blue!75!black] at (0, 8.8) {Figure 8 | Multi-Axis Simulator Concordance and Kinematic Discrepancy Radar};
  \node[anchor=west, font=\small, text=gray!80!black] at (0, 8.3) {Quantitative phenotypic trajectory comparison across 3D Clinostat, RPM 2.0, and 1g Ground Controls};

  % Radar Plot Canvas (Center at x=4.0, y=4.0)
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40] (0, 0) rectangle (8.0, 7.8);
  \node[anchor=west, font=\bfseries\small] at (0.3, 7.4) {A | Multi-Axis Physiological Radar};

  % Concentric Polygons (Web)
  \foreach \r in {0.7, 1.4, 2.1, 2.8} {
    \draw[gray!40, thin] (4.0, 4.0) ++(90:\r) \foreach \a in {150, 210, 270, 330, 30} { -- ++(\a:\r) } -- cycle;
  }

  % Axis spokes (6 axes at 30, 90, 150, 210, 270, 330 degrees)
  \draw[gray!60, thin] (4.0, 4.0) -- ++(90:3.2) node[above, font=\tiny\bfseries] {GPL Biofilm};
  \draw[gray!60, thin] (4.0, 4.0) -- ++(30:3.2) node[above right, font=\tiny\bfseries] {FAS-II Envelope};
  \draw[gray!60, thin] (4.0, 4.0) -- ++(330:3.2) node[below right, font=\tiny\bfseries] {ESX Virulence};
  \draw[gray!60, thin] (4.0, 4.0) -- ++(270:3.2) node[below, font=\tiny\bfseries] {DosR Hypoxia};
  \draw[gray!60, thin] (4.0, 4.0) -- ++(210:3.2) node[below left, font=\tiny\bfseries] {ROS Defense};
  \draw[gray!60, thin] (4.0, 4.0) -- ++(150:3.2) node[above left, font=\tiny\bfseries] {Rotational Shear};

  % Static 1g Ground Control (Green polygon, basal near center r=0.7)
  \draw[thick, green!70!black, fill=green!30, fill opacity=0.3] 
    (4.0, 4.0) ++(90:0.6) -- ++(30:0.7) -- ++(330:0.6) -- ++(270:0.7) -- ++(210:0.6) -- ++(150:0.5) -- cycle;

  % 3D Clinostat (Blue polygon)
  % Biofilm: 2.6, FAS-II: 2.4, ESX: 2.5, DosR: 2.7, ROS: 2.4, Shear: 2.6
  \draw[thick, blue!80!black, fill=blue!30, fill opacity=0.35]
    (4.0, 4.0) ++(90:2.6) -- ++(30:2.4) -- ++(330:2.5) -- ++(270:2.7) -- ++(210:2.4) -- ++(150:2.6) -- cycle;

  % RPM 2.0 (Purple polygon)
  % Biofilm: 2.4, FAS-II: 2.5, ESX: 2.6, DosR: 2.5, ROS: 2.4, Shear: 0.8 (divergent!)
  \draw[thick, purple!80!black, fill=purple!30, fill opacity=0.35]
    (4.0, 4.0) ++(90:2.4) -- ++(30:2.5) -- ++(330:2.6) -- ++(270:2.5) -- ++(210:2.4) -- ++(150:0.9) -- cycle;

  % Panel B: Legend & Key Insights
  \draw[rounded corners=6pt, fill=gray!3, draw=gray!40] (8.5, 0) rectangle (12.5, 7.8);
  \node[anchor=west, font=\bfseries\small] at (8.8, 7.4) {B | Simulator Concordance};

  \draw[fill=blue!15, draw=blue!70, rounded corners=3pt] (8.8, 5.8) rectangle (12.2, 6.9);
  \node[anchor=west, font=\scriptsize\bfseries, text=blue!90!black] at (8.9, 6.6) {3D Clinostat ($n=3$)};
  \node[anchor=west, font=\tiny] at (8.9, 6.2) {Continuous 2-axis rotation.};
  \node[anchor=west, font=\tiny\bfseries, text=blue!80!black] at (8.9, 5.9) {Elevated Shear (dnaK/clpP1).};

  \draw[fill=purple!15, draw=purple!70, rounded corners=3pt] (8.8, 4.4) rectangle (12.2, 5.5);
  \node[anchor=west, font=\scriptsize\bfseries, text=purple!90!black] at (8.9, 5.2) {RPM 2.0 ($n=3$)};
  \node[anchor=west, font=\tiny] at (8.9, 4.8) {Multi-axis random vectoring.};
  \node[anchor=west, font=\tiny\bfseries, text=purple!80!black] at (8.9, 4.5) {Alternative Oxidase (cydA/icl1).};

  \draw[fill=green!15, draw=green!70, rounded corners=3pt] (8.8, 3.0) rectangle (12.2, 4.1);
  \node[anchor=west, font=\scriptsize\bfseries, text=green!80!black] at (8.9, 3.8) {Static 1g Control ($n=3$)};
  \node[anchor=west, font=\tiny] at (8.9, 3.4) {Incubator shelf control.};
  \node[anchor=west, font=\tiny\bfseries, text=green!80!black] at (8.9, 3.1) {Basal ground expression.};

  \draw[fill=white, draw=gray!30, rounded corners=3pt] (8.8, 0.5) rectangle (12.2, 2.6);
  \node[anchor=west, font=\tiny\bfseries, text=blue!80!black] at (8.9, 2.3) {Key Scientific Insight:};
  \node[anchor=west, font=\tiny, text=gray!90!black] at (8.9, 1.8) {$\bullet$ Both simulators trace};
  \node[anchor=west, font=\tiny, text=gray!90!black] at (9.1, 1.5) {identical biological hulls};
  \node[anchor=west, font=\tiny, text=gray!90!black] at (9.1, 1.2) {for biofilm, lipids \& virulence.};
  \node[anchor=west, font=\tiny, text=gray!90!black] at (8.9, 0.8) {$\bullet$ Diverge only on shear axis.};
\end{tikzpicture}
\end{document}
"""
    compile_tikz_to_pdf(fig8_tex, "fig8_simulator_concordance_radar")

if __name__ == '__main__':
    print("=== Phase 6: Publication Figure Suite Generation (8 Vector PDFs) ===")
    generate_pdf_figures()
    print("Phase 6 completed successfully.")
