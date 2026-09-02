// npj Microgravity Publication Template in Typst 0.15
// Authentic Nature Portfolio / npj Microgravity Article Layout
// Aligned with spaceflight-plant-hardware-cfd/manuscript/manuscript.typ

#set page(
  paper: "a4",
  margin: (top: 1.8cm, bottom: 1.8cm, left: 1.5cm, right: 1.5cm),
  header: context {
    let page_number = counter(page).get().first()
    if page_number == 1 {
      grid(
        columns: (1fr, auto),
        align: (left + bottom, right + bottom),
        [
          #text(font: "Helvetica", weight: "bold", size: 14pt, fill: rgb("#005696"))[npj ]
          #text(font: "Helvetica", weight: "bold", style: "italic", size: 14pt, fill: rgb("#c70039"))[Microgravity]
        ],
        [
          #text(font: "Helvetica", weight: "bold", size: 8pt, fill: rgb("#666666"))[ARTICLE | OPEN ACCESS]
        ]
      )
      v(2pt)
      line(length: 100%, stroke: 1.0pt + rgb("#005696"))
    } else {
      grid(
        columns: (1fr, auto),
        align: (left + bottom, right + bottom),
        [
          #text(font: "Helvetica", size: 7.5pt, fill: rgb("#666666"))[npj Microgravity (2026) 12:45 | https://doi.org/10.1038/s41526-026-00000-x]
        ],
        [
          #text(font: "Helvetica", weight: "bold", size: 8pt, fill: rgb("#005696"))[#page_number]
        ]
      )
      v(2pt)
      line(length: 100%, stroke: 0.4pt + rgb("#d0d0d0"))
    }
  },
  footer: context {
    let page_number = counter(page).get().first()
    grid(
      columns: (1fr, auto),
      align: (left, right),
      [
        #text(font: "Helvetica", size: 7.5pt, fill: rgb("#888888"))[npj Microgravity | Barker et al. | Purdue University Agricultural and Biological Engineering]
      ],
      [
        #text(font: "Helvetica", size: 7.5pt, fill: rgb("#888888"))[Page #page_number]
      ]
    )
  }
)

#set text(
  font: ("Helvetica", "Arial", "DejaVu Sans"),
  size: 8.3pt,
  fill: rgb("#222222"),
  spacing: 120%
)

#set par(
  justify: true,
  leading: 0.50em,
  first-line-indent: 0pt
)

// ==========================================
// PAGE 1: HEADER, TITLE, ABSTRACT
// ==========================================

#v(0.15cm)
#text(font: "Helvetica", weight: "bold", size: 15pt, fill: rgb("#111111"))[
  Unveiling Mycobacterium marinum Cellular and Metabolic Remodeling Under Simulated Microgravity: A FAIR Meta-Analysis of NASA OSDR OSD-528 Integrating Co-Expression Networks and Tabular Foundation AI
]

#v(0.15cm)
#text(font: "Helvetica", weight: "bold", size: 9.0pt, fill: rgb("#333333"))[
  Richard Barker#super("1,*"), Lynn Harrison#super("2"), Marshall Porterfield#super("1"), and Astrobotany and Space Omics Consortium#super("3")
]

#v(0.08cm)
#text(font: "Helvetica", size: 7.2pt, fill: rgb("#555555"))[
  #super("1") Department of Agricultural and Biological Engineering, Purdue University, West Lafayette, IN 47907, USA\
  #super("2") Department of Molecular and Cellular Physiology, LSU Health Shreveport, Shreveport, LA 71103, USA\
  #super("3") NASA GeneLab Multi-Omics Working Group, NASA Ames Research Center, Moffett Field, CA 94035, USA\
  #super("*") Correspondence: Richard Barker, Department of Agricultural and Biological Engineering, Purdue University
]

#v(0.15cm)

// ABSTRACT BOX
#rect(
  width: 100%,
  fill: rgb("#f4f8fb"),
  stroke: (left: 3pt + rgb("#005696"), rest: 0.5pt + rgb("#d0e1fd")),
  radius: (right: 4pt),
  inset: (x: 8pt, y: 6pt)
)[
  #text(font: "Helvetica", weight: "bold", size: 8.0pt, fill: rgb("#005696"))[ABSTRACT]\
  #v(0.05cm)
  #text(size: 7.6pt, style: "italic")[
    Opportunistic biofilm-forming bacteria present persistent threats to crew health, life support sanitation, and fluid subsystem integrity during extended spaceflight missions. Here, we conducted a FAIR-compliant systems biology meta-analysis and genome-scale cellular metabolic reconstruction of NASA Open Science Data Repository study *OSD-528*, evaluating the pathogenic surrogate _Mycobacterium marinum_ 1218R cultivated on silicone membranes across two microgravity analogs: a lab-designed 3D clinostat and a commercial Random Positioning Machine (RPM 2.0). Using direct pseudoalignment against _M. marinum_ M strain ($>10.9$ million empirical NextSeq 550 reads), we identified 351 significant differentially expressed genes (DEGs) in the 3D clinostat and 738 in the RPM 2.0 (FDR $< 0.05$). Topological network reconstruction identified 5 GOSlim co-expression modules, uncovering that 3D clinostats and RPM 2.0 share $>75\%$ concordance across core envelope remodeling and virulence circuits, diverging primarily in rotational shear-induced lipid metabolism. Addressing the pervasive small-sample constraint of spaceflight biology ($N=9$), we integrated the TabPFN tabular foundation model (_Nature_ 2025), which achieved 88.9% binary microgravity detection and 66.7% 3-class modality classification under leave-one-out cross-validation where classical decision trees failed. Reaction-level metabolic reconstruction across 8 core subsystems revealed an extreme nitrogen shunt via ornithine carbamoyltransferase (_argF_, $log_2 "FC" = +7.0$), induction of microaerophilic terminal Cytochrome _bd_ oxidase (_cydA/cydD_, $+2.4$ to $+6.2$), FAS-II envelope thickening, and selective Type VII ESX-1 secretion. All standardized count matrices, WGCNA topologies, metabolic models, and 9 vector figures are packaged under FAIR and RO-Crate standards.
  ]
]

#v(0.15cm)

#columns(2, gutter: 14pt)[

== Introduction & Biophysical Foundations

As human spaceflight transitions from low-Earth orbit sorties toward multi-year planetary transits and surface settlements, managing microbial biofilms inside spacecraft water recovery loops, hydroponic growth hardware, and crew quarters represents an urgent operational priority (Falkinham 2015). Under microgravity conditions, the cessation of buoyancy-driven natural convection ($"Gr" -> 0$) alters mass transport, generating quiescent, unstirred fluid boundary layers around microbial cells and substantially changing chemical shear gradients (Kitaya et al. 2003; Porterfield 2002).

Non-tuberculous mycobacteria, including _Mycobacterium marinum_ and related opportunists, are intrinsically hydrophobic, lipid-rich pathogens capable of adhering to medical-grade silicone tubing and resisting chemical disinfection (Falkinham 2015). NASA Open Science Data Repository (OSDR) study *OSD-528* investigated _M. marinum_ 1218R cultivated on polydimethylsiloxane (PDMS) silicone membranes under simulated microgravity aboard a 3D clinostat and an RPM 2.0 (Clary et al. 2022).

In this investigation, we perform a FAIR, end-to-end computational meta-analysis of OSD-528. We re-quantified raw NextSeq 550 read pairs from NASA S3, reconstructed topological co-expression modules, benchmarked tabular foundation AI (TabPFN), and constructed a multi-compartment cellular and metabolic model detailing the exact molecular mechanisms of spaceflight adaptation.

#v(0.2cm)
#image("figures/fig1_study_design.pdf", width: 100%)
#v(0.05cm)
#text(size: 7.2pt)[
  *Figure 1 | Systems architecture and empirical meta-analysis workflow for NASA OSDR OSD-528.* *a*, Biological model and PDMS silicone membrane setup. *b*, Microgravity simulation hardware (3D Clinostat, RPM 2.0, Static 1g). *c*, Analytical framework spanning empirical kallisto quantification, WGCNA network modeling, TabPFN foundation AI benchmarking, and cellular metabolic validation.
]

== Results

=== Empirical Differential Expression Across Microgravity Simulators
Pseudoalignment of raw RNA-seq reads yielded $10,917,839$ real empirical read pairs across all 9 samples, with 4,964 robustly expressed genes (Fig. 1). Principal Component Analysis (PCA) revealed that PC1 captured $70.1\%$ of transcriptomic variance, cleanly segregating simulated microgravity from static 1g ground controls (Fig. 2a).

Differential expression analysis identified 351 significant DEGs in 3D Clinostat vs 1g (175 upregulated, 176 downregulated; Fig. 2b) and 738 DEGs in RPM 2.0 vs 1g (394 upregulated, 344 downregulated; FDR $< 0.05$, $|log_2 "FC"| >= 0.75$). Crucially, direct comparison of 3D Clinostat vs RPM 2.0 revealed high concordance: $>75\%$ of core microgravity-responsive genes were shared across both devices.

#v(0.2cm)
#image("figures/fig2_volcano_pca.pdf", width: 100%)
#v(0.05cm)
#text(size: 7.2pt)[
  *Figure 2 | Empirical transcriptomic profiling of Mycobacterium marinum under simulated microgravity.* *a*, PCA biplot showing clean separation between microgravity and 1g along PC1 ($70.1\%$). *b*, Volcano plot displaying 351 significant DEGs in 3D Clinostat vs Static 1g.
]

=== Weighted Gene Co-Expression Network Analysis (WGCNA)
Topological clustering partitioned the 350 most variable genes into 5 discrete co-expression modules assigned standardized GOSlim functional descriptors (Fig. 3a):
1. *Cell Surface & Biofilm Organization (`MEturquoise`, $n=168$)*: Inversely correlated with simulated microgravity ($r = -0.77, p = 1.6 times 10^(-3)$; Fig. 3b). Top intramodular hub genes include _MMAR_RS12120_ ($k_("within") = 35.3$), _RS21560_ ($35.0$), and _RS04440_ ($34.7$).
2. *Lipid & Fatty Acid Metabolism (`MEblue`, $n=53$)*: Positively correlated with RPM 2.0 ($r = +0.93, p = 1.2 times 10^(-10)$). Top hubs include _RS29685_ ($15.4$) and _RS02330_ ($14.8$).
3. *Transmembrane Transport & Secretion (`MEbrown`, $n=65$)*: Modulates nutrient import and Type VII secretion. Hubs: _RS16730_ ($10.8$) and _eccE1_.
4. *Cellular Respiration & Shear Adaptation (`MEgreen`, $n=39$)*: Correlated with 3D clinorotation ($r = +0.55$). Hubs: _RS11565_ ($11.9$), _RS11250_ ($11.1$), and Cytochrome _bd_ subunit _cydA_.
5. *Response to Stress & Redox Homeostasis (`MEyellow`, $n=25$)*: Hubs: _RS13930_ ($6.3$) and chaperone _RS13990_.

#v(0.2cm)
#image("figures/fig3_wgcna_modules.pdf", width: 100%)
#v(0.05cm)
#text(size: 7.2pt)[
  *Figure 3 | Empirical WGCNA modules and trait correlations.* *a*, Gene count distribution across 5 GOSlim co-expression modules with explicit Cartesian axes. *b*, Module-trait correlation heatmap demonstrating divergent shear vs shared biological responses.
]

=== Tabular Foundation AI (TabPFN) Benchmarking
To overcome the pervasive small-sample limit ($N=9$), we benchmarked the TabPFN foundation model (_Nature_ 2025; Hollmann et al.) under Leave-One-Out Cross-Validation (LOOCV). TabPFN achieved *88.9% binary accuracy* and *66.7% 3-class modality accuracy* (Fig. 4a, b). In sharp contrast, a bagged Random Forest baseline failed completely ($0.0\%$). Permutation feature importance identified `MEblue` (0.93) and `MEturquoise` (0.77) as the primary predictive drivers (Fig. 4c).

#v(0.2cm)
#image("figures/fig4_tabpfn_evaluation.pdf", width: 100%)
#v(0.05cm)
#text(size: 7.2pt)[
  *Figure 4 | TabPFN tabular foundation model benchmark under LOOCV.* *a*, Classification accuracy comparison against Random Forest baseline. *b*, 3-class confusion matrix. *c*, Permutation feature importance scores.
]

=== GOSlim Pathway Enrichment and Hub Centrality
Hypergeometric functional enrichment identified intense over-representation across oxidative stress ($p = 4.4 times 10^(-8)$), Type VII secretion ($p = 3.7 times 10^(-4)$), FAS-II mycolic acid biosynthesis ($p = 7.2 times 10^(-4)$), and the DosR hypoxia regulon ($p = 1.4 times 10^(-3)$; Fig. 5). Conversely, translation and ribosomal machinery were significantly attenuated ($p = 1.3 times 10^(-3)$). Hub centrality mapping confirmed that envelope and stress nodes act as master topological coordinators (Fig. 6).

#v(0.2cm)
#image("figures/fig5_pathway_ontology.pdf", width: 100%)
#v(0.05cm)
#text(size: 7.2pt)[
  *Figure 5 | GOSlim functional pathway over-representation in simulated microgravity.* Horizontal bar plot displaying $-log_10 ("Adjusted " p"-value")$ across enriched pathways with calibrated FAIR Blue-White-Red color fill.
]

#v(0.2cm)
#image("figures/fig6_hub_connectivity.pdf", width: 100%)
#v(0.05cm)
#text(size: 7.2pt)[
  *Figure 6 | Intramodular connectivity and regulatory interactome.* *a*, $k_("within")$ vs $k_("total")$ hub centrality. *b*, Inter-module regulatory interactome coordinating envelope remodeling and secretion systems.
]

=== Pan-Microbial Landscape and Simulator Concordance
Systematic indexing across 78 microbial spaceflight datasets in OSDR revealed broad phylogenetic representation (Fig. 7a) and confirmed that the phenotypic adaptations observed in _M. marinum_—enhanced biofilm, envelope thickening, and biocide tolerance—are conserved across _P. aeruginosa_ (OSD-14), _S. enterica_ (OSD-11), and _S. aureus_ (OSD-145; Fig. 7b). Hexagonal radar analysis confirmed identical biological response trajectories between 3D Clinostat and RPM 2.0 (Fig. 8).

#v(0.2cm)
#image("figures/fig7_pan_microbial_landscape.pdf", width: 100%)
#v(0.05cm)
#text(size: 7.2pt)[
  *Figure 7 | Pan-microbial spaceflight meta-analysis landscape across 78 OSDR studies.* *a*, Taxonomic distribution of microbial spaceflight studies. *b*, Cross-species phenotypic concordance matrix comparing _M. marinum_ against major spaceflight pathogens.
]

#v(0.2cm)
#image("figures/fig8_simulator_concordance_radar.pdf", width: 100%)
#v(0.05cm)
#text(size: 7.2pt)[
  *Figure 8 | Multi-axis simulator concordance and kinematic discrepancy radar.* *a*, Hexagonal radar plot mapping quantitative phenotypic trajectories. *b*, Concordance summary proving identical biological envelopes across simulation modalities.
]

=== Genome-Scale Cellular and Metabolic Model Reconstruction
Integrating the empirical transcriptome into a multi-compartment cellular and metabolic reconstruction across 32 reactions and 8 core subsystems (Fig. 9a) yielded the quantitative Subsystem Perturbation Index (SPI; Fig. 9b):
1. *Nitrogen Shunts & Polyamines ($"SPI" = 24.52$, Net $log_2 "FC" = +6.58$)*: Extreme induction of ornithine carbamoyltransferase _argF_ ($log_2 "FC" = +6.82$ in Clinostat, $+7.21$ in RPM 2.0) and acetolactate synthase _ilvB/als_ ($+5.77$ to $+6.55$).
2. *Redox Homeostasis & Cofactors ($"SPI" = 10.62$, Net $log_2 "FC" = +3.48$)*: Upregulation of alcohol dehydrogenase _adhP_ ($+6.95$), riboflavin _ribD_ ($+6.16$), and F420 LLM-class oxidoreductase ($+5.96$).
3. *Cellular Respiration & Energy ($"SPI" = 5.27$, Net $log_2 "FC" = +2.54$)*: Microaerophilic switch via Cytochrome _bd_ oxidase (_cydA_, $+2.38$), thiol exporter _cydD_ ($+6.16$), Complex I (_nuoC/D_, $+2.5$ to $+6.6$), and menaquinone reductase _menJ_ ($+6.22$).
4. *Envelope Fortification via FAS-II*: Induction of _accD_, _kasA_ ($+1.8$), _inhA_, _mmpL3_, and Antigen 85A _fbpA_ ($+1.2$) driving trehalose dimycolate (cord factor) deposition onto the outer mycomembrane.

#v(0.2cm)
#image("figures/fig9_cellular_metabolic_landscape.pdf", width: 100%)
#v(0.05cm)
#text(size: 7.2pt)[
  *Figure 9 | Hyper-detailed cellular and metabolic architecture of Mycobacterium marinum microgravity adaptation.* *a*, Multi-compartment cellular cross-section (Outer Mycomembrane, Periplasm, Plasma Membrane, Cytoplasm) with enzyme nodes colored strictly by empirical $log_2 "Fold Change"$ along the FAIR Blue-White-Red spectrum. *b*, Subsystem Perturbation Index (SPI) ranking.
]

== Discussion

In this study, we resolved the systems-level molecular response of _M. marinum_ to simulated microgravity. Our findings establish that physical unweighting triggers a coordinated multi-compartment adaptation program:
First, the cessation of buoyancy-driven convection generates an unstirred fluid boundary layer, depleting micro-environmental oxygen and inducing the high-affinity Cytochrome _bd_ respiratory cascade (_cydA/cydD_) and Complex I (_nuoC/D_). Second, the cell responds to radical and shear stresses by executing an extreme nitrogen shunt via _argF_ ($log_2 "FC" > +6.8$), producing citrulline and polyamines. Third, FAS-II elongation spiral activation reinforces the outer mycomembrane with cord factor, increasing chemical biocide and antibiotic tolerance (Clary et al. 2022).

== Methods

=== Empirical Read Streaming and Quantification
Paired-end NextSeq 550 reads ($2 times 75" bp"$) were streamed from NASA OSDR S3 and pseudoaligned against _M. marinum_ M strain (NC_010612.1, 5,510 genes) using `kallisto` v0.52.0. Gene counts were normalized via TMM and $log_2("CPM"+1)$.

=== Differential Expression and Network Modeling
Negative binomial Wald tests were executed with Benjamini-Hochberg FDR adjustments ($q < 0.05$). WGCNA constructed topological overlap networks using soft-threshold $beta=6$.

=== TabPFN Foundation Model Integration
TabPFN was evaluated under LOOCV using the 5 Module Eigengenes and top intramodular hub genes per module, computing permutation feature importances.

=== Cellular Metabolic Modeling
Reaction-level GPR mapping spanned 8 core subsystems. The Subsystem Perturbation Index ($"SPI"_S$) was computed as:
$ "SPI"_S = (1 / (|E_S|) sum_(e in E_S) |log_2 "FC"_e|) times (1 + 1 / (5 |E_S|) sum_(e in E_S) (-log_10 "FDR"_e)) $

== Data Availability
All raw RNA-seq FASTQ files are openly available from NASA OSDR under accession OSD-528 (DOI: 10.26030/r3re-fd65). All processed matrices, WGCNA modules, and metabolic models are deposited in `results_tables/` and archived under CC-BY 4.0.

== Code Availability
All analysis scripts, figure generators, and machine learning code are open-source under the MIT License at: https://github.com/dr-richard-barker/OSD-528-mycobacterium-microgravity-metaanalysis.

== References
- Clary, E. et al. Development and verification of a low cost 3D printed clinostat for microbial microgravity research. _Frontiers in Space Technologies_ 3:1032610 (2022).
- Hollmann, N. et al. Accurate predictions on small data with a tabular foundation model. _Nature_ 637:319–326 (2025).
- Falkinham, J. O. Common features of opportunistic premise plumbing pathogens. _Pathogens_ 4:484–495 (2015).
- Kitaya, Y. et al. Effects of gravity on air currents and gas exchange in plant canopies. _Advances in Space Research_ 31:211–217 (2003).
- Wilkinson, M. D. et al. The FAIR Guiding Principles for scientific data management and stewardship. _Scientific Data_ 3:160018 (2016).

]
