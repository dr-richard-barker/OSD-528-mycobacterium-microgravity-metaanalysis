#!/usr/bin/env python3
"""
07_cellular_metabolic_modeling.py

KEGG-anchored reaction and subsystem summary of the OSD-528 DESeq2 results
(Mycobacterium marinum strain M; 3D clinostat, RPM 2.0 and static 1g).

What this is, and what it is not
--------------------------------
A small, hand-curated set of reactions and named proteins grouped into eight
physiological subsystems. Every reaction is linked to M. marinum M genes only
through KEGG (organism code mmi, T00694). Each reaction and subsystem is then
summarised with statistics computed directly from the DESeq2 tables. It is not
a genome-scale model. It runs no flux balance analysis and makes no flux
predictions. It also asserts no gene-protein-reaction (AND/OR) logic, because a
KEGG EC or KO link does not say whether linked genes are subunits of one
complex or alternative isozymes.

This version replaces an earlier reaction list in which most of the 32 curated
reactions carried the wrong M. marinum locus. The superseded assignments are
kept in OLD_ASSIGNMENTS below only so that every one of them can be audited
(metabolic_model_locus_audit.tsv). None of them is used for any result.

Gene mapping (KEGG)
-------------------
1. Each reaction lists one or more KEGG queries:
     ec:<EC>  when the reaction has a complete, current EC number
              -> https://rest.kegg.jp/link/mmi/ec:<EC>
     ko:<KO>  when the entry is a named protein with no EC reaction (transporter,
              secretion-system component, chaperone, regulator), or when KEGG
              links no mmi gene to the EC and puts the named enzyme in a KO
              -> https://rest.kegg.jp/link/mmi/ko:<KO>
   The reaction's genes are the union of the mmi genes those queries return.
   Nothing is added by hand. A reaction whose queries return no mmi gene, or
   that has no valid query, is dropped and listed.
2. KEGG gene records (get/mmi:MMAR_xxxx), ENZYME records (get/ec:) and KO
   records (get/ko:) give names, orthology, positions and reaction equations.
3. Identifier bridge. KEGG names mmi genes by the original Sanger locus tags
   (MMAR_xxxx). The DESeq2 tables use RefSeq PGAP locus tags (MMAR_RSxxxxx) of
   assembly GCF_000018345.1. The two are joined through the /old_locus_tag
   qualifiers of the RefSeq GBFF, and every pair is checked by comparing the
   KEGG gene position (list/mmi) with the RefSeq gene coordinates (same
   replicon and strand, overlapping span). If the DESeq2 tables use MMAR_xxxx
   ids instead, the KEGG ids are used directly.

Expression statistics (all computed from the DESeq2 tables)
-----------------------------------------------------------
Contrasts: 3D clinostat vs static 1g, and RPM 2.0 vs static 1g (the two
simulator-vs-ground contrasts; the earlier version of this script did not use
clinostat-vs-RPM, so neither does this one).

  tested gene   a gene present in the contrast's DESeq2 table with non-missing padj
  DEG           significant == "YES" in the DESeq2 table (BH-adjusted Wald p < 0.05,
                as defined by 02_differential_expression.R)
  up / down     DEG with log2FoldChange > 0 / < 0

Reaction level, per contrast, over the reaction's tested KEGG genes:
  n_tested, n_deg, n_up, n_down
  lead gene      the tested gene with the smallest padj (ties: smaller pvalue,
                 then gene id), reported with its padj and its apeglm-shrunken
                 log2 fold change
  median_lfc     median of log2FoldChange_shrunk (apeglm) over the tested genes.
                 Shrunken values are used so that low-count genes with large,
                 unstable MLE fold changes do not dominate the summary.

Subsystem level, per contrast. The subsystem's genes are the de-duplicated
union of its reactions' KEGG genes.
  n_genes_mapped, n_tested, n_deg, n_up, n_down, median_lfc, lead gene (as above)
  deg_fraction       n_deg / n_tested
  fold_enrichment    (n_deg / n_tested) / (all DEGs / all tested genes)
  hypergeom_p        one-sided hypergeometric test for DEG over-representation,
                     P(X >= n_deg), X ~ Hypergeom(N = all tested genes in the
                     contrast, K = all DEGs in the contrast, n = n_tested)
  bh_fdr             Benjamini-Hochberg FDR over the subsystems, within each contrast
No weights, scaling constants or composite indices are used. The earlier
"Subsystem Perturbation Index" (mean |log2FC| x (1 + mean -log10 padj / 5)),
the +/-0.5 "direction" calls on the mean of the two simulators, and the
"predominant flux shift" labels are removed.

Reproducibility
---------------
Raw KEGG REST responses are cached in data/reference/kegg/ (git-ignored, like
the rest of data/reference/). The parts this script uses are written to a
tracked snapshot, data/processed/kegg_metabolic_snapshot.json, with the KEGG
release and retrieval dates. By default the script runs from that snapshot, so
a fresh clone reproduces the outputs with no network access and no GBFF. KEGG
requests are spaced at least KEGG_MIN_INTERVAL_S apart, batched where the API
allows, and sent with a plain User-Agent and no contact details. The outputs
contain no run timestamps.

Usage
-----
  python3 analysis/07_cellular_metabolic_modeling.py
      Run from the snapshot. Build the snapshot first if it is missing or does
      not cover the current reaction list.
  python3 analysis/07_cellular_metabolic_modeling.py --rebuild-kegg
      Rebuild the snapshot from the cached KEGG responses, downloading only
      what is missing. Needs the RefSeq GBFF in data/reference/.
  python3 analysis/07_cellular_metabolic_modeling.py --refresh-kegg
      Re-download every KEGG response, then rebuild the snapshot.

Outputs (data/processed/)
-------------------------
  metabolic_model_reactions.tsv        one row per retained reaction, with its KEGG
                                       genes and per-contrast summaries (Table S14)
  metabolic_subsystem_perturbation.tsv one row per subsystem x contrast (Table S15)
  metabolic_model_sbml.json            the model as JSON: definitions, reactions,
                                       KEGG gene links, statistics, audit (Table S17).
                                       The file name is historical; it is not SBML.
  metabolic_model_gene_expression.tsv  one row per reaction x KEGG gene, with the
                                       DESeq2 values for both contrasts
  metabolic_model_locus_audit.tsv      one row per superseded assignment (32): the old
                                       locus, what it really is, the KEGG genes, verdict
  kegg_metabolic_snapshot.json         the KEGG data used (see Reproducibility)
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import math
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PROCESSED = os.path.join(PROJECT_DIR, "data", "processed")
DATA_REFERENCE = os.path.join(PROJECT_DIR, "data", "reference")
KEGG_CACHE = os.path.join(DATA_REFERENCE, "kegg")
GBFF_FILE = os.path.join(DATA_REFERENCE, "GCF_000018345.1_ASM1834v1_genomic.gbff.gz")

SNAPSHOT = os.path.join(DATA_PROCESSED, "kegg_metabolic_snapshot.json")
OUT_REACTIONS = os.path.join(DATA_PROCESSED, "metabolic_model_reactions.tsv")
OUT_SUBSYSTEMS = os.path.join(DATA_PROCESSED, "metabolic_subsystem_perturbation.tsv")
OUT_JSON = os.path.join(DATA_PROCESSED, "metabolic_model_sbml.json")
OUT_GENES = os.path.join(DATA_PROCESSED, "metabolic_model_gene_expression.tsv")
OUT_AUDIT = os.path.join(DATA_PROCESSED, "metabolic_model_locus_audit.tsv")

# ---------------------------------------------------------------- parameters
USER_AGENT = "osd528-reanalysis/1.0"
KEGG_BASE = "https://rest.kegg.jp/"
KEGG_ORG = "mmi"
KEGG_MIN_INTERVAL_S = 1.0
KEGG_GET_BATCH = 10          # KEGG accepts at most 10 entries per get request
KEGG_RETRIES = 3
CHROMOSOME_ACC = "NC_010612.1"

CONTRASTS = [
    # (column prefix, contrast id, DESeq2 file, label)
    ("clinostat", "clinostat_vs_static1g", "deg_3dclinostat_vs_static1g.tsv", "3D clinostat vs static 1g"),
    ("rpm", "rpm2_vs_static1g", "deg_rpm2_vs_static1g.tsv", "RPM 2.0 vs static 1g"),
]

# --------------------------------------------------------- curated reactions
# Curated here: the reaction, its subsystem and compartment, and which KEGG
# identifier(s) define it. The genes come only from KEGG. "replaces" lists the
# superseded OLD_ASSIGNMENTS entries; "scope" says how specific the KEGG link is.
SUBSYSTEM_RENAMES = {
    "GPL Biofilm & Cell Surface": (
        "Peptidoglycan Biosynthesis",
        "KEGG has no mmi gene in the glycopeptidolipid KOs (mps1 K16228, mps2 K16229), "
        "so both GPL reactions are dropped; only MurA remains."),
    "Nitrogen Shunts & Polyamines": (
        "Arginine & Branched-Chain Amino Acid Biosynthesis",
        "The two reactions are ornithine carbamoyltransferase (arginine biosynthesis) and "
        "acetolactate synthase (valine/leucine/isoleucine); neither is a polyamine reaction."),
    "Type VII Secretion & Virulence": (
        "Type VII (ESX) Secretion",
        "The ESX KOs group the paralogues of every ESX locus, so these gene sets are not "
        "specific to the ESX-1 virulence system."),
}

S_RESP = "Cellular Respiration & Energy"
S_CARB = "Central Carbon & Glyoxylate Shunt"
S_MYCO = "Mycolic Acid & FAS-II Envelope"
S_PGN = "Peptidoglycan Biosynthesis"
S_ESX = "Type VII (ESX) Secretion"
S_AMINO = "Arginine & Branched-Chain Amino Acid Biosynthesis"
S_REDOX = "Redox Homeostasis & Cofactors"
S_DOSR = "DosR Hypoxia & Dormancy"
SUBSYSTEM_ORDER = [S_RESP, S_CARB, S_MYCO, S_PGN, S_ESX, S_AMINO, S_REDOX, S_DOSR]

REACTIONS = [
    # --- Cellular respiration and energy
    dict(reaction_id="R_NDH1", subsystem=S_RESP, compartment="Plasma Membrane",
         reaction_name="NADH:quinone oxidoreductase (Complex I, NDH-1)",
         kegg_queries=["ec:7.1.1.2"], scope="complex: every subunit KEGG links to the EC",
         replaces=["R_NDH1_C", "R_NDH1_D"],
         note="The two old entries (subunits C and D) are one reaction, EC 7.1.1.2; merged."),
    dict(reaction_id="R_MENJ", subsystem=S_RESP, compartment="Plasma Membrane",
         reaction_name="Menaquinone-9 beta-reductase (MenJ)",
         kegg_queries=["ec:1.3.99.38"], scope="specific", replaces=["R_MENJ"],
         note="Old EC 1.3.5.2 is dihydroorotate dehydrogenase (quinone) in KEGG; MenJ is EC 1.3.99.38."),
    dict(reaction_id="R_CYDA", subsystem=S_RESP, compartment="Plasma Membrane",
         reaction_name="Cytochrome bd quinol oxidase",
         kegg_queries=["ec:7.1.1.7"], scope="complex: every subunit KEGG links to the EC",
         replaces=["R_CYDA"], note="KEGG lists old EC 1.10.3.14 as transferred to 7.1.1.7."),
    dict(reaction_id="R_CYDD", subsystem=S_RESP, compartment="Plasma Membrane",
         reaction_name="Thiol reductant ABC exporter subunit CydD",
         kegg_queries=["ko:K16013"], scope="specific", replaces=["R_CYDD"],
         note="Old EC 7.6.2.2 (ABC-type xenobiotic transporter) links no mmi gene; CydD is KO K16013 (no EC)."),
    dict(reaction_id="R_ATPS", subsystem=S_RESP, compartment="Plasma Membrane",
         reaction_name="F1Fo ATP synthase (catalytic subunits)",
         kegg_queries=["ec:7.1.2.2"], scope="complex: only the subunits KEGG links to the EC",
         replaces=["R_ATPS"],
         note="KEGG links only the catalytic alpha and beta subunits to EC 7.1.2.2."),
    # --- Central carbon, TCA and glyoxylate shunt
    dict(reaction_id="R_ICL", subsystem=S_CARB, compartment="Cytoplasm",
         reaction_name="Isocitrate lyase (glyoxylate shunt)",
         kegg_queries=["ec:4.1.3.1"], scope="isozymes or paralogues linked to the EC",
         replaces=["R_ICL1"], note=""),
    dict(reaction_id="R_GLCB", subsystem=S_CARB, compartment="Cytoplasm",
         reaction_name="Malate synthase G",
         kegg_queries=["ec:2.3.3.9"], scope="specific", replaces=["R_GLCB"], note=""),
    dict(reaction_id="R_CS", subsystem=S_CARB, compartment="Cytoplasm",
         reaction_name="Citrate synthase",
         kegg_queries=["ec:2.3.3.1"], scope="isozymes or paralogues linked to the EC",
         replaces=["R_CS"], note=""),
    dict(reaction_id="R_SUCD", subsystem=S_CARB, compartment="Cytoplasm",
         reaction_name="Succinyl-CoA synthetase (ADP-forming)",
         kegg_queries=["ec:6.2.1.5"], scope="complex: every subunit KEGG links to the EC",
         replaces=["R_SUCD"], note="The old entry named only the alpha subunit."),
    # --- Mycolic acid and FAS-II envelope
    dict(reaction_id="R_FAS", subsystem=S_MYCO, compartment="Cytoplasm",
         reaction_name="Fatty acid synthase I (bacterial type FAS-I)",
         kegg_queries=["ko:K11533"], scope="specific", replaces=["R_FAS"],
         note="Old EC 2.3.1.85 (animal-type FAS) links no mmi gene; bacterial FAS-I is KO K11533 (EC 2.3.1.-)."),
    dict(reaction_id="R_ACCD", subsystem=S_MYCO, compartment="Cytoplasm",
         reaction_name="Acetyl-CoA carboxylase (malonyl-CoA synthesis)",
         kegg_queries=["ec:6.4.1.2"], scope="complex: every subunit KEGG links to the EC",
         replaces=["R_ACCD"], note=""),
    dict(reaction_id="R_KASA", subsystem=S_MYCO, compartment="Cytoplasm",
         reaction_name="Meromycolate beta-ketoacyl-ACP synthase KasA (FAS-II)",
         kegg_queries=["ec:2.3.1.293"], scope="specific", replaces=["R_KASA"],
         note="Old EC 2.3.1.41 links no mmi gene; KEGG assigns KasA (K11609) to EC 2.3.1.293."),
    dict(reaction_id="R_INHA", subsystem=S_MYCO, compartment="Cytoplasm",
         reaction_name="Meromycolate enoyl-ACP reductase InhA (FAS-II)",
         kegg_queries=["ec:1.3.1.118"], scope="specific", replaces=["R_INHA"],
         note="Old EC 1.3.1.9 links no mmi gene; KEGG assigns InhA (K11611) to EC 1.3.1.118."),
    dict(reaction_id="R_MMPL3", subsystem=S_MYCO, compartment="Plasma Membrane",
         reaction_name="Trehalose monomycolate transporter MmpL3",
         kegg_queries=["ko:K20470"], scope="specific", replaces=["R_MMPL3"],
         note="Old '2.A.1' is a TCDB transporter class, not an EC number; MmpL3 is KO K20470."),
    dict(reaction_id="R_FBPA", subsystem=S_MYCO, compartment="Periplasm / Cell Wall",
         reaction_name="Antigen 85 trehalose mycolyltransferases",
         kegg_queries=["ec:2.3.1.122"], scope="isozymes or paralogues linked to the EC",
         replaces=["R_FBPA"], note=""),
    dict(reaction_id="R_FACL", subsystem=S_MYCO, compartment="Cytoplasm",
         reaction_name="Long-chain-fatty-acid--CoA ligase",
         kegg_queries=["ec:6.2.1.3"], scope="isozymes or paralogues linked to the EC",
         replaces=["R_FADD7"],
         note="KEGG gene MMAR_0319 is named fadD7 (GenBank) but has no KO and is not linked to "
              "EC 6.2.1.3, so the reaction carries the EC-linked ligases, not FadD7."),
    # --- GPL / peptidoglycan
    dict(reaction_id="R_MPS1", subsystem=S_PGN, compartment="Cytoplasm / Inner Membrane",
         reaction_name="Glycopeptidolipid synthetase Mps1",
         kegg_queries=["ko:K16228"], scope="specific", replaces=["R_MPS1"],
         note="Old EC 6.3.2.- is incomplete, so the entry is mapped by the mps1 KO (K16228)."),
    dict(reaction_id="R_MPS2", subsystem=S_PGN, compartment="Cytoplasm / Inner Membrane",
         reaction_name="Glycopeptidolipid synthetase Mps2",
         kegg_queries=["ko:K16229"], scope="specific", replaces=["R_MPS2"],
         note="Old EC 6.3.2.- is incomplete, so the entry is mapped by the mps2 KO (K16229)."),
    dict(reaction_id="R_MURA", subsystem=S_PGN, compartment="Cytoplasm",
         reaction_name="UDP-N-acetylglucosamine 1-carboxyvinyltransferase MurA",
         kegg_queries=["ec:2.5.1.7"], scope="specific", replaces=["R_MURA"], note=""),
    # --- Type VII (ESX) secretion
    dict(reaction_id="R_ECCE", subsystem=S_ESX, compartment="Plasma Membrane",
         reaction_name="ESX secretion system protein EccE (all ESX loci)",
         kegg_queries=["ko:K27088"], scope="family: KO spans the ESX paralogues; ESX-1 not resolved",
         replaces=["R_ECCE1"],
         note="Old EC 3.6.3.- is an incomplete, obsolete class and EccE is not an ATPase; mapped by KO."),
    dict(reaction_id="R_ECCB", subsystem=S_ESX, compartment="Plasma Membrane",
         reaction_name="ESX secretion system protein EccB (all ESX loci)",
         kegg_queries=["ko:K27085"], scope="family: KO spans the ESX paralogues; ESX-1 not resolved",
         replaces=["R_ECCB1"], note="Old EC 3.6.3.- is an incomplete, obsolete class; mapped by KO."),
    dict(reaction_id="R_ESXA", subsystem=S_ESX, compartment="Extracellular / Capsule",
         reaction_name="ESAT-6 family secreted proteins (EsxA and paralogues)",
         kegg_queries=["ko:K14956"], scope="family: one KO covers every ESAT-6/CFP-10 family protein",
         replaces=["R_ESXA"],
         note="KEGG has no EsxA-specific KO; K14956 groups all ESAT-6 family proteins."),
    dict(reaction_id="R_ESPB", subsystem=S_ESX, compartment="Extracellular / Capsule",
         reaction_name="ESX secretion-associated protein EspB family",
         kegg_queries=["ko:K27099"], scope="family: KO spans EspB paralogues",
         replaces=["R_ESPB"], note=""),
    # --- Amino acid biosynthesis
    dict(reaction_id="R_ARGF", subsystem=S_AMINO, compartment="Cytoplasm",
         reaction_name="Ornithine carbamoyltransferase ArgF",
         kegg_queries=["ec:2.1.3.3"], scope="specific", replaces=["R_ARGF"], note=""),
    dict(reaction_id="R_ALS", subsystem=S_AMINO, compartment="Cytoplasm",
         reaction_name="Acetolactate synthase (branched-chain amino acids)",
         kegg_queries=["ec:2.2.1.6"], scope="complex or isozymes linked to the EC",
         replaces=["R_ALS"],
         note="The old locus (MMAR_2697) is named acetolactate synthase IlvG in GenBank, but KEGG gives "
              "it no KO or EC, so it is not among the EC 2.2.1.6 genes."),
    # --- Redox and cofactors
    dict(reaction_id="R_ADHP", subsystem=S_REDOX, compartment="Cytoplasm",
         reaction_name="Alcohol dehydrogenase",
         kegg_queries=["ec:1.1.1.1"], scope="isozymes or paralogues linked to the EC",
         replaces=["R_ADHP"], note=""),
    dict(reaction_id="R_RIBD", subsystem=S_REDOX, compartment="Cytoplasm",
         reaction_name="Bifunctional riboflavin biosynthesis deaminase/reductase (RibD/RibG)",
         kegg_queries=["ec:3.5.4.26", "ec:1.1.1.193"], scope="specific", replaces=["R_RIBD"], note=""),
    dict(reaction_id="R_F420", subsystem=S_REDOX, compartment="Cytoplasm",
         reaction_name="F420-dependent LLM class oxidoreductase",
         kegg_queries=[], scope="none", replaces=["R_F420"],
         note="Old EC 1.5.-.- is incomplete and the old locus has no KEGG KO, so KEGG offers no "
              "reaction or ortholog group to map."),
    dict(reaction_id="R_KATG", subsystem=S_REDOX, compartment="Cytoplasm / Periplasm",
         reaction_name="Catalase-peroxidase KatG",
         kegg_queries=["ec:1.11.1.21"], scope="specific", replaces=["R_KATG"], note=""),
    # --- DosR hypoxia and dormancy
    dict(reaction_id="R_HSPX", subsystem=S_DOSR, compartment="Cytoplasm / Cell Wall",
         reaction_name="HSP20 family (alpha-crystallin) small heat shock proteins",
         kegg_queries=["ko:K13993"], scope="family: KO covers all HSP20 family proteins",
         replaces=["R_HSPX"],
         note="KEGG has no HspX-specific KO. GenBank-named hspX (MMAR_0912) is HtpX (K03799) in KEGG; "
              "hspX_1 (MMAR_3484) has an HSP20 domain but no KO, so it is not included."),
    dict(reaction_id="R_DOSR", subsystem=S_DOSR, compartment="Cytoplasm",
         reaction_name="DevR (DosR) two-component response regulator",
         kegg_queries=["ko:K07695"], scope="specific", replaces=["R_DOSR"], note=""),
]

# ------------------------------------------- superseded assignments (audited)
# The earlier version's reaction list, kept only so that each assignment can be
# checked against KEGG. None of these loci is used for any result.
# (old reaction id, old subsystem, old reaction name, old EC, old gene id, old symbol)
OLD_ASSIGNMENTS = [
    ("R_NDH1_C", "Cellular Respiration & Energy", "NADH-quinone oxidoreductase (Complex I)", "7.1.1.2", "MMAR_RS07320", "nuoC"),
    ("R_NDH1_D", "Cellular Respiration & Energy", "NADH dehydrogenase subunit D (Complex I)", "7.1.1.2", "MMAR_RS07325", "nuoD"),
    ("R_MENJ", "Cellular Respiration & Energy", "Menaquinone reductase", "1.3.5.2", "MMAR_RS04470", "menJ"),
    ("R_CYDA", "Cellular Respiration & Energy", "Cytochrome bd microaerophilic terminal oxidase", "1.10.3.14", "MMAR_RS12095", "cydA"),
    ("R_CYDD", "Cellular Respiration & Energy", "Thiol reductant ABC exporter CydD (Cytochrome bd assembly)", "7.6.2.2", "MMAR_RS12080", "cydD"),
    ("R_ATPS", "Cellular Respiration & Energy", "ATP synthase F1F0", "7.1.2.2", "MMAR_RS07455", "atpA"),
    ("R_ICL1", "Central Carbon & Glyoxylate Shunt", "Isocitrate lyase (Glyoxylate bypass)", "4.1.3.1", "MMAR_RS05335", "icl1"),
    ("R_GLCB", "Central Carbon & Glyoxylate Shunt", "Malate synthase G", "2.3.3.9", "MMAR_RS13410", "glcB"),
    ("R_CS", "Central Carbon & Glyoxylate Shunt", "Citrate synthase", "2.3.3.1", "MMAR_RS04685", "gltA"),
    ("R_SUCD", "Central Carbon & Glyoxylate Shunt", "Succinyl-CoA synthetase subunit alpha", "6.2.1.5", "MMAR_RS13725", "sucD"),
    ("R_FAS", "Mycolic Acid & FAS-II Envelope", "Fatty acid synthase I (FAS-I de novo synthesis)", "2.3.1.85", "MMAR_RS13525", "fas"),
    ("R_ACCD", "Mycolic Acid & FAS-II Envelope", "Acyl-CoA carboxylase subunit beta (Malonyl-CoA synthesis)", "6.4.1.2", "MMAR_RS22330", "accD"),
    ("R_KASA", "Mycolic Acid & FAS-II Envelope", "Beta-ketoacyl-ACP synthase KasA (FAS-II initiation)", "2.3.1.41", "MMAR_RS10900", "kasA"),
    ("R_INHA", "Mycolic Acid & FAS-II Envelope", "Enoyl-ACP reductase InhA (FAS-II elongation)", "1.3.1.9", "MMAR_RS14560", "inhA"),
    ("R_MMPL3", "Mycolic Acid & FAS-II Envelope", "Trehalose monomycolate translocase MmpL3", "2.A.1", "MMAR_RS01245", "mmpL3"),
    ("R_FBPA", "Mycolic Acid & FAS-II Envelope", "Antigen 85A mycolyltransferase (Cord factor TDM synthesis)", "2.3.1.122", "MMAR_RS00720", "fbpA"),
    ("R_FADD7", "Mycolic Acid & FAS-II Envelope", "Long-chain fatty-acid-CoA ligase FadD7", "6.2.1.3", "MMAR_RS12495", "fadD7"),
    ("R_MPS1", "GPL Biofilm & Cell Surface", "GPL core non-ribosomal peptide synthetase Mps1", "6.3.2.-", "MMAR_RS18395", "mps1"),
    ("R_MPS2", "GPL Biofilm & Cell Surface", "GPL lipopeptide synthetase Mps2", "6.3.2.-", "MMAR_RS18400", "mps2"),
    ("R_MURA", "GPL Biofilm & Cell Surface", "UDP-N-acetylglucosamine 1-carboxyvinyltransferase (Peptidoglycan)", "2.5.1.7", "MMAR_RS06710", "murA"),
    ("R_ECCE1", "Type VII Secretion & Virulence", "Type VII ESX-1 pore complex subunit EccE1", "3.6.3.-", "MMAR_RS00995", "eccE1"),
    ("R_ECCB1", "Type VII Secretion & Virulence", "Type VII ESX-1 membrane core subunit EccB1", "3.6.3.-", "MMAR_RS00980", "eccB1"),
    ("R_ESXA", "Type VII Secretion & Virulence", "6 kDa early secretory antigenic target ESAT-6 (EsxA)", "N/A", "MMAR_RS01000", "esxA"),
    ("R_ESPB", "Type VII Secretion & Virulence", "ESX-1 secretion system-associated protein EspB", "N/A", "MMAR_RS00990", "espB"),
    ("R_ARGF", "Nitrogen Shunts & Polyamines", "Ornithine carbamoyltransferase ArgF", "2.1.3.3", "MMAR_RS12305", "argF"),
    ("R_ALS", "Nitrogen Shunts & Polyamines", "Acetolactate synthase (Branched-chain amino acids)", "2.2.1.6", "MMAR_RS13445", "ilvB"),
    ("R_ADHP", "Redox Homeostasis & Cofactors", "Alcohol dehydrogenase AdhP", "1.1.1.1", "MMAR_RS24640", "adhP"),
    ("R_RIBD", "Redox Homeostasis & Cofactors", "Diaminohydroxyphosphoribosylaminopyrimidine deaminase (RibD)", "3.5.4.26 / 1.1.1.193", "MMAR_RS11045", "ribD"),
    ("R_F420", "Redox Homeostasis & Cofactors", "F420-dependent LLM class oxidoreductase", "1.5.-.-", "MMAR_RS06790", "f420_red"),
    ("R_KATG", "Redox Homeostasis & Cofactors", "Catalase-peroxidase KatG", "1.11.1.21", "MMAR_RS09685", "katG"),
    ("R_HSPX", "DosR Hypoxia & Dormancy", "Alpha-crystallin family heat shock protein HspX (Acr)", "N/A", "MMAR_RS13905", "hspX"),
    ("R_DOSR", "DosR Hypoxia & Dormancy", "DosR response regulator (Hypoxic boundary layer switch)", "N/A", "MMAR_RS03930", "dosR"),
]

EC_COMPLETE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")


def rel(path: str) -> str:
    return os.path.relpath(path, PROJECT_DIR)


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def old_ec_list(old_ec: str) -> list[str]:
    """Complete EC numbers in an old EC field ('3.5.4.26 / 1.1.1.193' -> two)."""
    return [e for e in (p.strip() for p in old_ec.split("/")) if EC_COMPLETE.match(e)]


# ================================================================ KEGG client
class Kegg:
    """Cached, rate-limited KEGG REST client. One cache file per request or entry."""

    def __init__(self, cache_dir: str, refresh: bool = False):
        self.cache_dir = cache_dir
        self.refresh = refresh
        self.requests = 0
        self.used_files: set[str] = set()
        self._last = 0.0
        self._refreshed: set[str] = set()
        os.makedirs(cache_dir, exist_ok=True)

    def _cache_path(self, key: str) -> str:
        return os.path.join(self.cache_dir, re.sub(r"[^A-Za-z0-9._+-]", "_", key) + ".txt")

    def _fetch(self, path: str) -> str:
        url = KEGG_BASE + urllib.parse.quote(path, safe="/:+._-")
        for attempt in range(KEGG_RETRIES):
            wait = KEGG_MIN_INTERVAL_S - (time.time() - self._last)
            if wait > 0:
                time.sleep(wait)
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    text = r.read().decode("utf-8")
                self._last = time.time()
                self.requests += 1
                return text
            except urllib.error.HTTPError as e:
                self._last = time.time()
                self.requests += 1
                if e.code == 404:          # KEGG answers 404 when a query has no hits
                    return ""
                if e.code in (403, 429, 500, 502, 503, 504) and attempt + 1 < KEGG_RETRIES:
                    time.sleep(10 * (attempt + 1))
                    continue
                raise RuntimeError(f"KEGG request failed ({e.code}): {url}") from e
            except urllib.error.URLError as e:
                self._last = time.time()
                if attempt + 1 < KEGG_RETRIES:
                    time.sleep(10 * (attempt + 1))
                    continue
                raise RuntimeError(f"KEGG request failed: {url}: {e}") from e
        raise RuntimeError(f"KEGG request failed: {url}")

    def _need(self, fn: str) -> bool:
        return not os.path.exists(fn) or (self.refresh and fn not in self._refreshed)

    def text(self, path: str) -> str:
        fn = self._cache_path(path)
        if self._need(fn):
            txt = self._fetch(path)
            with open(fn, "w") as f:
                f.write(txt)
            self._refreshed.add(fn)
        self.used_files.add(fn)
        with open(fn) as f:
            return f.read()

    def link_genes(self, query: str) -> list[str]:
        """Genes of KEGG_ORG linked to an 'ec:..' or 'ko:..' query, as bare locus tags."""
        genes = []
        for line in self.text(f"link/{KEGG_ORG}/{query}").splitlines():
            parts = line.split("\t")
            if len(parts) == 2 and parts[1].startswith(KEGG_ORG + ":"):
                genes.append(parts[1].split(":", 1)[1])
        return sorted(set(genes))

    def get_entries(self, db: str, ids: list[str]) -> dict[str, str]:
        """Flat-file records for db:id entries; cached one file per entry, fetched in batches."""
        out, missing = {}, []
        for i in dict.fromkeys(ids):
            fn = self._cache_path(f"get/{db}:{i}")
            if self._need(fn):
                missing.append(i)
            else:
                self.used_files.add(fn)
                with open(fn) as f:
                    out[i] = f.read()
        for k in range(0, len(missing), KEGG_GET_BATCH):
            batch = missing[k:k + KEGG_GET_BATCH]
            txt = self._fetch("get/" + "+".join(f"{db}:{i}" for i in batch))
            found = {}
            for rec in txt.split("///"):
                rec = rec.strip("\n")
                if not rec.strip():
                    continue
                m = re.match(r"ENTRY\s+(?:EC\s+)?(\S+)", rec)
                if m:
                    found[m.group(1)] = rec + "\n///\n"
            for i in batch:
                fn = self._cache_path(f"get/{db}:{i}")
                rec = found.get(i, "")      # empty file = KEGG returned no entry
                with open(fn, "w") as f:
                    f.write(rec)
                self._refreshed.add(fn)
                self.used_files.add(fn)
                out[i] = rec
        return out


def parse_flat(record: str) -> dict[str, list[str]]:
    """KEGG flat file -> {FIELD: [lines]} with continuation lines folded in."""
    fields: dict[str, list[str]] = {}
    key = None
    for line in record.splitlines():
        if line.startswith("///"):
            break
        k = line[:12].strip()
        if k:
            key = k
            fields.setdefault(key, [])
        if key is not None:
            fields[key].append(line[12:].strip())
    return fields


def kegg_gene_info(record: str) -> dict:
    f = parse_flat(record)
    name = " ".join(f.get("NAME", []))
    return {
        "kegg_symbol": " ".join(f.get("SYMBOL", [])),
        "kegg_name": re.sub(r"^\(GenBank\)\s*", "", name),
        "kegg_orthology": "; ".join(re.sub(r"\s+", " ", o) for o in f.get("ORTHOLOGY", [])),
        "kegg_ko": sorted(set(re.findall(r"\bK\d{5}\b", " ".join(f.get("ORTHOLOGY", []))))),
        "kegg_position": " ".join(f.get("POSITION", [])),
        "kegg_record_found": bool(record.strip()),
    }


def kegg_enzyme_info(record: str) -> dict:
    f = parse_flat(record)
    names = [n.rstrip(";") for n in f.get("NAME", [])]
    reaction = re.sub(r"\s+", " ", " ".join(f.get("REACTION", []))).strip()
    obsolete = "Obsolete" in " ".join(f.get("ENTRY", [])) or any(n.startswith("Transferred to") for n in names)
    return {"kegg_name": names[0] if names else "", "kegg_equation": reaction,
            "obsolete": obsolete, "record_found": bool(record.strip())}


def kegg_ko_info(record: str) -> dict:
    f = parse_flat(record)
    sym = " ".join(f.get("SYMBOL", []))
    name = " ".join(f.get("NAME", []))
    return {"kegg_name": f"{sym}; {name}" if sym else name, "kegg_equation": "",
            "obsolete": False, "record_found": bool(record.strip())}


# ======================================================= RefSeq GBFF parsing
def parse_gbff(path: str) -> dict[str, dict]:
    """RefSeq locus_tag -> {old_locus_tag, gene, product, location, replicon, pseudo}."""
    loci: dict[str, dict] = {}
    skip = ("translation", "inference", "note", "db_xref", "experiment")
    acc = None
    in_features = False
    cur = None
    key = None

    def flush(feat):
        if not feat:
            return
        q = feat["q"]
        lt = q.get("locus_tag")
        if not lt:
            return
        d = loci.setdefault(lt, {"old_locus_tag": "", "gene": "", "product": "", "location": "",
                                 "replicon": "", "pseudo": False})
        if q.get("old_locus_tag"):
            d["old_locus_tag"] = q["old_locus_tag"]
        if q.get("gene") and not d["gene"]:
            d["gene"] = q["gene"]
        if q.get("product"):
            d["product"] = q["product"]
        if "pseudo" in q:
            d["pseudo"] = True
        if feat["type"] == "gene":
            d["location"] = feat["location"]
            d["replicon"] = acc

    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("VERSION"):
                acc = line.split()[1]
            elif line.startswith("FEATURES"):
                in_features = True
                continue
            elif line.startswith(("ORIGIN", "CONTIG", "//")):
                flush(cur)
                cur, in_features = None, False
                continue
            if not in_features or len(line) < 22:
                continue
            if line[5] != " ":
                flush(cur)
                cur = {"type": line[5:21].strip(), "location": line[21:].strip(), "q": {}, "in_loc": True}
                key = None
                continue
            s = line[21:].rstrip("\n")
            m = re.match(r'^/([A-Za-z_]+)(?:=(.*))?$', s)
            if m:
                cur["in_loc"] = False
                key = m.group(1)
                if key in skip or key.startswith("GO_"):
                    key = None          # not needed; its continuation lines are skipped too
                    continue
                cur["q"][key] = (m.group(2) or "").strip('"')
            elif cur["in_loc"]:
                cur["location"] += s.strip()
            elif key is not None:
                cur["q"][key] = (cur["q"][key] + " " + s.strip()).strip('"')
    return loci


def span(location: str) -> tuple[str, int, int] | None:
    nums = [int(x) for x in re.findall(r"\d+", location.split(":")[-1] if ":" in location else location)]
    if not nums:
        return None
    return ("-" if "complement" in location else "+", min(nums), max(nums))


def position_check(kegg_pos: str, refseq_loc: str, refseq_replicon: str) -> str:
    """Compare a KEGG gene position with RefSeq gene coordinates."""
    if not kegg_pos or not refseq_loc:
        return "not_compared"
    kegg_plasmid = ":" in kegg_pos and not kegg_pos.startswith(("complement", "join"))
    refseq_plasmid = refseq_replicon != CHROMOSOME_ACC
    a, b = span(kegg_pos), span(refseq_loc)
    if a is None or b is None or kegg_plasmid != refseq_plasmid:
        return "mismatch"
    if a == b:
        return "exact"
    if a[0] == b[0] and a[1] <= b[2] and b[1] <= a[2]:
        return "overlap_same_strand"
    return "mismatch"


# ============================================================ snapshot build
def needed_queries() -> tuple[list[str], list[str]]:
    """(reaction queries, extra EC queries for the old EC numbers)."""
    rq = list(dict.fromkeys(q for r in REACTIONS for q in r["kegg_queries"]))
    old_ecs = list(dict.fromkeys(f"ec:{e}" for a in OLD_ASSIGNMENTS for e in old_ec_list(a[3])))
    return rq, [q for q in old_ecs if q not in rq]


def build_snapshot(refresh: bool) -> dict:
    if not os.path.exists(GBFF_FILE):
        sys.exit(f"Error: {rel(GBFF_FILE)} is needed to build the KEGG snapshot "
                 "(00_quantify_real_osdr_rnaseq.py downloads it).")
    kegg = Kegg(KEGG_CACHE, refresh=refresh)
    print(f"  Building KEGG snapshot (cache: {rel(KEGG_CACHE)}; refresh={refresh})")

    # KEGG no longer prints a release number in info/kegg; record the update
    # dates of the databases used instead (and the release line if present).
    info_kegg = kegg.text("info/kegg").splitlines()
    release = next((re.sub(r"^\S+\s+", "", l).strip() for l in info_kegg if "Release" in l), "")
    db_dates = {p[0]: p[2] for p in (l.split() for l in info_kegg) if len(p) == 3 and p[0] in ("genes", "ko", "enzyme")}
    mmi_info = [" ".join(l.split()) for l in kegg.text(f"info/{KEGG_ORG}").splitlines() if l.strip()][:3]

    reaction_queries, old_ec_queries = needed_queries()
    links = {q: kegg.link_genes(q) for q in reaction_queries + old_ec_queries}

    # Query metadata: ENZYME and KO records.
    ec_recs = kegg.get_entries("ec", [q[3:] for q in links if q.startswith("ec:")])
    ko_recs = kegg.get_entries("ko", [q[3:] for q in links if q.startswith("ko:")])
    queries = {}
    for q in reaction_queries + old_ec_queries:
        info = kegg_enzyme_info(ec_recs[q[3:]]) if q.startswith("ec:") else kegg_ko_info(ko_recs[q[3:]])
        queries[q] = {**info, "genes": links[q], "used_by_reactions": [
            r["reaction_id"] for r in REACTIONS if q in r["kegg_queries"]]}

    # RefSeq annotation and the genome-wide identifier bridge.
    loci = parse_gbff(GBFF_FILE)
    old2rs: dict[str, list[str]] = {}
    for lt, d in loci.items():
        if d["old_locus_tag"]:
            old2rs.setdefault(d["old_locus_tag"], []).append(lt)
    kegg_list = {}
    for line in kegg.text(f"list/{KEGG_ORG}").splitlines():
        p = line.split("\t")
        if len(p) >= 4 and p[0].startswith(KEGG_ORG + ":"):
            kegg_list[p[0].split(":", 1)[1]] = {"type": p[1], "position": p[2]}
    bridge_counts = {"exact": 0, "overlap_same_strand": 0, "mismatch": 0, "not_compared": 0}
    no_refseq = []
    for g, rec in kegg_list.items():
        rs = old2rs.get(g, [])
        if len(rs) != 1:
            no_refseq.append(g)
            continue
        d = loci[rs[0]]
        bridge_counts[position_check(rec["position"], d["location"], d["replicon"])] += 1
    id_bridge = {
        "method": "KEGG mmi ids are original Sanger locus tags (MMAR_xxxx); DESeq2 ids are RefSeq "
                  "PGAP locus tags (MMAR_RSxxxxx). Joined through the RefSeq GBFF /old_locus_tag "
                  "qualifier; each pair checked by KEGG position (list/mmi) vs RefSeq gene coordinates.",
        "refseq_loci": len(loci),
        "refseq_loci_with_old_locus_tag": sum(1 for d in loci.values() if d["old_locus_tag"]),
        "refseq_old_locus_tags_shared_by_several_loci": sum(1 for v in old2rs.values() if len(v) > 1),
        "kegg_genes": len(kegg_list),
        "kegg_genes_with_one_refseq_locus": len(kegg_list) - len(no_refseq),
        "kegg_genes_without_refseq_locus": len(no_refseq),
        "refseq_old_locus_tags_absent_from_kegg": sum(1 for o in old2rs if o not in kegg_list),
        "position_check_counts": bridge_counts,
    }

    # Gene records for every linked gene and every old locus.
    old_loci = {}
    for a in OLD_ASSIGNMENTS:
        d = loci.get(a[4])
        old_loci[a[4]] = {"in_refseq": d is not None,
                          "old_locus_tag": d["old_locus_tag"] if d else "",
                          "refseq_gene": d["gene"] if d else "",
                          "refseq_product": d["product"] if d else ""}
    gene_ids = sorted({g for gs in links.values() for g in gs}
                      | {v["old_locus_tag"] for v in old_loci.values() if v["old_locus_tag"]})
    recs = kegg.get_entries(KEGG_ORG, gene_ids)
    genes = {}
    for g in gene_ids:
        info = kegg_gene_info(recs[g])
        rs = old2rs.get(g, [])
        d = loci[rs[0]] if len(rs) == 1 else None
        genes[g] = {
            "refseq_locus_tag": rs[0] if len(rs) == 1 else "",
            "refseq_gene": d["gene"] if d else "",
            "refseq_product": d["product"] if d else "",
            "refseq_location": d["location"] if d else "",
            "refseq_pseudo": d["pseudo"] if d else False,
            **info,
            "position_check": position_check(info["kegg_position"], d["location"], d["replicon"]) if d
            else "no_refseq_locus",
        }

    mtimes = sorted(os.path.getmtime(fn) for fn in kegg.used_files)
    day = lambda t: dt.datetime.fromtimestamp(t).date().isoformat()  # noqa: E731
    snap = {
        "description": "KEGG data used by analysis/07_cellular_metabolic_modeling.py. Raw responses are "
                       "cached in data/reference/kegg/ (git-ignored); this file is the tracked snapshot.",
        "source": {
            "kegg_rest": KEGG_BASE,
            "kegg_release": release or "not stated by info/kegg; see kegg_database_dates",
            "kegg_database_dates": db_dates,
            "kegg_organism": " | ".join(mmi_info),
            "retrieved_first": day(mtimes[0]) if mtimes else "",
            "retrieved_last": day(mtimes[-1]) if mtimes else "",
            "cache_files_used": len(kegg.used_files),
            "user_agent": USER_AGENT,
            "refseq_gbff": rel(GBFF_FILE),
            "refseq_gbff_sha256": sha256(GBFF_FILE),
        },
        "id_bridge": id_bridge,
        "queries": queries,
        "genes": genes,
        "old_loci": old_loci,
    }
    with open(SNAPSHOT, "w") as f:
        json.dump(snap, f, indent=1)
        f.write("\n")
    print(f"  KEGG requests sent: {kegg.requests}; snapshot written to {rel(SNAPSHOT)}")
    return snap


def snapshot_covers(snap: dict) -> bool:
    rq, oq = needed_queries()
    if any(q not in snap.get("queries", {}) for q in rq + oq):
        return False
    if any(a[4] not in snap.get("old_loci", {}) for a in OLD_ASSIGNMENTS):
        return False
    needed = {g for q in rq for g in snap["queries"][q]["genes"]}
    return needed <= set(snap.get("genes", {}))


# ================================================================ statistics
def load_deseq(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"gene_id": str})
    need = {"gene_id", "log2FoldChange", "log2FoldChange_shrunk", "pvalue", "padj", "significant"}
    missing = need - set(df.columns)
    if missing:
        sys.exit(f"Error: {rel(path)} lacks columns {sorted(missing)}")
    if df["gene_id"].duplicated().any():
        sys.exit(f"Error: duplicated gene_id values in {rel(path)}")
    df = df.set_index("gene_id")
    df["tested"] = df["padj"].notna()
    df["deg"] = df["tested"] & (df["significant"].astype(str).str.upper() == "YES")
    return df


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return []
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out.tolist()


def summarise(keys: list[str], de: pd.DataFrame) -> dict:
    """Per-contrast summary over the tested genes among `keys` (DESeq2 gene ids)."""
    present = [k for k in dict.fromkeys(keys) if k in de.index]
    sub = de.loc[present]
    sub = sub[sub["tested"]]
    s = {"n_in_table": len(present), "n_tested": int(len(sub)), "n_deg": int(sub["deg"].sum()),
         "n_up": int((sub["deg"] & (sub["log2FoldChange"] > 0)).sum()),
         "n_down": int((sub["deg"] & (sub["log2FoldChange"] < 0)).sum()),
         "median_lfc_shrunk": float(sub["log2FoldChange_shrunk"].median()) if len(sub) else None,
         "lead_gene": None, "lead_padj": None, "lead_lfc_shrunk": None}
    if len(sub):
        lead = sub.assign(_id=sub.index).sort_values(["padj", "pvalue", "_id"], kind="mergesort").iloc[0]
        s.update(lead_gene=lead["_id"], lead_padj=float(lead["padj"]),
                 lead_lfc_shrunk=float(lead["log2FoldChange_shrunk"]))
    return s


# ================================================================ formatting
def _missing(x) -> bool:
    return x is None or (isinstance(x, float) and math.isnan(x))


def f_lfc(x) -> str:
    return "NA" if _missing(x) else f"{x:.4f}"


def f_p(x) -> str:
    return "NA" if _missing(x) else f"{x:.3e}"


def f_num(x, nd=3) -> str:
    return "NA" if _missing(x) else f"{x:.{nd}f}"


def j_num(x, sig=6):
    return None if _missing(x) else float(f"{x:.{sig}g}")


def write_tsv(path: str, header: list[str], rows: list[list]) -> None:
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join("" if v is None else str(v) for v in r) + "\n")


def count(values) -> dict:
    out: dict = {}
    for v in values:
        out[v] = out.get(v, 0) + 1
    return out


# ====================================================================== main
def main() -> None:
    ap = argparse.ArgumentParser(description="KEGG-anchored reaction and subsystem summary (OSD-528)")
    ap.add_argument("--rebuild-kegg", action="store_true",
                    help="rebuild the KEGG snapshot from the cache, downloading what is missing")
    ap.add_argument("--refresh-kegg", action="store_true",
                    help="re-download every KEGG response, then rebuild the snapshot")
    args = ap.parse_args()

    print("=== Phase 7: KEGG-anchored reaction and subsystem summary (OSD-528) ===")

    # ---- 1. KEGG snapshot
    snap = None
    if not (args.rebuild_kegg or args.refresh_kegg) and os.path.exists(SNAPSHOT):
        with open(SNAPSHOT) as f:
            snap = json.load(f)
        if not snapshot_covers(snap):
            print("  Snapshot does not cover the current reaction list; rebuilding.")
            snap = None
    if snap is None:
        snap = build_snapshot(refresh=args.refresh_kegg)
    Q, GENES, OLD = snap["queries"], snap["genes"], snap["old_loci"]
    print(f"  KEGG database dates {snap['source'].get('kegg_database_dates')}; retrieved "
          f"{snap['source']['retrieved_first']}..{snap['source']['retrieved_last']}")

    # ---- 2. DESeq2 tables and identifier style
    des = {}
    for prefix, cid, fname, label in CONTRASTS:
        path = os.path.join(DATA_PROCESSED, fname)
        if not os.path.exists(path):
            sys.exit(f"Error: {rel(path)} not found")
        des[cid] = load_deseq(path)
    ids = set().union(*(set(d.index) for d in des.values()))
    n_rs = sum(1 for i in ids if re.match(r"^MMAR_RS\d+$", i))
    n_old = sum(1 for i in ids if re.match(r"^MMAR_(\d+|p\d+)$", i))
    id_style = "refseq_locus_tag" if n_rs >= n_old else "kegg_locus_tag"
    print(f"  DESeq2 gene ids: {n_rs} MMAR_RS-style, {n_old} MMAR_xxxx-style -> joining on {id_style}")

    def de_key(g: str) -> str:
        return GENES[g]["refseq_locus_tag"] if id_style == "refseq_locus_tag" else g

    def symbol(g: str) -> str:
        return GENES[g]["kegg_symbol"] or GENES[g]["refseq_gene"]

    def reaction_genes(r: dict) -> list[str]:
        return sorted({g for q in r["kegg_queries"] for g in Q[q]["genes"]})

    # ---- 3. Reactions: KEGG genes and per-contrast summaries
    kept, dropped = [], []
    for r in REACTIONS:
        genes = reaction_genes(r)
        if not genes:
            reason = r["note"] if not r["kegg_queries"] else \
                f"No mmi gene linked to {', '.join(r['kegg_queries'])}. {r['note']}"
            dropped.append({"reaction_id": r["reaction_id"], "reaction_name": r["reaction_name"],
                            "kegg_queries": r["kegg_queries"], "replaces": r["replaces"],
                            "reason": reason.strip()})
            continue
        ecs = [q[3:] for q in r["kegg_queries"] if q.startswith("ec:")]
        rx = dict(r)
        rx["genes"] = genes
        rx["ec_number"] = " / ".join(ecs)
        rx["mapping_basis"] = "EC" if len(ecs) == len(r["kegg_queries"]) else ("KO" if not ecs else "EC+KO")
        rx["kegg_name"] = " / ".join(Q[q]["kegg_name"] for q in r["kegg_queries"])
        eq = " | ".join(Q[q]["kegg_equation"] for q in r["kegg_queries"] if Q[q]["kegg_equation"])
        rx["equation"] = eq if eq else "n/a (mapped by KO; KEGG gives no reaction equation)"
        rx["stats"] = {cid: summarise([de_key(g) for g in genes if de_key(g)], des[cid])
                       for _, cid, _, _ in CONTRASTS}
        kept.append(rx)
    kept_ids = {r["reaction_id"] for r in kept}
    key2gene = {de_key(g): g for rx in kept for g in rx["genes"] if de_key(g)}

    # ---- 4. Audit of the superseded assignments
    new_by_old = {o: r for r in REACTIONS for o in r["replaces"]}
    audit = []
    for old_id, old_sub, old_name, old_ec, old_gid, old_sym in OLD_ASSIGNMENTS:
        r = new_by_old[old_id]
        new_genes = reaction_genes(r)
        o = OLD[old_gid]
        tag = o["old_locus_tag"]
        gi = GENES.get(tag, {}) if tag else {}
        syms = {s.lower() for s in re.split(r"[,\s]+", gi.get("kegg_symbol", "")) if s}
        if o["refseq_gene"]:
            syms.add(o["refseq_gene"].lower())
        sym_match = old_sym.lower() in syms
        in_set = bool(tag) and tag in new_genes
        if not o["in_refseq"]:
            verdict, correct = "locus_absent_from_refseq_annotation", "NO"
        elif in_set and sym_match:
            verdict, correct = "correct", "YES"
        elif in_set:
            verdict, correct = "wrong_subunit_of_same_reaction", "NO"
        elif new_genes:
            verdict, correct = "not_in_kegg_gene_set", "NO"
        elif r["kegg_queries"]:
            verdict, correct = "kegg_has_no_mmi_gene_for_query", "NO"
        else:
            verdict, correct = "unverifiable_no_kegg_reaction_or_ortholog", "UNVERIFIABLE"
        status = "dropped" if r["reaction_id"] not in kept_ids else ("kept" if correct == "YES" else "re-mapped")
        old_ec_names = []
        for e in old_ec_list(old_ec):
            qi = Q.get(f"ec:{e}", {})
            old_ec_names.append(f"{e}: {qi.get('kegg_name', '')} (mmi genes linked: {len(qi.get('genes', []))})")
        if not old_ec_list(old_ec):
            old_ec_names.append("no EC given" if old_ec in ("", "N/A") else f"{old_ec}: not a complete EC number")
        audit.append({
            "old_reaction_id": old_id, "old_subsystem": old_sub, "old_reaction_name": old_name,
            "old_ec": old_ec, "old_ec_in_kegg": "; ".join(old_ec_names),
            "old_gene_id": old_gid, "old_gene_symbol": old_sym,
            "old_locus_kegg_id": tag or "NA",
            "old_locus_refseq_gene": o["refseq_gene"],
            "old_locus_refseq_product": o["refseq_product"] if o["in_refseq"] else "not in RefSeq annotation",
            "old_locus_kegg_symbol": gi.get("kegg_symbol", ""), "old_locus_kegg_name": gi.get("kegg_name", ""),
            "old_locus_kegg_orthology": gi.get("kegg_orthology", ""),
            "old_locus_in_kegg_gene_set": "YES" if in_set else "NO",
            "old_symbol_matches_locus": "YES" if sym_match else "NO",
            "old_locus_correct": correct, "verdict": verdict, "status": status,
            "merged_into": r["reaction_id"] if len(r["replaces"]) > 1 else "",
            "new_reaction_id": r["reaction_id"] if r["reaction_id"] in kept_ids else "",
            "kegg_queries": ",".join(r["kegg_queries"]) or "none",
            "new_kegg_genes": ";".join(new_genes),
            "new_refseq_loci": ";".join(GENES[g]["refseq_locus_tag"] or "NA" for g in new_genes),
            "new_symbols": ";".join(symbol(g) or "-" for g in new_genes),
            "note": r["note"],
        })

    # ---- 5. Subsystems
    sub_genes = {s: sorted({g for rx in kept if rx["subsystem"] == s for g in rx["genes"]})
                 for s in SUBSYSTEM_ORDER}
    sub_rx = {s: [rx["reaction_id"] for rx in kept if rx["subsystem"] == s] for s in SUBSYSTEM_ORDER}
    subsystems = [s for s in SUBSYSTEM_ORDER if sub_rx[s]]
    sub_stats: dict[str, dict] = {s: {} for s in subsystems}
    backgrounds = {}
    for prefix, cid, fname, label in CONTRASTS:
        de = des[cid]
        N, K = int(de["tested"].sum()), int(de["deg"].sum())
        backgrounds[cid] = {"n_tested": N, "n_deg": K}
        pvals = []
        for s in subsystems:
            st = summarise([de_key(g) for g in sub_genes[s] if de_key(g)], de)
            n, k = st["n_tested"], st["n_deg"]
            st["n_genes_mapped"] = len(sub_genes[s])
            st["deg_fraction"] = k / n if n else None
            st["fold_enrichment"] = (k / n) / (K / N) if n and K else None
            st["hypergeom_p"] = float(hypergeom.sf(k - 1, N, K, n)) if n else None
            pvals.append(1.0 if st["hypergeom_p"] is None else st["hypergeom_p"])
            sub_stats[s][cid] = st
        for s, q in zip(subsystems, bh(pvals)):
            sub_stats[s][cid]["bh_fdr"] = None if sub_stats[s][cid]["hypergeom_p"] is None else q

    sub_rows = []
    for prefix, cid, fname, label in CONTRASTS:
        order = sorted(subsystems, key=lambda s: (2.0 if sub_stats[s][cid]["hypergeom_p"] is None
                                                  else sub_stats[s][cid]["hypergeom_p"], s))
        for s in order:
            st = sub_stats[s][cid]
            sub_rows.append([
                s, cid, len(sub_rx[s]), st["n_genes_mapped"], st["n_tested"], st["n_deg"], st["n_up"],
                st["n_down"], f_num(st["deg_fraction"]), backgrounds[cid]["n_tested"], backgrounds[cid]["n_deg"],
                f_num(st["fold_enrichment"]), f_p(st["hypergeom_p"]), f_p(st["bh_fdr"]),
                f_lfc(st["median_lfc_shrunk"]),
                st["lead_gene"] or "NA", (symbol(key2gene[st["lead_gene"]]) or "-") if st["lead_gene"] else "NA",
                f_lfc(st["lead_lfc_shrunk"]), f_p(st["lead_padj"]),
            ])

    # ---- 6. Write outputs
    # 6a. reactions
    rx_header = ["reaction_id", "subsystem", "compartment", "reaction_name", "mapping_basis", "kegg_queries",
                 "ec_number", "kegg_name", "equation", "kegg_scope", "n_kegg_genes", "kegg_gene_ids",
                 "gene_ids", "gene_symbols", "old_assignment", "old_locus_correct"]
    for prefix, *_ in CONTRASTS:
        rx_header += [f"{prefix}_n_tested", f"{prefix}_n_deg", f"{prefix}_n_up", f"{prefix}_n_down",
                      f"{prefix}_median_log2fc_shrunk", f"{prefix}_lead_gene_id", f"{prefix}_lead_gene_symbol",
                      f"{prefix}_lead_log2fc_shrunk", f"{prefix}_lead_padj"]
    audit_by_old = {a["old_reaction_id"]: a for a in audit}
    rx_rows = []
    for rx in kept:
        olds = [audit_by_old[o] for o in rx["replaces"]]
        row = [rx["reaction_id"], rx["subsystem"], rx["compartment"], rx["reaction_name"], rx["mapping_basis"],
               ",".join(rx["kegg_queries"]), rx["ec_number"] or "NA", rx["kegg_name"], rx["equation"], rx["scope"],
               len(rx["genes"]), ";".join(rx["genes"]),
               ";".join(GENES[g]["refseq_locus_tag"] or "NA" for g in rx["genes"]),
               ";".join(symbol(g) or "-" for g in rx["genes"]),
               "; ".join(f"{a['old_reaction_id']}:{a['old_gene_id']}({a['old_gene_symbol']})" for a in olds),
               "; ".join(f"{a['old_reaction_id']}:{a['old_locus_correct']}" for a in olds) if len(olds) > 1
               else olds[0]["old_locus_correct"]]
        for prefix, cid, *_ in CONTRASTS:
            st = rx["stats"][cid]
            row += [st["n_tested"], st["n_deg"], st["n_up"], st["n_down"], f_lfc(st["median_lfc_shrunk"]),
                    st["lead_gene"] or "NA",
                    (symbol(key2gene[st["lead_gene"]]) or "-") if st["lead_gene"] else "NA",
                    f_lfc(st["lead_lfc_shrunk"]), f_p(st["lead_padj"])]
        rx_rows.append(row)
    write_tsv(OUT_REACTIONS, rx_header, rx_rows)

    # 6b. subsystems
    write_tsv(OUT_SUBSYSTEMS, [
        "subsystem", "contrast", "n_reactions", "n_genes_mapped", "n_genes_tested", "n_deg", "n_up", "n_down",
        "deg_fraction", "background_n_tested", "background_n_deg", "fold_enrichment", "hypergeom_p", "bh_fdr",
        "median_log2fc_shrunk", "lead_gene_id", "lead_gene_symbol", "lead_log2fc_shrunk", "lead_padj"],
        sub_rows)

    # 6c. gene-level table
    g_header = ["reaction_id", "subsystem", "kegg_gene_id", "refseq_locus_tag", "kegg_symbol", "refseq_gene",
                "kegg_name", "kegg_orthology", "position_check"]
    for prefix, *_ in CONTRASTS:
        g_header += [f"{prefix}_log2fc", f"{prefix}_log2fc_shrunk", f"{prefix}_padj", f"{prefix}_significant"]
    g_rows = []
    for rx in kept:
        for g in rx["genes"]:
            gi = GENES[g]
            row = [rx["reaction_id"], rx["subsystem"], g, gi["refseq_locus_tag"] or "NA", gi["kegg_symbol"],
                   gi["refseq_gene"], gi["kegg_name"], gi["kegg_orthology"], gi["position_check"]]
            for prefix, cid, *_ in CONTRASTS:
                de, k = des[cid], de_key(g)
                if k and k in de.index:
                    d = de.loc[k]
                    row += [f_lfc(float(d["log2FoldChange"])), f_lfc(float(d["log2FoldChange_shrunk"])),
                            f_p(float(d["padj"])) if d["tested"] else "NA", "YES" if d["deg"] else "NO"]
                else:
                    row += ["NA", "NA", "NA", "not_in_table"]
            g_rows.append(row)
    write_tsv(OUT_GENES, g_header, g_rows)

    # 6d. audit
    a_header = list(audit[0].keys())
    write_tsv(OUT_AUDIT, a_header, [[a[k] for k in a_header] for a in audit])

    # 6e. JSON model
    model_genes = sorted({g for rx in kept for g in rx["genes"]})
    in_tables = {cid: sum(1 for g in model_genes if de_key(g) and de_key(g) in des[cid].index)
                 for _, cid, *_ in CONTRASTS}
    pos_counts = count(GENES[g]["position_check"] for g in model_genes)
    status_counts = count(a["status"] for a in audit)
    correct_counts = count(a["old_locus_correct"] for a in audit)

    def stats_json(st: dict, sub: bool = False) -> dict:
        out = {"n_tested": st["n_tested"], "n_deg": st["n_deg"], "n_up": st["n_up"], "n_down": st["n_down"],
               "median_log2fc_shrunk": j_num(st["median_lfc_shrunk"]),
               "lead_gene_id": st["lead_gene"], "lead_log2fc_shrunk": j_num(st["lead_lfc_shrunk"]),
               "lead_padj": j_num(st["lead_padj"])}
        if sub:
            out.update(n_genes_mapped=st["n_genes_mapped"], deg_fraction=j_num(st["deg_fraction"]),
                       fold_enrichment=j_num(st["fold_enrichment"]), hypergeom_p=j_num(st["hypergeom_p"]),
                       bh_fdr=j_num(st["bh_fdr"]))
        return out

    model = {
        "model_id": "M_marinum_OSD528_KEGG_reaction_sets",
        "organism": "Mycobacterium marinum M (RefSeq GCF_000018345.1; KEGG organism mmi, T00694)",
        "study_accession": "NASA OSDR OSD-528",
        "description": (f"{len(kept)} curated reactions and named proteins in {len(subsystems)} subsystems. "
                        "Each is linked to M. marinum M genes only through KEGG EC or KO links and summarised "
                        "with DESeq2 statistics. Not a genome-scale model; no flux balance analysis or flux "
                        "predictions; no gene-protein-reaction logic is asserted. The file name is historical; "
                        "the content is not SBML."),
        "generated_by": "analysis/07_cellular_metabolic_modeling.py",
        "kegg": snap["source"],
        "gene_identifier_bridge": {
            **snap["id_bridge"],
            "deseq2_id_style": id_style,
            "model_genes": len(model_genes),
            "model_genes_position_check": pos_counts,
            "model_genes_in_deseq2_table": in_tables,
        },
        "contrasts": [{"id": cid, "label": label, "deseq2_table": f"data/processed/{fname}",
                       "column_prefix": prefix, "background": backgrounds[cid]}
                      for prefix, cid, fname, label in CONTRASTS],
        "definitions": {
            "reaction_genes": "Union of mmi genes returned by KEGG link/mmi/<query> for the reaction's EC "
                              "(preferred) or KO queries. No gene is added by hand.",
            "tested_gene": "Gene present in the contrast's DESeq2 table with non-missing padj.",
            "deg": "significant == 'YES' in the DESeq2 table (BH-adjusted Wald p < 0.05).",
            "up_down": "DEG with log2FoldChange > 0 (up) or < 0 (down).",
            "reaction_lead_gene": "Tested gene with the smallest padj (ties: smaller pvalue, then gene id); "
                                  "reported with its padj and apeglm-shrunken log2FC.",
            "median_log2fc_shrunk": "Median of log2FoldChange_shrunk (apeglm) over the tested genes.",
            "subsystem_genes": "De-duplicated union of the subsystem's reaction genes.",
            "deg_fraction": "n_deg / n_tested.",
            "fold_enrichment": "(n_deg / n_tested) / (background DEGs / background tested genes).",
            "hypergeom_p": "One-sided hypergeometric P(X >= n_deg), X ~ Hypergeom(N = tested genes in the "
                           "contrast, K = DEGs in the contrast, n = subsystem tested genes).",
            "bh_fdr": "Benjamini-Hochberg FDR over the subsystems within each contrast.",
            "removed_metrics": "The earlier Subsystem Perturbation Index, the mean-of-simulators log2FC with "
                               "+/-0.5 direction calls, and the 'predominant flux shift' labels are removed.",
        },
        "subsystem_renames": [{"old": o, "new": n, "reason": why} for o, (n, why) in SUBSYSTEM_RENAMES.items()],
        "old_assignment_summary": {"assignments": len(audit), "status": status_counts,
                                   "old_locus_correct": correct_counts},
        "merged_reactions": [{"new_reaction_id": r["reaction_id"], "old_reaction_ids": r["replaces"],
                              "reason": r["note"]} for r in kept if len(r["replaces"]) > 1],
        "dropped_reactions": dropped,
        "subsystems": [{"subsystem": s, "reactions": sub_rx[s], "kegg_genes": sub_genes[s],
                        "statistics": {cid: stats_json(sub_stats[s][cid], sub=True) for _, cid, *_ in CONTRASTS}}
                       for s in subsystems],
        "reactions": [{
            "reaction_id": rx["reaction_id"], "subsystem": rx["subsystem"], "compartment": rx["compartment"],
            "reaction_name": rx["reaction_name"], "mapping_basis": rx["mapping_basis"],
            "kegg_queries": rx["kegg_queries"], "ec_number": rx["ec_number"] or None,
            "kegg_name": rx["kegg_name"], "equation": rx["equation"], "kegg_scope": rx["scope"],
            "gene_reaction_rule": None,
            "genes": [{"kegg_gene_id": g, "refseq_locus_tag": GENES[g]["refseq_locus_tag"] or None,
                       "symbol": symbol(g) or None, "kegg_name": GENES[g]["kegg_name"],
                       "kegg_orthology": GENES[g]["kegg_orthology"] or None,
                       "position_check": GENES[g]["position_check"]} for g in rx["genes"]],
            "replaces": rx["replaces"], "note": rx["note"] or None,
            "statistics": {cid: stats_json(rx["stats"][cid]) for _, cid, *_ in CONTRASTS},
        } for rx in kept],
        "old_assignment_audit": audit,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(model, f, indent=2)
        f.write("\n")

    # ---- 7. Console summary
    print(f"  Reactions: {len(kept)} retained, {len(dropped)} dropped "
          f"({', '.join(d['reaction_id'] for d in dropped)}); {len(model_genes)} KEGG genes")
    print(f"  Old assignments: {status_counts}; old locus correct: {correct_counts}")
    print(f"  Model genes found in DESeq2 tables: {in_tables}; position check: {pos_counts}")
    for prefix, cid, fname, label in CONTRASTS:
        print(f"  {label} (background {backgrounds[cid]['n_deg']}/{backgrounds[cid]['n_tested']} DEGs):")
        for row in [r for r in sub_rows if r[1] == cid]:
            print(f"    {row[0]:<50s} DEG {row[5]:>3}/{row[4]:<3} up {row[6]:>3} down {row[7]:>3} "
                  f"median {row[14]:>8}  p {row[12]}  FDR {row[13]}")
    for p in (OUT_REACTIONS, OUT_SUBSYSTEMS, OUT_JSON, OUT_GENES, OUT_AUDIT):
        print(f"  wrote {rel(p)}")


if __name__ == "__main__":
    main()
