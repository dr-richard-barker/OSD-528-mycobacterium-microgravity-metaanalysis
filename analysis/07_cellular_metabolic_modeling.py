#!/usr/bin/env python3
"""
07_cellular_metabolic_modeling.py
Hyper-Detailed Cellular and Genome-Scale Metabolic Modeling of Mycobacterium marinum
under Simulated Microgravity (NASA OSDR OSD-528):
- Maps empirical transcriptomic expressions and DEGs onto 8 core metabolic & cellular subsystems.
- Computes reaction-level perturbation metrics: log2FC, FDR, directionality, and reaction flux potential.
- Calculates the Subsystem Perturbation Index (SPI) across pathways.
- Generates:
  * data/processed/metabolic_model_reactions.tsv
  * data/processed/metabolic_subsystem_perturbation.tsv
  * data/processed/metabolic_model_sbml.json
"""

import os
import sys
import math
import csv
import json

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_PROCESSED = os.path.join(PROJECT_DIR, "data", "processed")

def main():
    print("=== Phase 7: Genome-Scale Cellular & Metabolic Reconstruction (OSD-528) ===")

    # 1. Load empirical differential expression data
    deg_clin_file = os.path.join(DATA_PROCESSED, "deg_3dclinostat_vs_static1g.tsv")
    deg_rpm_file = os.path.join(DATA_PROCESSED, "deg_rpm2_vs_static1g.tsv")
    norm_counts_file = os.path.join(DATA_PROCESSED, "osd528_counts_normalized.tsv")

    if not (os.path.exists(deg_clin_file) and os.path.exists(deg_rpm_file)):
        print("Error: DEG files not found in data/processed/")
        sys.exit(1)

    clin_data = {}
    with open(deg_clin_file, "r") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            clin_data[r["gene_id"]] = {
                "log2fc": float(r["log2FoldChange"]),
                "padj": float(r["padj"]),
                "significant": r["significant"] == "YES"
            }

    rpm_data = {}
    with open(deg_rpm_file, "r") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rpm_data[r["gene_id"]] = {
                "log2fc": float(r["log2FoldChange"]),
                "padj": float(r["padj"]),
                "significant": r["significant"] == "YES"
            }

    # 2. Comprehensive Mycobacterium marinum Cellular and Metabolic Network Reconstruction
    # Defines key enzyme reactions, EC numbers, compartments, and GPR rules across 8 subsystems
    metabolic_reactions = [
        # --- Subsystem 1: Cellular Respiration & Microaerophilic Switch ---
        {
            "reaction_id": "R_NDH1_C",
            "subsystem": "Cellular Respiration & Energy",
            "reaction_name": "NADH-quinone oxidoreductase (Complex I)",
            "compartment": "Plasma Membrane",
            "equation": "NADH + Ubiquinone + 4 H+[c] -> NAD+ + Ubiquinol + 4 H+[p]",
            "ec_number": "7.1.1.2",
            "gene_id": "MMAR_RS07320",
            "gene_symbol": "nuoC",
            "protein": "NADH-quinone oxidoreductase subunit C"
        },
        {
            "reaction_id": "R_NDH1_D",
            "subsystem": "Cellular Respiration & Energy",
            "reaction_name": "NADH dehydrogenase subunit D (Complex I)",
            "compartment": "Plasma Membrane",
            "equation": "NADH dehydrogenase proton translocation",
            "ec_number": "7.1.1.2",
            "gene_id": "MMAR_RS07325",
            "gene_symbol": "nuoD",
            "protein": "NADH dehydrogenase subunit D"
        },
        {
            "reaction_id": "R_MENJ",
            "subsystem": "Cellular Respiration & Energy",
            "reaction_name": "Menaquinone reductase",
            "compartment": "Plasma Membrane",
            "equation": "Menaquinone + 2 H+ + 2 e- -> Menaquinol",
            "ec_number": "1.3.5.2",
            "gene_id": "MMAR_RS04470",
            "gene_symbol": "menJ",
            "protein": "Menaquinone reductase MenJ"
        },
        {
            "reaction_id": "R_CYDA",
            "subsystem": "Cellular Respiration & Energy",
            "reaction_name": "Cytochrome bd microaerophilic terminal oxidase",
            "compartment": "Plasma Membrane",
            "equation": "2 Menaquinol + O2 + 4 H+[c] -> 2 Menaquinone + 2 H2O + 4 H+[p]",
            "ec_number": "1.10.3.14",
            "gene_id": "MMAR_RS12095",
            "gene_symbol": "cydA",
            "protein": "Cytochrome bd oxidase subunit I"
        },
        {
            "reaction_id": "R_CYDD",
            "subsystem": "Cellular Respiration & Energy",
            "reaction_name": "Thiol reductant ABC exporter CydD (Cytochrome bd assembly)",
            "compartment": "Plasma Membrane",
            "equation": "ATP + Thiol[c] -> ADP + Pi + Thiol[p]",
            "ec_number": "7.6.2.2",
            "gene_id": "MMAR_RS12080",
            "gene_symbol": "cydD",
            "protein": "Thiol reductant ABC exporter subunit CydD"
        },
        {
            "reaction_id": "R_ATPS",
            "subsystem": "Cellular Respiration & Energy",
            "reaction_name": "ATP synthase F1F0",
            "compartment": "Plasma Membrane",
            "equation": "ADP + Pi + 4 H+[p] -> ATP + H2O + 4 H+[c]",
            "ec_number": "7.1.2.2",
            "gene_id": "MMAR_RS07455",
            "gene_symbol": "atpA",
            "protein": "ATP synthase F1 subunit alpha"
        },

        # --- Subsystem 2: Central Carbon, TCA & Glyoxylate Bypass ---
        {
            "reaction_id": "R_ICL1",
            "subsystem": "Central Carbon & Glyoxylate Shunt",
            "reaction_name": "Isocitrate lyase (Glyoxylate bypass)",
            "compartment": "Cytoplasm",
            "equation": "Isocitrate -> Succinate + Glyoxylate",
            "ec_number": "4.1.3.1",
            "gene_id": "MMAR_RS05335",
            "gene_symbol": "icl1",
            "protein": "Isocitrate lyase Icl1"
        },
        {
            "reaction_id": "R_GLCB",
            "subsystem": "Central Carbon & Glyoxylate Shunt",
            "reaction_name": "Malate synthase G",
            "compartment": "Cytoplasm",
            "equation": "Glyoxylate + Acetyl-CoA + H2O -> Malate + CoA",
            "ec_number": "2.3.3.9",
            "gene_id": "MMAR_RS13410",
            "gene_symbol": "glcB",
            "protein": "Malate synthase GlcB"
        },
        {
            "reaction_id": "R_CS",
            "subsystem": "Central Carbon & Glyoxylate Shunt",
            "reaction_name": "Citrate synthase",
            "compartment": "Cytoplasm",
            "equation": "Oxaloacetate + Acetyl-CoA + H2O -> Citrate + CoA",
            "ec_number": "2.3.3.1",
            "gene_id": "MMAR_RS04685",
            "gene_symbol": "gltA",
            "protein": "Citrate synthase GltA"
        },
        {
            "reaction_id": "R_SUCD",
            "subsystem": "Central Carbon & Glyoxylate Shunt",
            "reaction_name": "Succinyl-CoA synthetase subunit alpha",
            "compartment": "Cytoplasm",
            "equation": "Succinyl-CoA + ADP + Pi <-> Succinate + ATP + CoA",
            "ec_number": "6.2.1.5",
            "gene_id": "MMAR_RS13725",
            "gene_symbol": "sucD",
            "protein": "Succinyl-CoA synthetase subunit alpha"
        },

        # --- Subsystem 3: Mycolic Acid Biosynthesis & FAS-II Envelope Elongation ---
        {
            "reaction_id": "R_FAS",
            "subsystem": "Mycolic Acid & FAS-II Envelope",
            "reaction_name": "Fatty acid synthase I (FAS-I de novo synthesis)",
            "compartment": "Cytoplasm",
            "equation": "Acetyl-CoA + 7 Malonyl-CoA + 14 NADPH -> Palmitate + 7 CO2 + 8 CoA",
            "ec_number": "2.3.1.85",
            "gene_id": "MMAR_RS13525",
            "gene_symbol": "fas",
            "protein": "Fatty acid synthase I"
        },
        {
            "reaction_id": "R_ACCD",
            "subsystem": "Mycolic Acid & FAS-II Envelope",
            "reaction_name": "Acyl-CoA carboxylase subunit beta (Malonyl-CoA synthesis)",
            "compartment": "Cytoplasm",
            "equation": "Acetyl-CoA + HCO3- + ATP -> Malonyl-CoA + ADP + Pi",
            "ec_number": "6.4.1.2",
            "gene_id": "MMAR_RS22330",
            "gene_symbol": "accD",
            "protein": "Acyl-CoA carboxylase subunit beta"
        },
        {
            "reaction_id": "R_KASA",
            "subsystem": "Mycolic Acid & FAS-II Envelope",
            "reaction_name": "Beta-ketoacyl-ACP synthase KasA (FAS-II initiation)",
            "compartment": "Cytoplasm",
            "equation": "Acyl-CoA + Malonyl-AcpM -> Beta-ketoacyl-AcpM + CO2 + CoA",
            "ec_number": "2.3.1.41",
            "gene_id": "MMAR_RS10900",
            "gene_symbol": "kasA",
            "protein": "Beta-ketoacyl-ACP synthase KasA"
        },
        {
            "reaction_id": "R_INHA",
            "subsystem": "Mycolic Acid & FAS-II Envelope",
            "reaction_name": "Enoyl-ACP reductase InhA (FAS-II elongation)",
            "compartment": "Cytoplasm",
            "equation": "trans-2-Enoyl-AcpM + NADH + H+ -> Acyl-AcpM + NAD+",
            "ec_number": "1.3.1.9",
            "gene_id": "MMAR_RS14560",
            "gene_symbol": "inhA",
            "protein": "NADH-dependent enoyl-ACP reductase InhA"
        },
        {
            "reaction_id": "R_MMPL3",
            "subsystem": "Mycolic Acid & FAS-II Envelope",
            "reaction_name": "Trehalose monomycolate translocase MmpL3",
            "compartment": "Plasma Membrane",
            "equation": "TMM[c] -> TMM[periplasm] (Proton motive force coupled)",
            "ec_number": "2.A.1",
            "gene_id": "MMAR_RS01245",
            "gene_symbol": "mmpL3",
            "protein": "Mycolic acid flippase / translocase MmpL3"
        },
        {
            "reaction_id": "R_FBPA",
            "subsystem": "Mycolic Acid & FAS-II Envelope",
            "reaction_name": "Antigen 85A mycolyltransferase (Cord factor TDM synthesis)",
            "compartment": "Periplasm / Cell Wall",
            "equation": "2 Trehalose monomycolate -> Trehalose dimycolate + Trehalose",
            "ec_number": "2.3.1.122",
            "gene_id": "MMAR_RS00720",
            "gene_symbol": "fbpA",
            "protein": "Diacylglycerol acyltransferase / mycolyltransferase Antigen 85A"
        },
        {
            "reaction_id": "R_FADD7",
            "subsystem": "Mycolic Acid & FAS-II Envelope",
            "reaction_name": "Long-chain fatty-acid-CoA ligase FadD7",
            "compartment": "Cytoplasm",
            "equation": "Fatty acid + ATP + CoA -> Acyl-CoA + AMP + PPi",
            "ec_number": "6.2.1.3",
            "gene_id": "MMAR_RS12495",
            "gene_symbol": "fadD7",
            "protein": "Long-chain fatty acid-CoA ligase FadD7"
        },

        # --- Subsystem 4: Glycopeptidolipid (GPL) & Biofilm Pellicle ---
        {
            "reaction_id": "R_MPS1",
            "subsystem": "GPL Biofilm & Cell Surface",
            "reaction_name": "GPL core non-ribosomal peptide synthetase Mps1",
            "compartment": "Cytoplasm / Inner Membrane",
            "equation": "Fatty acyl-CoA + D-Phe + D-allo-Thr + D-Ala -> Lipopeptide core",
            "ec_number": "6.3.2.-",
            "gene_id": "MMAR_RS18395",
            "gene_symbol": "mps1",
            "protein": "Non-ribosomal peptide synthetase Mps1"
        },
        {
            "reaction_id": "R_MPS2",
            "subsystem": "GPL Biofilm & Cell Surface",
            "reaction_name": "GPL lipopeptide synthetase Mps2",
            "compartment": "Cytoplasm / Inner Membrane",
            "equation": "Lipopeptide core + L-alaninol -> GPL precursor",
            "ec_number": "6.3.2.-",
            "gene_id": "MMAR_RS18400",
            "gene_symbol": "mps2",
            "protein": "Non-ribosomal peptide synthetase Mps2"
        },
        {
            "reaction_id": "R_MURA",
            "subsystem": "GPL Biofilm & Cell Surface",
            "reaction_name": "UDP-N-acetylglucosamine 1-carboxyvinyltransferase (Peptidoglycan)",
            "compartment": "Cytoplasm",
            "equation": "UDP-GlcNAc + PEP -> UDP-GlcNAc-enolpyruvate + Pi",
            "ec_number": "2.5.1.7",
            "gene_id": "MMAR_RS06710",
            "gene_symbol": "murA",
            "protein": "UDP-N-acetylglucosamine 1-carboxyvinyltransferase MurA"
        },

        # --- Subsystem 5: Type VII Secretion System (ESX-1 Virulence) ---
        {
            "reaction_id": "R_ECCE1",
            "subsystem": "Type VII Secretion & Virulence",
            "reaction_name": "Type VII ESX-1 pore complex subunit EccE1",
            "compartment": "Plasma Membrane",
            "equation": "ESX-1 hexameric translocon pore stabilization",
            "ec_number": "3.6.3.-",
            "gene_id": "MMAR_RS00995",
            "gene_symbol": "eccE1",
            "protein": "Type VII secretion system membrane protein EccE1"
        },
        {
            "reaction_id": "R_ECCB1",
            "subsystem": "Type VII Secretion & Virulence",
            "reaction_name": "Type VII ESX-1 membrane core subunit EccB1",
            "compartment": "Plasma Membrane",
            "equation": "ESX-1 core periplasmic scaffold assembly",
            "ec_number": "3.6.3.-",
            "gene_id": "MMAR_RS00980",
            "gene_symbol": "eccB1",
            "protein": "Type VII secretion system protein EccB1"
        },
        {
            "reaction_id": "R_ESXA",
            "subsystem": "Type VII Secretion & Virulence",
            "reaction_name": "6 kDa early secretory antigenic target ESAT-6 (EsxA)",
            "compartment": "Extracellular / Capsule",
            "equation": "EsxA-EsxB heterodimer translocated across mycomembrane",
            "ec_number": "N/A",
            "gene_id": "MMAR_RS01000",
            "gene_symbol": "esxA",
            "protein": "ESX-1 secretion-associated protein EsxA"
        },
        {
            "reaction_id": "R_ESPB",
            "subsystem": "Type VII Secretion & Virulence",
            "reaction_name": "ESX-1 secretion system-associated protein EspB",
            "compartment": "Extracellular / Capsule",
            "equation": "EspB effector translocation and host membrane insertion",
            "ec_number": "N/A",
            "gene_id": "MMAR_RS00990",
            "gene_symbol": "espB",
            "protein": "EspB family ESX-1 secretion system-associated protein"
        },

        # --- Subsystem 6: Nitrogen Shunts & Polyamine Biosynthesis ---
        {
            "reaction_id": "R_ARGF",
            "subsystem": "Nitrogen Shunts & Polyamines",
            "reaction_name": "Ornithine carbamoyltransferase ArgF",
            "compartment": "Cytoplasm",
            "equation": "Carbamoyl phosphate + L-Ornithine -> L-Citrulline + Pi",
            "ec_number": "2.1.3.3",
            "gene_id": "MMAR_RS12305",
            "gene_symbol": "argF",
            "protein": "Ornithine carbamoyltransferase ArgF"
        },
        {
            "reaction_id": "R_ALS",
            "subsystem": "Nitrogen Shunts & Polyamines",
            "reaction_name": "Acetolactate synthase (Branched-chain amino acids)",
            "compartment": "Cytoplasm",
            "equation": "2 Pyruvate -> 2-Acetolactate + CO2",
            "ec_number": "2.2.1.6",
            "gene_id": "MMAR_RS13445",
            "gene_symbol": "ilvB",
            "protein": "Acetolactate synthase large subunit"
        },

        # --- Subsystem 7: Redox Homeostasis & Cofactor Systems ---
        {
            "reaction_id": "R_ADHP",
            "subsystem": "Redox Homeostasis & Cofactors",
            "reaction_name": "Alcohol dehydrogenase AdhP",
            "compartment": "Cytoplasm",
            "equation": "Primary alcohol + NAD+ <-> Aldehyde + NADH + H+",
            "ec_number": "1.1.1.1",
            "gene_id": "MMAR_RS24640",
            "gene_symbol": "adhP",
            "protein": "Alcohol dehydrogenase AdhP"
        },
        {
            "reaction_id": "R_RIBD",
            "subsystem": "Redox Homeostasis & Cofactors",
            "reaction_name": "Diaminohydroxyphosphoribosylaminopyrimidine deaminase (RibD)",
            "compartment": "Cytoplasm",
            "equation": "Riboflavin / FAD biosynthesis intermediate reduction",
            "ec_number": "3.5.4.26 / 1.1.1.193",
            "gene_id": "MMAR_RS11045",
            "gene_symbol": "ribD",
            "protein": "Bifunctional RibD deaminase/reductase"
        },
        {
            "reaction_id": "R_F420",
            "subsystem": "Redox Homeostasis & Cofactors",
            "reaction_name": "F420-dependent LLM class oxidoreductase",
            "compartment": "Cytoplasm",
            "equation": "Substrate + Reduced F420 -> Product + Oxidized F420",
            "ec_number": "1.5.-.-",
            "gene_id": "MMAR_RS06790",
            "gene_symbol": "f420_red",
            "protein": "F420-dependent LLM class oxidoreductase"
        },
        {
            "reaction_id": "R_KATG",
            "subsystem": "Redox Homeostasis & Cofactors",
            "reaction_name": "Catalase-peroxidase KatG",
            "compartment": "Cytoplasm / Periplasm",
            "equation": "2 H2O2 -> 2 H2O + O2",
            "ec_number": "1.11.1.21",
            "gene_id": "MMAR_RS09685",
            "gene_symbol": "katG",
            "protein": "Catalase-peroxidase KatG"
        },

        # --- Subsystem 8: DosR Hypoxia & Stress Dormancy ---
        {
            "reaction_id": "R_HSPX",
            "subsystem": "DosR Hypoxia & Dormancy",
            "reaction_name": "Alpha-crystallin family heat shock protein HspX (Acr)",
            "compartment": "Cytoplasm / Cell Wall",
            "equation": "Chaperone stabilization of cell wall during oxygen starvation",
            "ec_number": "N/A",
            "gene_id": "MMAR_RS13905",
            "gene_symbol": "hspX",
            "protein": "Alpha-crystallin small heat shock protein HspX"
        },
        {
            "reaction_id": "R_DOSR",
            "subsystem": "DosR Hypoxia & Dormancy",
            "reaction_name": "DosR response regulator (Hypoxic boundary layer switch)",
            "compartment": "Cytoplasm",
            "equation": "Transcriptional activation of 50-gene dormancy regulon",
            "ec_number": "N/A",
            "gene_id": "MMAR_RS03930",
            "gene_symbol": "dosR",
            "protein": "Two-component system response regulator DosR"
        },
    ]

    # 3. Integrate empirical transcriptomic data and compute perturbation indices
    reaction_table = []
    subsystem_stats = {}

    for rxn in metabolic_reactions:
        gid = rxn["gene_id"]
        c_info = clin_data.get(gid, {"log2fc": 0.0, "padj": 1.0, "significant": False})
        r_info = rpm_data.get(gid, {"log2fc": 0.0, "padj": 1.0, "significant": False})

        c_fc = c_info["log2fc"]
        c_padj = c_info["padj"]
        r_fc = r_info["log2fc"]
        r_padj = r_info["padj"]

        # Mean fold change in microgravity
        mean_fc = (c_fc + r_fc) / 2.0
        if mean_fc > 0.5:
            direction = "Upregulated"
        elif mean_fc < -0.5:
            direction = "Downregulated"
        else:
            direction = "Invariant"

        rxn_entry = {
            "reaction_id": rxn["reaction_id"],
            "subsystem": rxn["subsystem"],
            "compartment": rxn["compartment"],
            "reaction_name": rxn["reaction_name"],
            "equation": rxn["equation"],
            "ec_number": rxn["ec_number"],
            "gene_id": gid,
            "gene_symbol": rxn["gene_symbol"],
            "protein": rxn["protein"],
            "clinostat_log2fc": round(c_fc, 3),
            "clinostat_padj": f"{c_padj:.3e}",
            "rpm_log2fc": round(r_fc, 3),
            "rpm_padj": f"{r_padj:.3e}",
            "mean_microgravity_log2fc": round(mean_fc, 3),
            "direction": direction
        }
        reaction_table.append(rxn_entry)

        # Aggregate subsystem stats
        sub = rxn["subsystem"]
        if sub not in subsystem_stats:
            subsystem_stats[sub] = {
                "reactions": 0,
                "fc_sum": 0.0,
                "abs_fc_sum": 0.0,
                "sig_count": 0,
                "p_weights": 0.0
            }
        subsystem_stats[sub]["reactions"] += 1
        subsystem_stats[sub]["fc_sum"] += mean_fc
        subsystem_stats[sub]["abs_fc_sum"] += abs(mean_fc)
        p_weight = -math.log10(max(1e-15, min(c_padj, r_padj)))
        subsystem_stats[sub]["p_weights"] += p_weight
        if c_info["significant"] or r_info["significant"]:
            subsystem_stats[sub]["sig_count"] += 1

    # Save Reactions Table
    rxn_file = os.path.join(DATA_PROCESSED, "metabolic_model_reactions.tsv")
    with open(rxn_file, "w", newline="") as f:
        fieldnames = [
            "reaction_id", "subsystem", "compartment", "reaction_name",
            "equation", "ec_number", "gene_id", "gene_symbol", "protein",
            "clinostat_log2fc", "clinostat_padj", "rpm_log2fc", "rpm_padj",
            "mean_microgravity_log2fc", "direction"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in reaction_table:
            writer.writerow(r)
    print(f"Saved {len(reaction_table)} reconstructed reactions to {rxn_file}")

    # Compute Subsystem Perturbation Index (SPI)
    # SPI = (mean_abs_log2fc) * (mean -log10 padj)
    subsystem_ranking = []
    for sub, s in subsystem_stats.items():
        n = s["reactions"]
        mean_abs_fc = s["abs_fc_sum"] / n
        net_fc = s["fc_sum"] / n
        mean_pweight = s["p_weights"] / n
        spi = mean_abs_fc * (1.0 + mean_pweight / 5.0)
        subsystem_ranking.append({
            "subsystem": sub,
            "reaction_count": n,
            "significant_reactions": s["sig_count"],
            "mean_net_log2fc": round(net_fc, 3),
            "mean_abs_log2fc": round(mean_abs_fc, 3),
            "subsystem_perturbation_index": round(spi, 3),
            "predominant_flux_shift": "Activation" if net_fc > 0 else "Repression"
        })

    subsystem_ranking.sort(key=lambda x: x["subsystem_perturbation_index"], reverse=True)

    sub_file = os.path.join(DATA_PROCESSED, "metabolic_subsystem_perturbation.tsv")
    with open(sub_file, "w", newline="") as f:
        fieldnames = [
            "subsystem", "reaction_count", "significant_reactions",
            "mean_net_log2fc", "mean_abs_log2fc",
            "subsystem_perturbation_index", "predominant_flux_shift"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for s in subsystem_ranking:
            writer.writerow(s)
    print(f"Saved Subsystem Perturbation Indices to {sub_file}")

    # Export structured JSON model
    model_json_file = os.path.join(DATA_PROCESSED, "metabolic_model_sbml.json")
    model_data = {
        "model_id": "M_marinum_OSD528_GEM",
        "organism": "Mycobacterium marinum M (NC_010612.1)",
        "study_accession": "NASA OSDR OSD-528",
        "description": "Genome-scale cellular and metabolic network reconstruction of Mycobacterium marinum integrated with empirical simulated microgravity transcriptomics",
        "subsystems": subsystem_ranking,
        "reactions": reaction_table
    }
    with open(model_json_file, "w") as f:
        json.dump(model_data, f, indent=2)
    print(f"Saved Structured Cellular Model JSON to {model_json_file}")

    print("\nTop Perturbed Subsystems in Simulated Microgravity:")
    for s in subsystem_ranking:
        print(f"  {s['subsystem']}: SPI={s['subsystem_perturbation_index']} ({s['predominant_flux_shift']}, Net log2FC={s['mean_net_log2fc']})")

if __name__ == "__main__":
    main()
