"""
GenoStruct Interpreter — Genomic Variant Dataset Generator
Generates 100,000 synthetic genomic variant records modeled after ClinVar/gnomAD/COSMIC
Reference datasets: ClinVar (https://www.ncbi.nlm.nih.gov/clinvar/), gnomAD v3.1, COSMIC v97
"""

import json
import random
import csv
import os
from datetime import datetime, timedelta

random.seed(42)

# ── Reference data modeled after real databases ─────────────────────────────

CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]

GENES_ONCOLOGY = [
    "TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "PIK3CA", "APC", "PTEN",
    "RB1", "VHL", "MLH1", "MSH2", "MSH6", "PMS2", "CDH1", "STK11",
    "PALB2", "ATM", "CHEK2", "NBN", "RAD51C", "RAD51D", "BARD1", "BRIP1",
    "MUTYH", "NF1", "NF2", "RET", "MEN1", "SDHA", "SDHB", "SDHC", "SDHD",
    "TSC1", "TSC2", "WT1", "DICER1", "FLCN", "FH", "BAP1", "SMAD4",
    "BMPR1A", "POLD1", "POLE", "AXIN2", "GREM1", "EPCAM", "CDK4",
    "CDKN2A", "MITF", "TMEM127", "MAX", "PRKAR1A", "AIP", "CASR",
    "CTNNA1", "MAP3K1", "RECQL", "BRAF", "ALK", "MET", "ROS1", "FGFR1",
    "FGFR2", "FGFR3", "IDH1", "IDH2", "DNMT3A", "FLT3", "NPM1", "CEBPA",
    "RUNX1", "ASXL1", "TET2", "JAK2", "MPL", "CALR", "CSF3R", "SETBP1",
    "EZH2", "KMT2A", "KMT2D", "CREBBP", "EP300", "ARID1A", "ARID1B",
    "SMARCA4", "SMARCB1", "SMARCC1", "KDM6A", "KDM5C", "BCOR", "BCORL1",
    "KDM2B", "PHF6", "NOTCH1", "NOTCH2", "FBXW7", "PTPN11", "CBL",
    "NRAS", "HRAS", "NF1", "RASA1", "RASA2", "SPRY4", "RAF1", "MAP2K1",
]

VARIANT_TYPES = ["SNV", "InDel", "CNV", "Fusion", "Frameshift", "Splice", "Nonsense", "Missense", "Silent"]
VARIANT_WEIGHTS = [0.45, 0.20, 0.08, 0.04, 0.08, 0.05, 0.04, 0.04, 0.02]

PATHOGENICITY = ["Benign", "Likely Benign", "Variant of Uncertain Significance", "Likely Pathogenic", "Pathogenic"]
PATH_WEIGHTS   = [0.28, 0.22, 0.20, 0.14, 0.16]

AMINO_ACIDS = ["Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile",
               "Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val"]
AA_CODES    = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

PROTEIN_DOMAINS = [
    "DNA-binding domain","Kinase domain","SH2 domain","SH3 domain","PH domain",
    "BRCT domain","RING finger domain","WD40 repeat","Ankyrin repeat","Leucine zipper",
    "Helix-loop-helix","Zinc finger domain","RAS-binding domain","Death domain",
    "Tudor domain","Chromodomain","Bromodomain","PHD finger","CARD domain","SAM domain",
    "Transmembrane domain","Extracellular domain","Ligand-binding domain","Coiled-coil",
    "OB fold","VHL box","BRCA2 OB fold","RAD51-binding domain","NLS","NES",
]

STRUCTURAL_EFFECTS = [
    "Disrupts alpha-helix formation",
    "Destabilizes beta-sheet",
    "Blocks ATP-binding pocket",
    "Prevents dimerization",
    "Alters substrate recognition loop",
    "Disrupts disulfide bond",
    "Impairs nuclear localization",
    "Reduces thermal stability",
    "Alters electrostatic surface",
    "Steric clash with co-factor",
    "Loss of hydrophobic core packing",
    "Disrupts pi-pi stacking",
    "Impairs zinc coordination",
    "Alters allosteric regulation",
    "Disrupts protein-protein interface",
    "No significant structural change",
    "Minor backbone perturbation",
    "Conservative substitution",
]

DRUGS_BY_TARGET = {
    "EGFR":   ["Erlotinib","Gefitinib","Afatinib","Osimertinib","Dacomitinib"],
    "BRCA1":  ["Olaparib","Niraparib","Rucaparib","Talazoparib"],
    "BRCA2":  ["Olaparib","Niraparib","Rucaparib","Talazoparib"],
    "BRAF":   ["Vemurafenib","Dabrafenib","Encorafenib"],
    "ALK":    ["Crizotinib","Alectinib","Brigatinib","Lorlatinib"],
    "ROS1":   ["Crizotinib","Entrectinib","Lorlatinib"],
    "KRAS":   ["Sotorasib","Adagrasib"],
    "PIK3CA": ["Alpelisib","Copanlisib","Idelalisib"],
    "JAK2":   ["Ruxolitinib","Fedratinib","Pacritinib"],
    "FLT3":   ["Midostaurin","Gilteritinib","Quizartinib"],
    "IDH1":   ["Ivosidenib","Olutasidenib"],
    "IDH2":   ["Enasidenib"],
    "MET":    ["Capmatinib","Tepotinib","Crizotinib"],
    "RET":    ["Selpercatinib","Pralsetinib"],
    "FGFR2":  ["Pemigatinib","Infigratinib","Futibatinib"],
    "FGFR3":  ["Erdafitinib"],
    "NTRK1":  ["Larotrectinib","Entrectinib"],
}
DEFAULT_DRUGS = ["Pembrolizumab","Nivolumab","Atezolizumab","Durvalumab","Ipilimumab",
                 "Trastuzumab","Bevacizumab","Cetuximab","Panitumumab","Ramucirumab"]

CANCER_TYPES = [
    "Lung Adenocarcinoma","Colorectal Cancer","Breast Cancer","Glioblastoma",
    "Pancreatic Ductal Adenocarcinoma","Ovarian Cancer","Prostate Cancer",
    "Melanoma","Acute Myeloid Leukemia","Chronic Lymphocytic Leukemia",
    "Non-Hodgkin Lymphoma","Hepatocellular Carcinoma","Cervical Cancer",
    "Bladder Cancer","Renal Cell Carcinoma","Thyroid Cancer","Gastric Cancer",
    "Endometrial Cancer","Head and Neck SCC","Multiple Myeloma",
    "Diffuse Large B-cell Lymphoma","Follicular Lymphoma","Mantle Cell Lymphoma",
    "Chronic Myeloid Leukemia","Essential Thrombocythemia","Polycythemia Vera",
]

EVIDENCE_LEVELS = ["Level 1A", "Level 1B", "Level 2A", "Level 2B", "Level 3", "Level 4"]
EVIDENCE_WEIGHTS = [0.10, 0.15, 0.20, 0.20, 0.20, 0.15]

ZYGOSITY = ["Heterozygous", "Homozygous", "Hemizygous", "Compound Heterozygous"]

FUNCTIONAL_IMPACT = ["Loss of Function", "Gain of Function", "Dominant Negative",
                     "Altered Binding Affinity", "Unknown", "No Impact"]

# ── Helper functions ─────────────────────────────────────────────────────────

def rand_position():
    return random.randint(1_000, 250_000_000)

def rand_hgvs_c():
    pos = random.randint(1, 10000)
    from_nt = random.choice("ACGT")
    to_nt   = random.choice([n for n in "ACGT" if n != from_nt])
    return f"c.{pos}{from_nt}>{to_nt}"

def rand_hgvs_p(vtype):
    if vtype in ("SNV", "Missense"):
        aa1 = random.choice(AMINO_ACIDS)
        pos = random.randint(1, 1500)
        aa2 = random.choice([a for a in AMINO_ACIDS if a != aa1])
        return f"p.{aa1}{pos}{aa2}"
    elif vtype == "Nonsense":
        aa1 = random.choice(AMINO_ACIDS)
        pos = random.randint(1, 1500)
        return f"p.{aa1}{pos}*"
    elif vtype == "Frameshift":
        aa1 = random.choice(AMINO_ACIDS)
        pos = random.randint(1, 1500)
        return f"p.{aa1}{pos}fs*{random.randint(1,50)}"
    elif vtype == "Silent":
        aa1 = random.choice(AMINO_ACIDS)
        pos = random.randint(1, 1500)
        return f"p.{aa1}{pos}="
    else:
        return "p.?"

def rand_allele_freq():
    r = random.random()
    if r < 0.5:
        return round(random.uniform(0, 0.001), 6)
    elif r < 0.8:
        return round(random.uniform(0.001, 0.01), 5)
    else:
        return round(random.uniform(0.01, 0.5), 4)

def rand_cadd_score(pathogenicity):
    base = {"Benign": (0,10), "Likely Benign": (5,15),
            "Variant of Uncertain Significance": (10,25),
            "Likely Pathogenic": (20,30), "Pathogenic": (25,40)}
    lo, hi = base[pathogenicity]
    return round(random.uniform(lo, hi), 2)

def rand_revel_score(pathogenicity):
    base = {"Benign": (0,0.3), "Likely Benign": (0.1,0.4),
            "Variant of Uncertain Significance": (0.3,0.7),
            "Likely Pathogenic": (0.5,0.8), "Pathogenic": (0.7,1.0)}
    lo, hi = base[pathogenicity]
    return round(random.uniform(lo, hi), 4)

def rand_drug(gene):
    pool = DRUGS_BY_TARGET.get(gene, DEFAULT_DRUGS)
    return random.choice(pool)

def rand_date():
    start = datetime(2010, 1, 1)
    delta = datetime.now() - start
    return (start + timedelta(days=random.randint(0, delta.days))).strftime("%Y-%m-%d")

def rand_clinvar_id():
    return f"VCV{random.randint(100000,999999):09d}"

def rand_cosmic_id():
    return f"COSM{random.randint(10000,9999999)}"

def rand_dbsnp_id():
    return f"rs{random.randint(1000000,2000000000)}"

def rand_tumor_content():
    return round(random.uniform(0.05, 0.98), 2)

def rand_read_depth():
    return random.randint(50, 5000)

def rand_vaf():
    return round(random.uniform(0.02, 0.95), 3)

def rand_sift(pathogenicity):
    if pathogenicity in ("Pathogenic","Likely Pathogenic"):
        return random.choice(["Deleterious","Deleterious (low confidence)"])
    return random.choice(["Tolerated","Tolerated (low confidence)","Deleterious"])

def rand_polyphen(pathogenicity):
    if pathogenicity in ("Pathogenic","Likely Pathogenic"):
        return random.choice(["Probably Damaging","Possibly Damaging"])
    return random.choice(["Benign","Possibly Damaging","Probably Damaging"])

# ── Main generation ──────────────────────────────────────────────────────────

def generate_row(idx):
    gene       = random.choice(GENES_ONCOLOGY)
    chrom      = random.choice(CHROMOSOMES[:22]) if random.random() > 0.05 else random.choice(CHROMOSOMES)
    vtype      = random.choices(VARIANT_TYPES, weights=VARIANT_WEIGHTS)[0]
    pathog     = random.choices(PATHOGENICITY, weights=PATH_WEIGHTS)[0]
    domain     = random.choice(PROTEIN_DOMAINS)
    struct_eff = random.choice(STRUCTURAL_EFFECTS)
    cancer     = random.choice(CANCER_TYPES)
    evid       = random.choices(EVIDENCE_LEVELS, weights=EVIDENCE_WEIGHTS)[0]
    drug       = rand_drug(gene)
    pos        = rand_position()
    hgvs_c     = rand_hgvs_c()
    hgvs_p     = rand_hgvs_p(vtype)
    cadd       = rand_cadd_score(pathog)
    revel      = rand_revel_score(pathog)
    af_gnomad  = rand_allele_freq()
    af_cosmic  = rand_allele_freq() if vtype in ("SNV","Missense","Nonsense") else 0.0
    vaf        = rand_vaf()
    depth      = rand_read_depth()
    tumor_pct  = rand_tumor_content()

    # Confidence score (0-100)
    base_conf = {"Benign": 60, "Likely Benign": 55, "Variant of Uncertain Significance": 40,
                 "Likely Pathogenic": 70, "Pathogenic": 85}
    confidence = min(100, max(0, base_conf[pathog] + random.randint(-15, 15)))

    # Structural score (higher = more disrupted)
    struct_score = round(random.uniform(0.1, 0.9)
                         if "No significant" not in struct_eff and "Conservative" not in struct_eff and "Minor" not in struct_eff
                         else random.uniform(0.01, 0.25), 3)

    ref_allele = random.choice("ACGT")
    alt_allele = random.choice([a for a in "ACGT" if a != ref_allele])

    return {
        "variant_id":          f"GS{idx:07d}",
        "gene":                gene,
        "chromosome":          chrom,
        "position":            pos,
        "ref_allele":          ref_allele,
        "alt_allele":          alt_allele,
        "hgvs_c":              hgvs_c,
        "hgvs_p":              hgvs_p,
        "variant_type":        vtype,
        "pathogenicity":       pathog,
        "pathogenicity_score": round(
            {"Benign":0.05,"Likely Benign":0.2,"Variant of Uncertain Significance":0.5,
             "Likely Pathogenic":0.78,"Pathogenic":0.94}[pathog] + random.uniform(-0.08,0.08), 3),
        "cadd_phred":          cadd,
        "revel_score":         revel,
        "sift_prediction":     rand_sift(pathog),
        "polyphen2_prediction":rand_polyphen(pathog),
        "gnomad_af":           af_gnomad,
        "cosmic_af":           af_cosmic,
        "vaf":                 vaf,
        "read_depth":          depth,
        "tumor_content":       tumor_pct,
        "zygosity":            random.choice(ZYGOSITY),
        "protein_domain":      domain,
        "structural_effect":   struct_eff,
        "structural_score":    struct_score,
        "functional_impact":   random.choice(FUNCTIONAL_IMPACT),
        "cancer_type":         cancer,
        "recommended_drug":    drug,
        "drug_class":          random.choice(["Targeted Therapy","Immunotherapy","Chemotherapy","PARP Inhibitor","Hormone Therapy"]),
        "evidence_level":      evid,
        "confidence_score":    confidence,
        "clinvar_id":          rand_clinvar_id(),
        "cosmic_id":           rand_cosmic_id() if vtype in ("SNV","Missense","Nonsense","Frameshift") else "",
        "dbsnp_id":            rand_dbsnp_id() if random.random() > 0.3 else "",
        "transcript_id":       f"NM_{random.randint(100000,999999)}.{random.randint(1,5)}",
        "exon_number":         f"Exon {random.randint(1,30)}",
        "inheritance_pattern": random.choice(["Autosomal Dominant","Autosomal Recessive","X-linked","de novo","Somatic","Germline"]),
        "reported_date":       rand_date(),
        "review_status":       random.choice(["Reviewed by expert panel","Multiple submitters","Single submitter","No assertion criteria"]),
        "lab_source":          random.choice(["MSKCC","Mayo Clinic","Broad Institute","UCSF","Johns Hopkins","MD Anderson","Dana-Farber","Stanford","NCI","Ambry Genetics"]),
        "population":          random.choice(["European","African","East Asian","South Asian","Latino","Ashkenazi Jewish","Finnish","Other"]),
        "sex":                 random.choice(["Male","Female","Unknown"]),
        "age_at_diagnosis":    random.randint(18, 90),
        "msi_status":          random.choice(["MSI-H","MSI-L","MSS","Unknown"]),
        "tmb":                 round(random.uniform(0.5, 60.0), 1),
        "copy_number":         random.choice(["Normal","Gain","Loss","Amplification","Deletion"]) if vtype == "CNV" else "Normal",
        "allele_origin":       random.choice(["Somatic","Germline","Unknown"]),
        "hotspot":             random.choice([True, False, False, False]),
        "actionable":          pathog in ("Pathogenic","Likely Pathogenic") and random.random() > 0.2,
        "tier":                random.choice(["Tier I","Tier II","Tier III","Tier IV"]),
    }

# ── Write outputs ─────────────────────────────────────────────────────────────

print("Generating 100,000 genomic variant records...")
N = 100_000
rows = [generate_row(i+1) for i in range(N)]

# CSV
csv_path = os.path.join(os.path.dirname(__file__), "genomic_variants_100k.csv")
with open(csv_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
print(f"  ✓ CSV  → {csv_path}  ({N:,} rows)")

# Sample JSON (first 500 rows for the live app)
sample = rows[:500]
json_path = os.path.join(os.path.dirname(__file__), "sample_variants.json")
with open(json_path, "w", encoding="utf-8") as f:
    json.dump(sample, f, indent=2)
print(f"  ✓ JSON sample → {json_path}  (500 rows)")

# Summary stats
from collections import Counter
path_counts = Counter(r["pathogenicity"] for r in rows)
gene_counts  = Counter(r["gene"] for r in rows)
print("\n── Pathogenicity distribution ──")
for k, v in sorted(path_counts.items(), key=lambda x: -x[1]):
    print(f"  {k:<42} {v:>7,}  ({v/N*100:.1f}%)")
print("\n── Top 10 genes ──")
for g, c in gene_counts.most_common(10):
    print(f"  {g:<12} {c:>7,}")
print("\nDone.")
