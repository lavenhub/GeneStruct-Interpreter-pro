# GenoStruct Interpreter — Precision Genomics Platform

A complete, deployable clinical genomics decision-support system with AI-powered variant analysis, 3D protein visualization, and precision medicine reporting.

---

## 🗂 Project Structure

```
genocstruct/
├── index.html                    ← Main application (single-file, deployable)
├── data/
│   ├── generate_dataset.py       ← Dataset generator script
│   ├── genomic_variants_100k.csv ← Pre-generated dataset (100,000 rows, 44 MB)
│   └── sample_variants.json      ← Sample data for the live app (500 rows)
└── README.md                     ← This file
```

---

## 🚀 Quick Start

### Option 1: Open directly in browser (no server needed)
```bash
open index.html
# or double-click index.html in your file explorer
```

### Option 2: Local server (recommended for file uploads)
```bash
# Python 3
python3 -m http.server 8080
# then open http://localhost:8080

# Node.js
npx serve .
```

### Option 3: Deploy to any static host
Upload `index.html` to:
- **Netlify**: drag & drop into netlify.com/drop
- **Vercel**: `npx vercel deploy`
- **GitHub Pages**: push to `gh-pages` branch
- **AWS S3**: static website hosting

---

## 🧬 Dataset

### Pre-generated Dataset
- **File**: `data/genomic_variants_100k.csv`
- **Rows**: 100,000
- **Size**: ~44 MB
- **Reference Databases**: ClinVar, gnomAD v3.1, COSMIC v97

### Dataset Schema (44 columns)

| Column | Type | Description |
|--------|------|-------------|
| variant_id | string | Unique identifier (GS0000001…) |
| gene | string | HGNC gene symbol (107 cancer genes) |
| chromosome | string | Chromosome (1–22, X, Y, MT) |
| position | int | Genomic position (GRCh38) |
| ref_allele | string | Reference nucleotide |
| alt_allele | string | Alternate nucleotide |
| hgvs_c | string | HGVS coding DNA notation |
| hgvs_p | string | HGVS protein notation |
| variant_type | enum | SNV, InDel, CNV, Frameshift, Splice, Nonsense, Missense, Silent |
| pathogenicity | enum | Benign / Likely Benign / VUS / Likely Pathogenic / Pathogenic |
| pathogenicity_score | float | AI confidence score (0–1) |
| cadd_phred | float | CADD PHRED score (0–40+) |
| revel_score | float | REVEL pathogenicity score (0–1) |
| sift_prediction | enum | Deleterious / Tolerated |
| polyphen2_prediction | enum | Probably Damaging / Possibly Damaging / Benign |
| gnomad_af | float | gnomAD allele frequency |
| cosmic_af | float | COSMIC allele frequency |
| vaf | float | Variant allele frequency (tumor) |
| read_depth | int | Sequencing read depth |
| tumor_content | float | Tumor purity (0–1) |
| zygosity | enum | Heterozygous / Homozygous / Hemizygous |
| protein_domain | string | Affected protein domain |
| structural_effect | string | Predicted 3D structural consequence |
| structural_score | float | Structural disruption score (0–1) |
| functional_impact | enum | Loss of Function / Gain of Function / etc. |
| cancer_type | string | Primary cancer type (26 types) |
| recommended_drug | string | FDA-approved targeted therapy |
| drug_class | string | Drug class category |
| evidence_level | string | Level 1A – Level 4 (CIViC-style) |
| confidence_score | int | AI confidence (0–100) |
| clinvar_id | string | ClinVar accession |
| cosmic_id | string | COSMIC mutation ID |
| dbsnp_id | string | dbSNP rs identifier |
| transcript_id | string | RefSeq transcript |
| exon_number | string | Affected exon |
| inheritance_pattern | enum | Somatic / Germline / Autosomal Dominant / etc. |
| reported_date | date | Date first reported |
| review_status | string | ClinVar review status |
| lab_source | string | Reporting laboratory |
| population | string | Population ancestry |
| sex | string | Biological sex |
| age_at_diagnosis | int | Age at diagnosis |
| msi_status | enum | MSI-H / MSI-L / MSS |
| tmb | float | Tumor mutational burden (mutations/Mb) |
| copy_number | enum | Normal / Gain / Loss / Amplification / Deletion |
| allele_origin | enum | Somatic / Germline / Unknown |
| hotspot | bool | Recurrent hotspot mutation |
| actionable | bool | Clinically actionable variant |
| tier | enum | Tier I–IV (AMP/ASCO guidelines) |

### Regenerating the Dataset
```bash
cd data/
python3 generate_dataset.py
```
This will regenerate:
- `genomic_variants_100k.csv` (100,000 rows)
- `sample_variants.json` (first 500 rows for the app)

---

## 🖥️ Application Features

### Core Diagnostic Engine

#### 1. Dashboard
- Real-time stats: 100K variants, 15,750 pathogenic, 20,253 VUS, 47 drug targets
- Pathogenicity donut chart
- Top 10 genes by variant count (bar chart)
- Variant type breakdown (polar area chart)
- Recently analyzed pathogenic variants

#### 2. Universal Genomic Uploader
- Drag-and-drop for `.vcf`, `.fastq`, `.json`, `.gz`
- Real-time upload progress simulation with step-by-step status
- HIPAA-compliant local processing disclosure
- Max file size: 5 GB

#### 3. Variant Explorer
- Paginated table of all 500 sample variants
- Search by gene, HGVS, variant ID, cancer type
- Filter by pathogenicity, variant type, gene
- Sortable columns
- Click-through to 3D Viewer

#### 4. 3D Protein Visualizer (WebGL / Three.js)
- Interactive protein backbone with:
  - Blue alpha-helices
  - Green beta-sheets
  - Orange loop regions
  - Red mutation site (pulsing animation)
  - Magenta drug-binding pocket (torus)
  - Side-chain atoms
  - Dark starfield background
- Mouse drag to rotate, scroll to zoom
- Touch support for mobile
- Controls: Zoom ±, Wireframe toggle, Auto-spin
- Structural analysis panel with binding affinity
- Drug binding analysis with alternative therapy

#### 5. Pathogenicity Scoring Dashboard
- Visual 5-tier scale with live percentages
- CADD score histogram
- REVEL score by pathogenicity class
- Top pathogenic hotspot variants table

### Clinical Support

#### 6. XAI Heatmap (Explainable AI)
- 200-residue attention heatmap
- Color-coded by AI attention score (white → deep red)
- Hover tooltips showing residue, position, attention %
- Top 5 contributing residues ranked by score
- Radar chart: 5 evidence dimensions

#### 7. Precision Medicine Report
- Structured clinical report layout
- Auto-populated from selected variant
- Variant identification, AI assessment, structural impact
- Drug recommendations with evidence levels
- Population genetics (gnomAD, ClinVar, COSMIC)
- PDF download via browser print

#### 8. Scientific Justification Feed (Digital Medical Representative)
- 10 curated clinical literature entries
- Filter by: All / Clinical / Drug / Structure / Review
- Full-text search
- Drug reimbursement status panel

#### 9. Comparative Variant Analysis
- Side-by-side wild-type vs mutant 3D viewers
- Structural difference bar chart
- Key metrics: stability, binding, DNA contact, thermal stability, RMSF

---

## 🎨 Design System

| Property | Value |
|----------|-------|
| Background | #FAFAF8 (warm white primer) |
| Font Display | DM Serif Display |
| Font Body | Plus Jakarta Sans |
| Font Mono | DM Mono |
| Accent | #1A1A6E / #3949AB |
| Pathogenic Red | #B71C1C |
| Benign Green | #2E7D32 |
| VUS Orange | #F57C00 |

---

## 🔬 Genes Covered (107 cancer genes)

TP53, BRCA1, BRCA2, EGFR, KRAS, PIK3CA, APC, PTEN, RB1, VHL, MLH1, MSH2, MSH6, PMS2, CDH1, STK11, PALB2, ATM, CHEK2, NBN, RAD51C, RAD51D, BARD1, BRIP1, MUTYH, NF1, NF2, RET, MEN1, SDHA, SDHB, SDHC, SDHD, TSC1, TSC2, WT1, DICER1, FLCN, FH, BAP1, SMAD4, BMPR1A, POLD1, POLE, AXIN2, GREM1, EPCAM, CDK4, CDKN2A, MITF, TMEM127, MAX, PRKAR1A, AIP, CASR, CTNNA1, MAP3K1, RECQL, BRAF, ALK, MET, ROS1, FGFR1, FGFR2, FGFR3, IDH1, IDH2, DNMT3A, FLT3, NPM1, CEBPA, RUNX1, ASXL1, TET2, JAK2, MPL, CALR, CSF3R, SETBP1, EZH2, KMT2A, KMT2D, CREBBP, EP300, ARID1A, ARID1B, SMARCA4, SMARCB1, SMARCC1, KDM6A, KDM5C, BCOR, BCORL1, NOTCH1, NOTCH2, FBXW7, PTPN11, CBL, NRAS, HRAS, RASA1, RAF1, MAP2K1

---

## 📋 Clinical Standards

- Pathogenicity classification follows **ACMG/AMP 2015 guidelines**
- Evidence levels follow **CIViC** and **OncoKB** framework (1A–4)
- Variant tiers follow **AMP/ASCO/CAP 2017** guidelines (Tier I–IV)
- Drug recommendations sourced from **FDA-approved indications**

---

## ⚠️ Disclaimer

This platform is a **clinical decision support tool** intended for use by qualified healthcare professionals. All AI-generated findings must be validated by a board-certified clinical geneticist or oncologist. This software does not constitute a final medical diagnosis.

---

*GenoStruct Interpreter v1.0 — Built with Three.js, Chart.js, and vanilla web technologies. No backend required.*
