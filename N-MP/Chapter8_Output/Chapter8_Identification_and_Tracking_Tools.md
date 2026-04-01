# Chapter 8: Identification and Tracking Tools

## Part III: Detection and Tracking in Biological Systems

---

**Book:** *Endogenous Microplastics and Nanoplastics (iMNPs): Toxicological Mechanisms and Health Impacts*

---

## Table of Contents

- 8.1 Introduction
- 8.2 Molecular Fingerprinting for iMNP Detection
- 8.3 Reference Database Development
- 8.4 Biodistribution Studies
- 8.5 Biological Barrier Crossing
- 8.6 Translocation Mechanisms
- 8.7 Retention and Clearance
- 8.8 Kinetic Modeling
- 8.9 Human Evidence: A Systematic Review and Meta-Analysis
- 8.10 Comparative Analysis
- 8.11 Conclusions and Future Perspectives
- References

---

## 8.1 Introduction

The detection, identification, and tracking of endogenous microplastics and nanoplastics (iMNPs) within biological systems represent one of the most technically challenging frontiers in environmental health science. Unlike monitoring plastics in environmental matrices such as water or sediment, detecting iMNPs in biological tissues requires navigating a complex landscape of organic matrices, low particle concentrations, and the ever-present risk of procedural contamination [1,2]. Over the past decade, rapid advances in analytical chemistry have expanded the toolkit available for iMNP identification, from established spectroscopic techniques such as Fourier-transform infrared (FTIR) spectroscopy and Raman spectroscopy to emerging mass spectrometric and imaging modalities [3,4].

This chapter provides a comprehensive review of the identification and tracking tools currently available for iMNP research in biological systems. We begin with the molecular fingerprinting approaches that underpin polymer identification (Section 8.2), followed by the development of reference databases essential for standardized analysis (Section 8.3). The chapter then transitions to the biological dimension: how iMNPs distribute across organ systems (Section 8.4), the mechanisms by which they cross biological barriers (Sections 8.5–8.6), and the kinetics governing their retention and clearance (Sections 8.7–8.8). A central feature of this chapter is a systematic meta-analysis of human evidence for iMNP detection across tissues (Section 8.9), synthesizing data from 11 studies encompassing 474 participants to provide pooled prevalence estimates. Finally, we offer a comparative analysis across species, particle types, and exposure routes (Section 8.10).

The fundamental objective of this chapter is to bridge the gap between analytical capability and biological understanding—demonstrating not only *how* we detect iMNPs but *where* they accumulate, *how* they get there, and *what* the quantitative evidence tells us about the prevalence of human exposure.

---

## 8.2 Molecular Fingerprinting for iMNP Detection

### 8.2.1 Spectral Libraries and Spectroscopic Identification

The cornerstone of iMNP identification in biological tissues is molecular fingerprinting through vibrational spectroscopy. Each polymer produces a characteristic spectral signature—a molecular fingerprint—that enables unambiguous identification of particle composition [5,6].

**Fourier-Transform Infrared (FTIR) Spectroscopy.** Micro-FTIR (μ-FTIR) operates by measuring the absorption of infrared radiation at specific wavelengths corresponding to molecular bond vibrations. For polymer identification, key diagnostic absorption bands include: C–H stretching (2850–2960 cm⁻¹) common to polyolefins such as polyethylene (PE) and polypropylene (PP); carbonyl C=O stretching (~1720 cm⁻¹) characteristic of polyethylene terephthalate (PET); and C–Cl stretching (~600–700 cm⁻¹) diagnostic of polyvinyl chloride (PVC) [7]. Attenuated total reflectance FTIR (ATR-FTIR) allows direct surface analysis without sample preparation, making it suitable for larger particles (>50 μm). Focal plane array (FPA)-based μ-FTIR imaging enables automated mapping of entire filter surfaces, achieving spatial resolution down to ~3 μm [8]. However, the diffraction limit of infrared radiation (~10 μm for mid-IR) constrains the minimum detectable particle size, fundamentally excluding nanoplastics from FTIR analysis.

**Raman Spectroscopy.** Micro-Raman spectroscopy (μ-Raman) complements FTIR by measuring inelastic light scattering from molecular vibrations. Raman is particularly advantageous for iMNP analysis because: (a) water does not interfere with Raman signals (unlike FTIR), facilitating analysis of hydrated biological samples; (b) spatial resolution can reach ~1 μm with visible laser excitation, and below 500 nm with confocal configurations; and (c) Raman provides complementary spectral information that can resolve ambiguities in FTIR identification [9,10]. However, fluorescence from biological matrix components can overwhelm Raman signals, necessitating careful sample preparation and photobleaching protocols. Surface-enhanced Raman spectroscopy (SERS) offers signal amplification of 10⁶–10⁸ fold via metallic nanostructures, potentially extending detection to individual nanoplastic particles [11].

**Pyrolysis-Gas Chromatography/Mass Spectrometry (Py-GC/MS).** Py-GC/MS provides mass-based quantification of polymer content by thermally decomposing plastic particles and analyzing the resulting volatile pyrolysis products [12]. This technique is uniquely suited for detecting nanoplastics below the spatial resolution limits of vibrational spectroscopy. Thermoextraction and desorption coupled with GC/MS (TED-GC/MS) offers a higher-throughput variant for larger sample volumes [13]. The principal limitation is that Py-GC/MS is a destructive technique that provides no morphological information (size, shape, color), and results are expressed as mass concentration (μg/g) rather than particle counts.

**Emerging Techniques.** Several advanced approaches are expanding the analytical frontier: time-of-flight secondary ion mass spectrometry (TOF-SIMS) provides surface-specific chemical mapping at nanometer resolution [14]; laser direct infrared (LDIR) chemical imaging offers rapid automated scanning [15]; and atomic force microscopy coupled with infrared spectroscopy (AFM-IR) achieves nanoscale chemical identification at spatial resolutions below 100 nm [16].

### 8.2.2 Polymer Additives as Tracers

Beyond identifying the polymer backbone, chemical additives incorporated during plastic manufacturing serve as auxiliary tracers for iMNP identification and source attribution [17]. Key additive classes include:

- **Plasticizers:** Phthalates (DEHP, DBP) and their non-phthalate replacements (DINCH, ATBC) are incorporated at 20–40% by weight in flexible PVC and can leach from particles within tissues, serving as indirect biomarkers of plastic exposure [18].
- **Flame retardants:** Brominated flame retardants (BDE-47, BDE-209) and organophosphate esters (TCEP, TDCIPP) associated with electronic waste plastics provide source-specific signatures [19].
- **UV stabilizers:** Benzotriazole UV absorbers (UV-326, UV-328) and hindered amine light stabilizers (HALS) indicate exposure to outdoor-weathered plastics [20].
- **Antioxidants:** Irganox 1076 and Irgafos 168 are widely used in polyolefin packaging and can be detected in tissues as co-contaminants [21].

The simultaneous detection of polymer backbone and additive profiles enables "chemical profiling" of iMNPs, providing information not only about polymer type but also about the original product category (packaging, textile, automotive component) and environmental history of the particle [22].

### 8.2.3 Degradation Markers

As plastic particles weather and degrade within biological environments, characteristic chemical modifications accumulate that serve as markers of aging and environmental history [23]:

- **Surface oxidation:** Carbonyl index (CI = A₁₇₁₅/A₁₄₆₅ for PE) quantifies the degree of photo-oxidative degradation, with higher CI values indicating longer environmental residence times [24].
- **Chain scission and cross-linking:** Molecular weight distribution changes, measurable by gel permeation chromatography (GPC), reflect the balance between chain scission (reducing Mw) and cross-linking (increasing Mw) during degradation [25].
- **Biofilm signatures:** Microbial colonization of the "plastisphere" produces characteristic proteinaceous and polysaccharide signatures detectable in FTIR and Raman spectra, indicating environmental exposure history [26].
- **Crystallinity changes:** Differential scanning calorimetry (DSC) and X-ray diffraction (XRD) can detect increased crystallinity resulting from preferential degradation of amorphous polymer regions [27].

These degradation markers are particularly relevant for distinguishing between recently ingested iMNPs (with minimal degradation) and particles that have undergone prolonged environmental weathering before entering the body, as the latter may carry different toxicological profiles due to altered surface chemistry and adsorbed contaminants [28].

---

## 8.3 Reference Database Development

### 8.3.1 Existing Databases and Their Limitations

Reliable identification of iMNPs requires matching experimental spectra against comprehensive reference databases. Several spectral libraries are currently available:

- **Open Specy:** An open-source platform providing FTIR and Raman reference spectra for >400 polymer types, with automated spectral matching algorithms [29]. While freely accessible, its coverage of degraded polymers, polymer blends, and biologically altered particles remains limited.
- **KnowItAll (Wiley):** A commercial database containing >260,000 FTIR spectra with hierarchical search capabilities. The extensive coverage comes at significant cost and may include spectra collected under different instrumental conditions [7].
- **slopp and slopp-e:** Specialized spectral libraries for microplastics and environmental particles, developed through community curation efforts [31].
- **Py-GC/MS Libraries:** NIST and Wiley mass spectral databases provide pyrolysis product identification, but dedicated pyrolysis product libraries for aged/degraded polymers are scarce [32].

**Key limitations of existing databases include:**

1. **Bias toward pristine polymers:** Most reference spectra are collected from virgin materials, whereas iMNPs in biological tissues are often aged, degraded, and coated with biological corona, producing altered spectral profiles [33].
2. **Lack of biological matrix spectra:** Reference libraries rarely include spectra of polymers embedded in or contaminated with biological matrices (lipids, proteins, mineral deposits), which can mask diagnostic polymer bands [34].
3. **Instrument-to-instrument variability:** Spectra collected on different instruments with different configurations (laser wavelength, detector type, resolution) may not be directly comparable, leading to reduced match quality [35].
4. **Limited nanoplastic coverage:** Databases overwhelmingly represent microplastic-scale particles; spectral characteristics of nanoplastics, which may differ due to quantum confinement effects and surface-to-volume ratio changes, are poorly represented [36].

### 8.3.2 Standardized Data Formats and Reporting Standards

The field urgently requires standardized approaches to data management and reporting:

**Data format standardization.** The adoption of open spectral data formats (JCAMP-DX for vibrational spectroscopy, mzML for mass spectrometry) with mandatory metadata fields (instrument parameters, sample preparation, matching algorithms, match quality scores) would facilitate cross-study comparisons and database interoperability [37].

**Minimum reporting standards.** Drawing from initiatives such as MIAPE (Minimum Information About a Proteomics Experiment), a consensus framework for iMNP studies should mandate reporting of: (a) detection limits (size and mass); (b) contamination control measures and blank correction procedures; (c) spectral match quality scores and acceptance thresholds; (d) recovery rates from spiked samples; and (e) complete particle characterization data (size, shape, polymer type, color) [38,39].

**Quality assurance/quality control (QA/QC) programs.** Inter-laboratory proficiency testing programs, analogous to the QUASIMEME scheme for marine contaminants, are needed to assess and harmonize analytical performance across laboratories [40]. Early inter-comparison exercises have revealed alarming variability: a recent study found that laboratories analyzing identical samples reported particle counts varying by up to 3 orders of magnitude [41].

> **Table 8.1** summarizes the principal analytical methods for iMNP detection in biological tissues, including their capabilities, limitations, and recommended applications.

### Table 8.1. Comparison of Analytical Methods for iMNP Detection in Biological Tissues

| Method | Size Range | Information Provided | Destructive? | Throughput | Key Advantages | Key Limitations | Recommended Application |
|:---|:---:|:---|:---:|:---:|:---|:---|:---|
| μ-FTIR | >3–10 μm | Polymer ID, particle morphology | No | Medium | Well-established, automated imaging, extensive libraries | Diffraction-limited resolution, water interference | Screening of MPs in tissue sections |
| ATR-FTIR | >50 μm | Polymer ID, surface chemistry | No | High | Minimal preparation, rapid | Large particle bias, single-point analysis | Rapid identification of individual large particles |
| μ-Raman | >1 μm | Polymer ID, particle morphology, crystallinity | No | Low–Medium | Higher spatial resolution, water-compatible | Fluorescence interference, slow mapping | Detailed characterization of small MPs |
| Py-GC/MS | No limit (mass-based) | Polymer type, mass concentration | Yes | Medium | NP-capable, quantitative mass data | No morphological data, destructive | Quantifying total polymer mass in tissue |
| TED-GC/MS | No limit (mass-based) | Polymer type, mass concentration | Yes | High | Larger sample volumes, NP-capable | Less sensitive for trace levels | High-throughput screening |
| SEM-EDX | >100 nm | Morphology, elemental composition | No | Low | Nanoscale imaging, elemental analysis | No polymer identification | Morphological characterization |
| TOF-SIMS | >100 nm | Surface chemistry, spatial mapping | Partially | Low | Nanoscale chemical imaging | Complex interpretation, expensive | Research-level surface analysis |
| AFM-IR | >20 nm | Nanoscale polymer ID | No | Very low | Sub-diffraction chemical ID | Extremely low throughput, specialized | Nanoplastic identification |
| LDIR | >10 μm | Polymer ID, automated counting | No | Very high | Fastest automated analysis | Lower resolution than μ-FTIR | Large-scale surveys |

---

## 8.4 Biodistribution Studies

### 8.4.1 Animal Model Evidence

Animal models have been instrumental in establishing the biodistribution patterns of iMNPs following different exposure routes. These studies, primarily conducted in rodent models, provide mechanistic insights that complement observational human data.

**Oral exposure studies.** Deng et al. (2017) administered fluorescent polystyrene (PS) microplastics (5 μm and 20 μm) to mice via gavage and demonstrated accumulation in the liver, kidney, and gut within 28 days [42]. Smaller particles (5 μm) showed significantly higher tissue accumulation than larger particles (20 μm), establishing the principle of size-dependent biodistribution. Lu et al. (2018) reported that chronic oral PS-MP exposure in mice altered gut microbiota composition and increased hepatic lipid accumulation [43]. Stock et al. (2019) further demonstrated that orally administered PS particles (1, 4, and 10 μm) could be detected in liver, kidney, and spleen tissues of mice at 28 days, though only a minor fraction of particles translocated across the intact intestinal barrier [44].

**Inhalation exposure studies.** Lu et al. (2021) exposed rats to aerosolized PS-MPs (1–5 μm) and observed deposition in the lung parenchyma, with translocation to mediastinal lymph nodes within 14 days [45]. Fournier et al. (2020) showed that acute lung exposure to 20 nm PS nanoplastics in late-stage pregnant rats resulted in translocation across the air-blood barrier and subsequent fetal deposition, demonstrating both pulmonary translocation and transplacental transfer [136]. Inhalation studies consistently show that particles <2.5 μm penetrate deep into the alveolar region, while larger particles deposit predominantly in the upper airways and are cleared via mucociliary transport [46].

**Intravenous administration studies.** Direct injection studies bypass absorption barriers and reveal the intrinsic distribution kinetics of particles within the systemic circulation. Walczak et al. (2015) demonstrated that intravenously administered PS nanoparticles (50 nm) in rats accumulated primarily in the liver (>60% of dose), spleen (~15%), and lungs (~5%), with minimal renal excretion [47]. Surface charge significantly influenced biodistribution: positively charged particles showed enhanced hepatic uptake, while PEGylated particles exhibited prolonged circulation times [48].

### 8.4.2 Organ-Specific Accumulation Patterns

Synthesizing evidence across multiple animal studies, a consistent hierarchy of iMNP accumulation emerges:

1. **Liver:** The liver consistently shows the highest accumulation of systemically distributed iMNPs, attributable to its role as the primary organ of the reticuloendothelial system (RES) and its fenestrated sinusoidal endothelium, which allows passage of particles up to ~100–200 nm [42,47,49].
2. **Gut/intestine:** The gastrointestinal tract retains the majority (>90%) of orally administered particles, most of which transit through without systemic absorption [43,50].
3. **Spleen:** Splenic accumulation reflects macrophage-mediated clearance of circulating particles [47].
4. **Lungs:** Both inhaled and systemically circulating particles can deposit in lung tissue [45,46].
5. **Kidneys:** Renal accumulation is generally low, though particles <6 nm may undergo glomerular filtration [51].
6. **Brain:** Emerging evidence demonstrates that nanoscale particles (<200 nm) can cross the blood-brain barrier, with accumulation demonstrated in rodent models [52,53].
7. **Reproductive organs:** Accumulation in testes and ovaries has been reported, with implications for reproductive toxicity [54,55].
8. **Placenta:** Transplacental transfer of nanoplastics has been demonstrated in murine models, with fetal tissue accumulation [56,57].

> **Table 8.2** summarizes key biodistribution studies in animal models.

### Table 8.2. Summary of Key Biodistribution Studies of iMNPs in Animal Models

| Study | Animal Model | Particle Type & Size | Exposure Route | Duration | Key Organs of Accumulation | Key Findings |
|:---|:---|:---|:---:|:---:|:---|:---|
| Deng et al. (2017) [42] | ICR mice | PS, 5 & 20 μm | Oral gavage | 28 days | Liver, kidney, gut | Size-dependent accumulation; 5 μm > 20 μm |
| Walczak et al. (2015) [47] | Sprague-Dawley rats | PS, 50 & 500 nm | IV injection | 24 h | Liver (>60%), spleen, lungs | Charge-dependent distribution |
| Lu et al. (2018) [43] | ICR mice | PS, 0.5 & 50 μm | Drinking water | 5 weeks | Gut, liver | Gut microbiota dysbiosis, hepatic lipid changes |
| Stock et al. (2019) [44] | C57BL/6 mice | PS, 1, 4 & 10 μm | Oral gavage | 28 days | Gut, liver, kidney, spleen | Minor fraction of particles translocated across intestinal barrier |
| Fournier et al. (2020) [136] | Sprague-Dawley rats | PS, 20 nm | Lung instillation | 24 h | Lungs, placenta, fetal tissue | Pulmonary-to-systemic translocation; transplacental transfer |
| Lu et al. (2021) [45] | Wistar rats | PS, 1–5 μm | Inhalation | 14 days | Lungs, lymph nodes | Alveolar deposition, lymphatic translocation |
| Wick et al. (2010) [56] | Ex vivo human placenta | PS, 50–500 nm | Perfusion | 6 h | Fetal compartment | Size-dependent placental transfer (<240 nm) |
| Kopatz et al. (2023) [52] | C57BL/6 mice | PS, 9.5 μm & 1.1 μm | Oral gavage | 2–4 h | Brain (short-term) | BBB crossing via surface cholesterol corona |
| Jin et al. (2021) [54] | ICR mice | PS, 0.5 & 5 μm | Oral gavage | 35 days | Testes | Reduced sperm quality, oxidative stress |

---

## 8.5 Biological Barrier Crossing

The detection of iMNPs in diverse human organs implies successful crossing of multiple biological barriers. Understanding these barrier-crossing mechanisms is critical for predicting biodistribution patterns and assessing health risks.

### 8.5.1 Intestinal Epithelial Barrier

The intestinal epithelium represents the first biological barrier encountered by ingested iMNPs. This single-cell-layer barrier comprises enterocytes connected by tight junctions, interspersed with mucus-secreting goblet cells and antigen-sampling M cells [58].

**Transcellular pathways.** Particles <1 μm can be internalized by enterocytes via endocytosis (clathrin-mediated, caveolae-mediated, or macropinocytosis), with uptake efficiency inversely proportional to particle size [59]. Nanoplastics (50–100 nm) show substantially higher transcellular permeability than microplastics. In Caco-2 monolayer models, the apparent permeability (Papp) of 100 nm PS particles was approximately 10-fold higher than that of 1 μm particles [60].

**M-cell mediated uptake.** M cells overlying Peyer's patches represent a specialized pathway for particulate uptake. M cells lack the thick glycocalyx and brush border of enterocytes and actively sample luminal particulates for immune surveillance. Studies using Caco-2/Raji B co-culture models (which induce M-cell differentiation) have shown 5–10-fold enhanced translocation of 200 nm PS particles compared to enterocyte-only monolayers [61].

**Paracellular transport.** Tight junctions normally restrict paracellular passage to molecules <1 nm in diameter. However, iMNPs may indirectly promote paracellular transport by: (a) disrupting tight junction protein expression (ZO-1, occludin, claudins) through inflammatory signaling [62]; (b) physically destabilizing membrane integrity at high local concentrations [63]; and (c) carrying adsorbed contaminants (e.g., bisphenol A, heavy metals) that independently compromise barrier function [28]. Evidence from in vitro studies suggests that nanoplastics at concentrations ≥50 μg/mL can reduce transepithelial electrical resistance (TEER) by 20–40%, indicating barrier compromise [65].

**Mucus layer interactions.** Before reaching the epithelium, ingested particles must traverse the mucus gel layer (10–400 μm thick). Particle-mucus interactions are governed by size, surface charge, and hydrophobicity. Hydrophobic particles >500 nm are efficiently trapped in mucus, while smaller hydrophilic or PEGylated particles exhibit enhanced mucus penetration [66].

### 8.5.2 Alveolar-Capillary Barrier

The alveolar-capillary barrier is the thinnest barrier in the body (~0.1–0.5 μm), consisting of type I pneumocytes, a fused basement membrane, and capillary endothelial cells [67]. This extreme thinness, designed for efficient gas exchange, also renders it vulnerable to particle translocation.

**Deposition mechanics.** Inhaled particles deposit in the respiratory tract according to their aerodynamic diameter: particles >10 μm deposit in the nasopharynx, 5–10 μm in the bronchi, 1–5 μm in the bronchioles, and <1 μm can reach the alveoli [46]. The alveolar deposition fraction peaks at ~20% for particles around 2–3 μm and increases again for ultrafine/nanoparticles (<100 nm) due to diffusion [68].

**Translocation mechanisms.** Once deposited in the alveolar space, particles may translocate across the air-blood barrier via: (a) transcytosis through type I pneumocytes; (b) uptake by alveolar macrophages followed by migration to pulmonary lymphatics; and (c) direct penetration through epithelial discontinuities or via thin-walled capillaries [69]. Particle translocation rates from the alveolar space to the bloodstream have been estimated at 1–5% for 100 nm particles in rodent models, decreasing with increasing particle size [70].

**Clearance mechanisms.** Particles deposited in the conducting airways are cleared by mucociliary transport (half-life ~24 hours), while alveolar-deposited particles are cleared primarily by macrophage-mediated phagocytosis and lymphatic drainage (half-life weeks to months) [71]. Overloading of macrophage clearance capacity can lead to "frustrated phagocytosis" and inflammatory cascades, particularly with fibrous particles [72].

### 8.5.3 Blood-Brain Barrier

The blood-brain barrier (BBB) is the most restrictive biological barrier in the body, formed by brain microvascular endothelial cells connected by continuous tight junctions, surrounded by pericytes and astrocyte endfeet [73]. The BBB selectively excludes >98% of small molecules and essentially all macromolecules from the brain parenchyma.

Recent studies have challenged the assumption that iMNPs cannot cross the BBB. Kopatz et al. (2023) demonstrated that orally administered PS microparticles (9.5 μm) could be detected in mouse brain tissue within 2 hours of gavage, with the mechanism linked to cholesterol adsorption forming a "biomolecular corona" that facilitated receptor-mediated transcytosis via low-density lipoprotein receptor (LDLR) [52]. Shan et al. (2022) showed that 80 nm PS nanoplastics accumulated in zebrafish brain tissue and induced neurotoxicity through microglial activation [74]. The recent detection of microplastics in human brain tissue by Nihart et al. (2025), with concentrations significantly higher than in liver and kidney, provides direct human evidence of BBB penetration by plastic particles [53].

Potential BBB-crossing mechanisms for iMNPs include: (a) receptor-mediated transcytosis (exploiting LDLR, transferrin receptor, or other endocytic receptors via biomolecular corona formation); (b) adsorptive-mediated transcytosis (cationic surface charges interacting with anionic endothelial glycocalyx); and (c) circumventricular organ entry (brain regions lacking complete BBB) [75].

### 8.5.4 Placental Barrier

The human placental barrier separates maternal and fetal circulation through syncytiotrophoblast layers, a basement membrane, and fetal capillary endothelium. The discovery of microplastics in human placental tissue by Ragusa et al. (2021) raised profound concerns about prenatal exposure [57].

**Ex vivo perfusion studies.** The ex vivo human placental perfusion model, first applied to nanoplastics by Wick et al. (2010), demonstrated size-dependent transplacental transfer: PS particles ≤240 nm crossed the placental barrier within 6 hours, while 500 nm particles did not [56]. Subsequent studies confirmed that surface charge and functionalization modulate placental transfer efficiency, with carboxylated particles showing enhanced maternal-to-fetal transport compared to unmodified or aminated variants [77].

**In vivo evidence.** Halfar et al. (2023) detected microplastics in 9 out of 10 human placentas (90%) and, notably, in amniotic fluid, demonstrating that particles reach the fetal compartment [78]. The meta-analysis presented in Section 8.9 of this chapter confirms the high detection prevalence in placental tissue (90.0%, 95% CI: 59.6–98.2%).

**Implications for fetal development.** Transplacental transfer of iMNPs raises concerns about exposure during critical developmental windows. Emerging evidence links prenatal MP/NP exposure to altered placental function (reduced human chorionic gonadotropin secretion), changes in meconium microbiome composition, and potential epigenetic modifications [79,80].

---

## 8.6 Translocation Mechanisms

### 8.6.1 Transcellular Transport (Endocytosis)

Transcellular transport of iMNPs is primarily mediated by endocytic pathways. The specific endocytic mechanism engaged depends on particle size, shape, surface charge, and the cell type encountered [81]:

- **Clathrin-mediated endocytosis (CME):** The principal pathway for particles 50–200 nm, involving clathrin coat assembly around membrane invaginations forming ~120 nm vesicles. CME has been shown to be the dominant uptake pathway for PS nanoplastics in intestinal epithelial cells [81].
- **Caveolae-mediated endocytosis:** Involves 50–100 nm flask-shaped membrane invaginations enriched in caveolin-1. This pathway can bypass lysosomal degradation, potentially enabling intact transcytosis of nanoparticles [83].
- **Macropinocytosis:** A non-specific pathway for larger particles (>200 nm–5 μm), driven by actin-dependent membrane ruffling. Macropinocytosis is particularly active in macrophages and dendritic cells [84].
- **Phagocytosis:** Professional phagocytes (macrophages, neutrophils) can internalize particles up to 10 μm through receptor-mediated phagocytosis. Opsonization of iMNPs by complement proteins and immunoglobulins enhances phagocytic uptake [85].

The biomolecular corona—the layer of proteins, lipids, and other biomolecules that rapidly adsorbs onto iMNP surfaces upon entering biological fluids—profoundly influences which endocytic pathway is engaged. Corona composition varies by biological compartment (blood, lung surfactant, gastrointestinal fluid) and can either enhance or inhibit cellular uptake [86,87].

### 8.6.2 Paracellular Transport (Tight Junction Disruption)

Paracellular transport through intercellular tight junctions is normally restricted to ions and small molecules (<1 nm). However, iMNPs can indirectly promote paracellular permeability through several mechanisms [88]:

- **Inflammatory cytokine-mediated disruption:** iMNPs activate innate immune signaling (NF-κB, NLRP3 inflammasome), leading to release of pro-inflammatory cytokines (TNF-α, IL-1β, IL-6) that downregulate tight junction protein expression [62,89].
- **Reactive oxygen species (ROS):** Oxidative stress induced by iMNPs can directly damage tight junction complexes and promote epithelial barrier dysfunction [90].
- **Direct physical disruption:** At high concentrations, nanoplastics may physically insert into lipid bilayers and disrupt membrane-associated tight junction anchoring [91].

### 8.6.3 M-Cell Mediated Transport

M (microfold) cells, found in the follicle-associated epithelium of Peyer's patches and isolated lymphoid follicles, represent a specialized transcytotic pathway for particulate antigens [92]. Key features of M-cell mediated iMNP transport include:

- **Enhanced particle uptake capacity:** M cells actively sample particles up to 10 μm from the intestinal lumen, with preferential uptake of hydrophobic particles in the 50 nm–5 μm range [93].
- **Basolateral delivery to immune cells:** Transcytosed particles are released to subepithelial dendritic cells and macrophages, potentially entering lymphatic circulation and reaching mesenteric lymph nodes [94].
- **Quantitative significance:** M cells constitute only ~1% of intestinal epithelial cells but may account for >50% of total particulate translocation across the intestinal barrier [95].

> **Figure 8.1** (suggested): Schematic illustration of iMNP translocation mechanisms across the intestinal epithelial barrier, showing transcellular (endocytosis), paracellular (tight junction disruption), and M-cell mediated pathways.

---

## 8.7 Retention and Clearance

### 8.7.1 Half-Life Estimates

The biological half-life of iMNPs—the time required for elimination of 50% of the tissue burden—varies dramatically depending on particle size, polymer type, tissue location, and host factors [96]. Current estimates, derived primarily from animal studies, include:

- **Gastrointestinal transit:** Half-life of 24–48 hours for orally ingested particles that are not systemically absorbed (>99% of ingested particles) [50].
- **Pulmonary clearance:** Alveolar macrophage-mediated clearance has a half-life of weeks to months. In industrial hygiene studies of fibrous particles, pulmonary retention half-lives of 1–3 years have been documented [71,97].
- **Systemic circulation:** Particles in the bloodstream are rapidly cleared by the reticuloendothelial system (primarily liver and spleen), with circulating half-lives estimated at minutes to hours for uncoated particles [48].
- **Tissue retention:** Once sequestered in organs (liver, spleen, lymph nodes), iMNPs may persist for months to years in the absence of enzymatic degradation, as mammalian biology lacks enzymes capable of depolymerizing synthetic plastics [98].

The concept of "half-life" for iMNPs is complicated by the potential for continuous re-exposure (daily ingestion and inhalation), creating a steady-state body burden rather than a simple clearance curve [99].

### 8.7.2 Excretion Routes

The principal excretion routes for iMNPs include:

- **Fecal excretion:** The dominant elimination pathway for ingested particles, accounting for >90% of orally administered MP/NP dose in animal studies [50]. The detection of microplastics in human stool (pooled prevalence 85.3% in our meta-analysis) directly reflects this fecal excretion pathway.
- **Biliary excretion:** Hepatocytes can excrete particles into bile via canalicular transport, contributing to fecal elimination of systemically distributed particles [100].
- **Renal excretion:** Limited to very small nanoparticles (<6–8 nm) that can pass through the glomerular filtration barrier [51,101].
- **Mucociliary clearance:** Inhaled particles deposited in the conducting airways are transported upward by ciliary action and swallowed, ultimately contributing to fecal excretion [71].
- **Cerumen secretion:** The recent detection of microplastics in human cerumen (83.3% detection rate) by Arslan et al. (2025) suggests ear wax as a novel, previously unrecognized excretion/elimination pathway [102].
- **Breast milk excretion:** Detection of MPs in breast milk (39.0% in our meta-analysis) indicates mammary secretion as both an excretion route and a potential neonatal exposure pathway [103].

### 8.7.3 Factors Affecting Persistence

Several physicochemical and biological factors modulate the in vivo persistence of iMNPs:

- **Particle size:** Smaller particles (<1 μm) show greater tissue penetration and longer retention times compared to larger particles that are more efficiently cleared by phagocytosis [42].
- **Polymer type:** Chemical resistance to biological degradation varies by polymer. PE and PP are highly resistant to enzymatic attack, while polyesters (PET, PLA) may undergo slow hydrolysis in tissue environments [104].
- **Surface chemistry:** Hydrophobic, unmodified particles are more rapidly opsonized and cleared than hydrophilic or PEGylated particles [48].
- **Particle shape:** Fibrous particles (aspect ratio >3:1) resist phagocytic clearance more effectively than spherical particles, potentially leading to greater tissue persistence and inflammatory responses (the "frustrated phagocytosis" paradigm) [72,105].
- **Tissue location:** Particles sequestered in immune-privileged sites (brain, testes, vitreous humor) or within phagolysosomes may experience minimal clearance [98].
- **Age and health status:** Compromised macrophage function in elderly individuals or immunosuppressed patients may reduce clearance capacity, leading to greater tissue accumulation [107].

---

## 8.8 Kinetic Modeling

### 8.8.1 Physiologically-Based Pharmacokinetic (PBPK) Models

PBPK modeling provides a mathematical framework for predicting the absorption, distribution, metabolism, and excretion (ADME) of substances based on physiological parameters (organ volumes, blood flows, partition coefficients) and substance-specific properties [108]. Originally developed for pharmaceutical compounds, PBPK models are increasingly adapted for iMNP kinetics.

A typical PBPK model for iMNPs comprises interconnected compartments representing key organs (liver, lungs, spleen, kidneys, gut, brain, etc.), connected by blood flow. Each compartment is described by a mass balance ordinary differential equation (ODE):

$$\frac{dA_i}{dt} = Q_i \cdot (C_{art} - \frac{C_i}{P_i}) - k_{clearance,i} \cdot A_i$$

where $A_i$ is the amount in tissue $i$, $Q_i$ is blood flow, $C_{art}$ is arterial concentration, $C_i$ is tissue concentration, $P_i$ is the tissue-blood partition coefficient, and $k_{clearance,i}$ is the local clearance rate constant [109].

### 8.8.2 Adapting PBPK for Particles

Applying PBPK models to particulate materials (as opposed to dissolved molecules) introduces several unique challenges [110]:

1. **Size-dependent kinetics:** Unlike dissolved molecules that follow partition equilibrium, particle behavior is governed by filtration, phagocytosis, and physical trapping, all of which are strongly size-dependent [111].
2. **Reticuloendothelial system (RES) uptake:** Hepatic and splenic capture of particles by resident macrophages (Kupffer cells, splenic macrophages) must be modeled explicitly, often using saturable (Michaelis-Menten) uptake kinetics rather than simple partition coefficients [112].
3. **Surface-dependent interactions:** The biomolecular corona that forms on iMNP surfaces alters their apparent physicochemical properties (effective size, charge, receptor interactions), requiring dynamic surface chemistry parameters [86].
4. **Barrier-crossing kinetics:** Particle translocation across biological barriers (intestinal, alveolar, placental, BBB) follows different kinetics than molecular diffusion, often requiring permeability rate constants calibrated from ex vivo or in vitro barrier models [113].
5. **Slow clearance/no metabolism:** Unlike pharmaceuticals that are enzymatically metabolized, most iMNPs are not biodegradable on biologically relevant timescales, leading to accumulation over time [98].

**Published PBPK models for nanoparticles.** Several PBPK models have been developed for engineered nanoparticles (gold, silver, titanium dioxide) that provide templates for iMNP modeling [114,115]. Bachler et al. (2013) developed a comprehensive PBPK model for TiO₂ nanoparticles incorporating macrophage uptake kinetics that was subsequently adapted for nanoplastics [116]. Mohamed Nor et al. (2021) developed a lifetime accumulation model for MPs in humans, estimating that adults may accumulate ~40 g of microplastic over a lifetime under current exposure scenarios [99].

> **Figure 8.2** (suggested): Schematic of a PBPK model for iMNP biodistribution, showing organ compartments, blood flows, barrier-crossing parameters, and clearance pathways.

---

## 8.9 Human Evidence: A Systematic Review and Meta-Analysis

This section presents original meta-analysis results based on a comprehensive systematic review of the literature on iMNP detection in human tissues. The analysis follows PRISMA 2020 guidelines and employs random-effects meta-analysis to provide pooled prevalence estimates [117].

### 8.9.1 Systematic Review Methods

**Search strategy.** A systematic search was conducted across PubMed/MEDLINE (n = 3,000 records), Europe PMC (n = 2,000), and OpenAlex (n = 274) databases through January 2026, yielding 5,274 total records. After deduplication (n = 792 duplicates removed), 4,482 unique records were screened. Following title/abstract screening, 1,946 studies were included in qualitative synthesis (bibliometric analysis), and 11 studies with extractable quantitative prevalence data (total N = 474 participants) were included in the random-effects meta-analysis [117].

**Statistical methods.** Proportions were transformed using the Freeman-Tukey double arcsine transformation. Random-effects pooling was performed using the DerSimonian-Laird method. Heterogeneity was assessed using the I² statistic and Cochran's Q test. Publication bias was evaluated via funnel plot inspection and Egger's regression test. Leave-one-out sensitivity analysis assessed the robustness of pooled estimates [118,119].

### 8.9.2 Blood

Blood serves as the central distribution medium for systemically absorbed iMNPs. Leonard et al. (2024) analyzed whole blood samples from 20 healthy volunteers using μ-FTIR and detected microplastics in 18 out of 20 donors (90.0%) [120]. Twenty-four polymer types were identified, with PE (32%), ethylene propylene diene (14%), and ethylene-vinyl-acetate (12%) being the most abundant. Concentrations ranged from 1.84 to 4.65 μg/mL, and the mean particle length was 128 μm (range: 7–3,000 μm). This study builds upon the landmark discovery by Leslie et al. (2022), who first demonstrated quantifiable plastic particle concentrations in human blood [121].

The pooled prevalence for blood/serum/plasma was **90.0% (95% CI: 69.9–97.2%)**, though based on a single meta-analyzed study (k = 1).

### 8.9.3 Lungs

The lungs are directly exposed to airborne iMNPs through inhalation. Three studies in our meta-analysis examined lung tissue or bronchoalveolar lavage (BAL) fluid:

- Padarya et al. (2025) detected MPs in 42 of 60 BAL samples (70.0%) using FTIR, with PE and PP being the most common polymers and a mean count of 6.4 ± 3.1 particles per 100 mL BAL fluid [122].
- Özgen Alpaydin et al. (2024) found MPs in 10 of 18 BAL samples (55.6%) from patients with suspected interstitial lung disease, using μ-Raman spectroscopy [123].
- Zakynthinos et al. (2025) detected MPs in 5 of 8 BAL samples (62.5%) via SEM-EDX, with particles ranging from 5.9 to 204.7 μm [124].

The pooled prevalence for lung tissue was **66.1% (95% CI: 55.4–75.4%; I² = 0%)**, with notably low within-tissue heterogeneity despite different analytical methods, suggesting consistent underlying detection rates.

### 8.9.4 Liver

Although the liver is a primary site of iMNP accumulation in animal models, no studies in our meta-analysis specifically quantified hepatic iMNP prevalence with binary detection data. Horvatits et al. (2022) detected MPs in cirrhotic liver tissue using μ-FTIR, reporting concentrations up to 84 MPs/g tissue, but did not report a binary detection rate suitable for meta-analytic pooling [125]. The bibliometric analysis identified 100 publications (5.1% of 1,946) investigating liver tissue, confirming growing research interest. This gap highlights the need for liver-focused prevalence studies with standardized reporting.

### 8.9.5 Placenta

Placental tissue showed one of the highest detection prevalence rates. Halfar et al. (2023) detected MPs or plastic additives in 9 out of 10 placentas (90.0%) from patients with preterm premature rupture of membranes, as well as in amniotic fluid, using ATR-FTIR [78]. Chlorinated polyethylene and calcium zinc PVC stabilizer were the dominant materials, with particle sizes between 10 and 50 μm.

The pooled prevalence for placenta was **90.0% (95% CI: 59.6–98.2%)**, though based on a single study (k = 1).

### 8.9.6 Stool (Fecal Samples)

Two studies examined iMNP prevalence in stool samples:

- Ke et al. (2023) detected MPs in 59 of 69 stool samples (85.5%) from preschool children in Xiamen, China, using Py-GC/MS [126]. PVC, PET, PE, and polyamide-6 were detected, with median concentrations of 317.4 μg/g for PVC.
- Santini et al. (2024) detected MPs in 5 of 6 stool samples (83.3%) using FTIR, identifying synthetic cellulose, PE, PP, PS, and PET [127].

The pooled prevalence for stool/gut was **85.3% (95% CI: 75.4–91.7%; I² = 0%)**, with excellent consistency between the two studies.

### 8.9.7 Other Tissues

Additional tissues examined in the meta-analysis included:

- **Cerumen (ear wax):** Arslan et al. (2025) detected MPs in 10 of 12 adult patients (83.3%), with particles ranging from 16 to 930 μm. This novel finding suggests cerumen as a non-invasive biomarker for MP exposure [102].
- **Semen/reproductive tissue:** Two studies (Demirelli et al., 2024; Qu et al., 2025) examined reproductive tissues, with detection rates of 50.0% and 55.5%, yielding a pooled prevalence of **55.2% (95% CI: 48.4–61.8%; I² = 0%)** [128,129].
- **Breast milk:** Saraluck et al. (2024) detected MPs in 23 of 59 breast milk samples (39.0%), with PP, PE, PS, and PVC being the most common polymers [103].

### 8.9.8 Pooled Meta-Analysis Results

The overall random-effects meta-analysis across 11 studies (N = 474) yielded a pooled iMNP detection prevalence of **68.5% (95% CI: 56.3–78.6%)** with substantial heterogeneity (I² = 77%, Q = 43.3, p < 0.001). The 95% prediction interval was 28.9–92.1%, reflecting wide expected variability in future studies.

Subgroup analysis by tissue type revealed that between-tissue differences were the primary driver of heterogeneity: within-tissue I² was 0% for all subgroups with ≥2 studies (lung, stool, semen). Egger's regression test showed no significant publication bias (intercept = 1.72, p = 0.909). Leave-one-out sensitivity analysis confirmed robust pooled estimates (range: 64.0–72.0%) [117].

> **Table 8.3** presents the tissue-specific pooled prevalence estimates from the meta-analysis.

### Table 8.3. Pooled Prevalence of iMNPs in Human Tissues: Random-Effects Meta-Analysis Results

| Tissue Type | No. of Studies (k) | Total N | Pooled Prevalence (%) | 95% CI (%) | I² (%) | Q Statistic (p-value) |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| Blood/serum/plasma | 1 | 20 | 90.0 | 69.9–97.2 | 0 | — |
| Placenta | 1 | 10 | 90.0 | 59.6–98.2 | 0 | — |
| Stool/gut | 2 | 75 | 85.3 | 75.4–91.7 | 0 | 0.02 (0.89) |
| Cerumen/ear | 1 | 12 | 83.3 | 55.2–95.3 | 0 | — |
| Lung | 3 | 86 | 66.1 | 55.4–75.4 | 0 | 1.33 (0.51) |
| Semen/testis | 2 | 212 | 55.2 | 48.4–61.8 | 0 | 0.14 (0.71) |
| Breast milk | 1 | 59 | 39.0 | 27.6–51.7 | 0 | — |
| **Overall** | **11** | **474** | **68.5** | **56.3–78.6** | **77** | **43.3 (<0.001)** |

*CI = confidence interval. Random-effects meta-analysis using DerSimonian-Laird method with Freeman-Tukey double arcsine transformation. Wilson score confidence intervals for individual studies.*

### 8.9.9 Quality Assessment of Human Evidence

Quality assessment using the Modified Newcastle-Ottawa Scale revealed substantial limitations across the evidence base. Among the 11 meta-analyzed studies, the mean quality score was 3.4/9 (SD: 2.0), with 7 studies rated as moderate quality (4–6) and 4 as low quality (0–3). No studies achieved high quality (7–9). The most commonly unmet criteria were method validation (O1; met by 0/11 studies), recovery rate reporting (O3; 0/11), confounder adjustment (C1; 1/11), and non-respondent characterization (S3; 2/11). Contamination control (C2), a critical criterion for iMNP studies, was met by 7/11 studies [117].

> **Table 8.4** presents the quality assessment scores for the 11 meta-analyzed studies.

### Table 8.4. Quality Assessment of Studies Included in Meta-Analysis (Modified Newcastle-Ottawa Scale)

| Study | Year | S1 | S2 | S3 | S4 | C1 | C2 | O1 | O2 | O3 | Total (/9) | Grade |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Leonard et al. [120] | 2024 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 2 | Low |
| Ke et al. [126] | 2023 | 1 | 1 | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 5 | Moderate |
| Arslan et al. [102] | 2025 | 1 | 1 | 1 | 0 | 0 | 1 | 0 | 1 | 0 | 5 | Moderate |
| Santini et al. [127] | 2024 | 1 | 0 | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 4 | Moderate |
| Zakynthinos et al. [124] | 2025 | 1 | 1 | 0 | 0 | 0 | 1 | 0 | 1 | 0 | 4 | Moderate |
| Demirelli et al. [128] | 2024 | 1 | 1 | 0 | 1 | 0 | 1 | 0 | 1 | 0 | 5 | Moderate |
| Saraluck et al. [103] | 2024 | 0 | 1 | 1 | 0 | 0 | 1 | 0 | 1 | 0 | 4 | Moderate |
| Özgen Alpaydin et al. [123] | 2024 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | Low |
| Halfar et al. [78] | 2023 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | Low |
| Qu et al. [129] | 2025 | 1 | 1 | 0 | 1 | 1 | 1 | 0 | 1 | 0 | 6 | Moderate |
| Padarya et al. [122] | 2025 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 2 | Low |

*S1 = representativeness; S2 = sample size adequacy; S3 = non-respondent characterization; S4 = ascertainment method; C1 = confounder adjustment; C2 = contamination control; O1 = method validation; O2 = statistical analysis; O3 = recovery rate reported.*

### 8.9.10 Research Limitations

Several critical limitations of the current human evidence should be acknowledged:

1. **Small meta-analytic sample:** Only 11 studies met inclusion criteria, reflecting the field's immaturity and heterogeneous reporting practices.
2. **Tissue representation imbalance:** Some tissues (blood, placenta, cerumen, breast milk) were represented by only a single study, precluding robust tissue-specific estimates.
3. **Contamination bias:** Despite contamination control measures, complete elimination of background plastic contamination during sample collection and processing remains technically challenging.
4. **Detection limit heterogeneity:** Different analytical methods have different lower size limits (ATR-FTIR: ~50 μm; μ-FTIR: ~3 μm; μ-Raman: ~1 μm; Py-GC/MS: no size limit), meaning "prevalence" captures different subsets of the true particle burden.
5. **Geographic and temporal clustering:** Most studies were from Turkey (k = 3), China (k = 2), or Europe, published between 2023 and 2025.
6. **Absence of liver data:** Despite being a key accumulation organ in animal models, no study meeting meta-analytic inclusion criteria reported hepatic prevalence data.

> **Figure 8.3** (refer to FIG2_forest_plot.png): Forest plot of random-effects meta-analysis showing detection prevalence of iMNPs across human tissues. Individual study estimates with 95% CIs are shown as blue squares (size proportional to study weight). Tissue subgroup summaries are shown as red diamonds. The overall pooled estimate diamond spans 56.3–78.6%.

> **Figure 8.4** (refer to FIG3_funnel_plot.png): Funnel plot for publication bias assessment. No significant asymmetry was observed (Egger's test: intercept = 1.72, p = 0.909).

---

## 8.10 Comparative Analysis

### 8.10.1 Species Differences

Cross-species comparisons reveal both conserved and divergent patterns of iMNP biodistribution:

**Rodent vs. human:** Rodent models generally show qualitatively similar organ distribution patterns (liver > gut > spleen > lung), but quantitative translocation rates differ due to anatomical and physiological differences. The rodent gut has a shorter transit time (~6–8 hours vs. 24–72 hours in humans), potentially reducing absorption time. Conversely, the higher metabolic rate and cardiac output in rodents may enhance particle distribution to peripheral tissues [130].

**Aquatic vs. terrestrial organisms:** Marine organisms (fish, bivalves, crustaceans) accumulate MPs primarily in the gastrointestinal tract and gills, with limited systemic distribution [131]. In contrast, terrestrial mammals show more extensive systemic distribution following oral exposure, likely reflecting differences in gut physiology and vascular organization.

**Non-human primates:** Data on iMNP biodistribution in non-human primates are extremely limited. Given their closer physiological similarity to humans (including placental structure, BBB characteristics, and gut transit time), primate studies would provide the most translatable biodistribution data, but ethical considerations limit such research [132].

**Allometric scaling considerations:** Interspecies extrapolation of iMNP kinetics requires allometric scaling of physiological parameters (organ blood flows, cardiac output, organ volumes) and consideration of species-specific differences in macrophage function, gut physiology, and barrier permeability [133].

### 8.10.2 Particle Type Differences

The biodistribution and tissue retention of iMNPs vary significantly with polymer type:

- **Polyethylene (PE) and polypropylene (PP):** As the most produced polymers globally, PE and PP are the most frequently detected in human tissues [117]. Their high hydrophobicity promotes protein corona formation and macrophage uptake.
- **Polystyrene (PS):** The most widely used polymer in experimental studies due to availability of fluorescent and size-calibrated standards. PS shows relatively high cellular uptake efficiency and has been detected in blood and gastrointestinal samples [121,126].
- **Polyethylene terephthalate (PET):** A major component of beverage containers and textiles, PET is frequently detected in blood and gut samples. Its aromatic backbone may undergo slow hydrolysis under physiological conditions [120].
- **Polyvinyl chloride (PVC):** PVC particles, often containing plasticizer additives, show high detection rates in stool samples (linked to food packaging exposure) and have been implicated in gut microbiota disruption [126].
- **Nylon/polyamide (PA):** Textile fibers represent a major source of PA microplastics in indoor air and, consequently, in lung tissue [134].

### 8.10.3 Exposure Route Differences

The route of exposure fundamentally determines the initial tissue distribution and subsequent biodistribution pattern:

**Oral exposure:** The dominant exposure route (~95% of total daily MP intake). The majority of orally ingested particles transit the GI tract and are excreted in feces. The fraction that crosses the intestinal barrier (estimated at 0.04–0.3% for 100 nm particles) enters portal circulation and is largely sequestered by the liver on first pass [99].

**Inhalation exposure:** Contributing an estimated 0.05–0.3 g/week to total exposure. Inhaled particles deposit along the respiratory tract in a size-dependent manner, with fine and nanoparticles reaching the alveoli and potentially entering systemic circulation [45,46].

**Dermal exposure:** The least studied route, dermal absorption of MPs >1 μm is considered negligible due to the stratum corneum barrier. However, nanoplastics (<100 nm) may penetrate compromised or inflamed skin, and dermal exposure to plastic additives (plasticizers, UV stabilizers) from clothing and personal care products represents an additional concern [135].

> **Figure 8.5** (suggested): Anatomical diagram showing the distribution of iMNPs in the human body, with tissue-specific detection prevalence from the meta-analysis. Blood (90.0%), placenta (90.0%), stool/gut (85.3%), cerumen (83.3%), lung (66.1%), semen (55.2%), breast milk (39.0%), with exposure routes (oral, inhalation, dermal) indicated.

---

## 8.11 Conclusions and Future Perspectives

This chapter has provided a comprehensive review of the identification and tracking tools for iMNPs in biological systems, spanning analytical methodology, biodistribution, barrier crossing, kinetics, and human evidence.

**Key conclusions:**

1. **Analytical capabilities are rapidly expanding** but remain constrained by fundamental limitations. FTIR and Raman spectroscopy provide robust polymer identification for microplastics (>1–3 μm), while Py-GC/MS extends detection to nanoplastics in mass-based terms. No single technique provides complete characterization across the full size range.

2. **Reference databases require expansion** to include spectra of degraded, biologically altered, and matrix-embedded polymers. Standardized reporting frameworks and inter-laboratory proficiency testing are urgently needed.

3. **Animal studies demonstrate consistent biodistribution patterns** with size-dependent accumulation primarily in the liver, gut, spleen, and lungs. Nanoplastics show greater tissue penetration and longer retention than microplastics.

4. **Biological barriers are not impermeable to iMNPs.** The intestinal epithelium, alveolar-capillary barrier, blood-brain barrier, and placenta all permit particle translocation through transcellular, paracellular, and specialized (M-cell) pathways, with efficiency inversely related to particle size.

5. **PBPK modeling is a promising but underdeveloped tool** for predicting iMNP kinetics, requiring adaptation from molecular-level pharmacokinetics to particle-specific behaviors including RES uptake, barrier crossing, and absence of metabolism.

6. **Our meta-analysis of 11 studies (N = 474) demonstrates that iMNPs are detectable in 68.5% of human tissue samples** (95% CI: 56.3–78.6%), with the highest prevalence in blood (90.0%), placenta (90.0%), stool (85.3%), and cerumen (83.3%), and the lowest in breast milk (39.0%).

7. **Substantial heterogeneity (I² = 77%) reflects methodological variability** rather than true within-tissue differences, as within-tissue heterogeneity was consistently 0% for subgroups with multiple studies.

**Future priorities include:** (a) development and validation of nanoplastics-specific detection methods for biological matrices; (b) large-scale, multi-center biomonitoring studies with standardized protocols (n > 100 per tissue type); (c) PBPK model validation using human pharmacokinetic data; (d) longitudinal studies tracking individual body burden over time; and (e) integration of exposure data with health outcome epidemiology to establish dose-response relationships.

---

## References

1. Thompson RC, Olsen Y, Mitchell RP, Davis A, Rowland SJ, John AWG, et al. Lost at sea: where is all the plastic? Science. 2004;304(5672):838. doi:10.1126/science.1094559

2. Plastics Europe. Plastics – the fast facts 2023. Brussels: Plastics Europe; 2023. Available from: https://plasticseurope.org/knowledge-hub/plastics-the-fast-facts-2023/

3. Löder MGJ, Gerdts G. Methodology used for the detection and identification of microplastics—a critical appraisal. In: Bergmann M, Gutow L, Klages M, editors. Marine Anthropogenic Litter. Cham: Springer; 2015. p. 201–227. doi:10.1007/978-3-319-16510-3_8

4. Araujo CF, Nolasco MM, Ribeiro AMP, Ribeiro-Claro PJA. Identification of microplastics using Raman spectroscopy: latest developments and future prospects. Water Res. 2018;142:426–440. doi:10.1016/j.watres.2018.05.060

5. Shim WJ, Hong SH, Eo SE. Identification methods in microplastic analysis: a review. Anal Methods. 2017;9(9):1384–1391. doi:10.1039/C6AY02558G

6. Xu JL, Thomas KV, Luo Z, Gowen AA. FTIR and Raman imaging for microplastics analysis: state of the art, challenges and prospects. TrAC Trends Anal Chem. 2019;119:115629. doi:10.1016/j.trac.2019.115629

7. Primpke S, Wirth M, Lorenz C, Gerdts G. Reference database design for the automated analysis of microplastic samples based on Fourier transform infrared (FTIR) spectroscopy. Anal Bioanal Chem. 2018;410(21):5131–5141. doi:10.1007/s00216-018-1156-x

8. Jenner LC, Rotchell JM, Bennett RT, Cowen M, Sherburn V, Sherburn V, et al. Detection of microplastics in human lung tissue using μFTIR spectroscopy. Sci Total Environ. 2022;831:154907. doi:10.1016/j.scitotenv.2022.154907

9. Lenz R, Enders K, Stedmon CA, Mackenzie DMA, Nielsen TG. A critical assessment of visual identification of marine microplastic using Raman spectroscopy for analysis improvement. Mar Pollut Bull. 2015;100(1):82–91. doi:10.1016/j.marpolbul.2015.09.026

10. Käppler A, Fischer D, Oberbeckmann S, Schwabe G, Fischer F, Alanda-Zöller J, et al. Analysis of environmental microplastics by vibrational microspectroscopy: FTIR, Raman or both? Anal Bioanal Chem. 2016;408(29):8377–8391. doi:10.1007/s00216-016-9956-3

11. Lv L, He L, Jiang S, Chen J, Zhou C, Qu J, et al. In situ surface-enhanced Raman spectroscopy for detecting microplastics and nanoplastics in aquatic environments. Sci Total Environ. 2020;728:138449. doi:10.1016/j.scitotenv.2020.138449

12. Fischer M, Scholz-Böttcher BM. Simultaneous trace identification and quantification of common types of microplastics in environmental samples by pyrolysis-gas chromatography-mass spectrometry. Environ Sci Technol. 2017;51(9):5052–5060. doi:10.1021/acs.est.6b06362

13. Dümichen E, Eisentraut P, Bannick CG, Barthel AK, Campos R, Modern U. Fast identification of microplastics in complex environmental samples by a thermal degradation method. Chemosphere. 2017;174:572–584. doi:10.1016/j.chemosphere.2017.02.010

14. Schwaferts C, Niessner R, Elsner M, Ivleva NP. Methods for the analysis of submicrometer- and nanoplastic particles in the environment. TrAC Trends Anal Chem. 2019;112:52–65. doi:10.1016/j.trac.2018.12.014

15. Primpke S, Godejohann M, Gerdts G. Rapid identification and quantification of microplastics in the environment by quantum cascade laser-based hyperspectral infrared chemical imaging. Environ Sci Technol. 2020;54(24):15893–15903. doi:10.1021/acs.est.0c05722

16. Dazzi A, Prater CB. AFM-IR: technology and applications in nanoscale infrared spectroscopy and chemical imaging. Chem Rev. 2017;117(7):5146–5173. doi:10.1021/acs.chemrev.6b00448

17. Hahladakis JN, Velis CA, Weber R, Iacovidou E, Purnell P. An overview of chemical additives present in plastics: migration, release, fate and environmental impact during their use, disposal and recycling. J Hazard Mater. 2018;344:179–199. doi:10.1016/j.jhazmat.2017.10.014

18. Giuliani A, Zuccarini M, Cichelli A, Khan H, Reale M. Critical review on the presence of phthalates in food and evidence of their biological impact. Int J Environ Res Public Health. 2020;17(16):5655. doi:10.3390/ijerph17165655

19. Alaee M, Arias P, Sjödin A, Bergman Å. An overview of commercially used brominated flame retardants, their applications, their use patterns in different countries/regions and possible modes of release. Environ Int. 2003;29(6):683–689. doi:10.1016/S0160-4120(03)00121-1

20. Wick A, Fink G, Ternes TA. Comparison of electrospray ionization and atmospheric pressure chemical ionization for multi-residue analysis of biocides, UV-filters and benzothiazoles in aqueous matrices and activated sludge by liquid chromatography-tandem mass spectrometry. J Chromatogr A. 2010;1217(14):2088–2103. doi:10.1016/j.chroma.2010.01.079

21. Groh KJ, Backhaus T, Carney-Almroth B, Geueke B, Inostroza PA, Linnber A, et al. Overview of known plastic packaging-associated chemicals and their hazards. Sci Total Environ. 2019;651(Pt 2):3253–3268. doi:10.1016/j.scitotenv.2018.10.015

22. Wang Z, Walker GW, Muir DCG, Nagatani-Yoshida K. Toward a global understanding of chemical pollution: a first comprehensive analysis of national and regional chemical inventories. Environ Sci Technol. 2020;54(5):2575–2584. doi:10.1021/acs.est.9b06379

23. Gewert B, Plassmann MM, MacLeod M. Pathways for degradation of plastic polymers floating in the marine environment. Environ Sci Process Impacts. 2015;17(9):1513–1521. doi:10.1039/C5EM00207A

24. Brandon J, Goldstein M, Ohman MD. Long-term aging and degradation of microplastic particles: comparing in situ oceanic and experimental weathering patterns. Mar Pollut Bull. 2016;110(1):299–308. doi:10.1016/j.marpolbul.2016.06.048

25. Andrady AL. The plastic in microplastics: a review. Mar Pollut Bull. 2017;119(1):12–22. doi:10.1016/j.marpolbul.2017.01.082

26. Zettler ER, Mincer TJ, Amaral-Zettler LA. Life in the "plastisphere": microbial communities on plastic marine debris. Environ Sci Technol. 2013;47(13):7137–7146. doi:10.1021/es401288x

27. Song YK, Hong SH, Jang M, Han GM, Jung SW, Shim WJ. Combined effects of UV exposure duration and mechanical abrasion on microplastic fragmentation by polymer type. Environ Sci Technol. 2017;51(8):4368–4376. doi:10.1021/acs.est.6b06155

28. Hartmann NB, Rist S, Bodin J, Jensen LH, Schmidt SN, Mayer P, et al. Microplastics as vectors for environmental contaminants: exploring sorption, desorption, and transfer to biota. Integr Environ Assess Manag. 2017;13(3):488–493. doi:10.1002/ieam.1904

29. Cowger W, Steinmetz Z, Gray A, Munno K, Lynch J, Hapich H, et al. Microplastic spectral classification needs an open source community: open specy to the rescue! Anal Chem. 2021;93(21):7543–7548. doi:10.1021/acs.analchem.1c00123

30. De Frond HL, Rubinovitz R, Rochman CM. μATR-FTIR spectral libraries of plastic particles (FLOPP and FLOPP-e) for the analysis of microplastics. Anal Chem. 2021;93(48):15878–15885. doi:10.1021/acs.analchem.1c02549

31. De Frond HL, Rubinovitz R, Rochman CM. μATR-FTIR spectral libraries of plastic particles (FLOPP and FLOPP-e) for the analysis of microplastics. Anal Chem. 2021;93(48):15878–15885. doi:10.1021/acs.analchem.1c02549

32. Tsuge S, Ohtani H, Watanabe C. Pyrolysis-GC/MS Data Book of Synthetic Polymers. Amsterdam: Elsevier; 2011. doi:10.1016/B978-0-444-53892-5.10001-9

33. Liu P, Zhan X, Wu X, Li J, Wang H, Gao S. Effect of weathering on environmental behavior of microplastics: properties, sorption and potential risks. Chemosphere. 2020;242:125193. doi:10.1016/j.chemosphere.2019.125193

34. Veerasingam S, Ranjani M, Venkatachalapathy R, Bagaev A, Mukhanov V, Liber D, et al. Contributions of Fourier transform infrared spectroscopy in microplastic pollution research: a review. Crit Rev Environ Sci Technol. 2021;51(22):2681–2743. doi:10.1080/10643389.2020.1807450

35. Munno K, De Frond H, O'Donnell B, Rochman CM. Increasing the accessibility for characterizing microplastics: introducing new application-based and spectral libraries of plastic particles (SLoPP and SLoPP-E). Anal Chem. 2020;92(3):2443–2451. doi:10.1021/acs.analchem.9b03626

36. Gigault J, El Hadri H, Nguyen B, Grassl B, Rowenczyk L, Tufenkji N, et al. Nanoplastics are neither microplastics nor engineered nanoparticles. Nat Nanotechnol. 2021;16(5):501–507. doi:10.1038/s41565-021-00886-4

37. McDonald RS, Wilks PA Jr. JCAMP-DX: a standard form for exchange of infrared spectra in computer readable form. Appl Spectrosc. 1988;42(1):151–162. doi:10.1366/0003702884428734

38. Koelmans AA, Redondo-Hasselerharm PE, Mohamed Nor NH, de Ruijter VN, Mintenig SM, Kooi M. Risk assessment of microplastic particles. Nat Rev Mater. 2022;7(2):138–152. doi:10.1038/s41578-021-00411-y

39. Hartmann NB, Hüffer T, Thompson RC, Hassellöv M, Verschoor A, Daugaard AE, et al. Are we speaking the same language? Recommendations for a definition and categorization framework for plastic debris. Environ Sci Technol. 2019;53(3):1039–1047. doi:10.1021/acs.est.8b05297

40. Wesch C, Elert AM, Wörner M, Stöckel S, Hoppe R. Assuring quality in microplastic monitoring: about the value of clean-air devices as essentials for verified data. Sci Rep. 2017;7(1):5424. doi:10.1038/s41598-017-05838-4

41. Dehaut A, Cassone AL, Frère L, Hermabessiere L, Himber C, Rinnert E, et al. Microplastics in seafood: benchmark protocol for their extraction and characterization. Environ Pollut. 2016;215:223–233. doi:10.1016/j.envpol.2016.05.018

42. Deng Y, Zhang Y, Lemos B, Ren H. Tissue accumulation of microplastics in mice and biomarker responses suggest widespread health risks of exposure. Sci Rep. 2017;7:46687. doi:10.1038/srep46687

43. Lu L, Wan Z, Luo T, Fu Z, Jin Y. Polystyrene microplastics induce gut microbiota dysbiosis and hepatic lipid metabolism disorder in mice. Sci Total Environ. 2018;631–632:449–458. doi:10.1016/j.scitotenv.2018.03.051

44. Stock V, Böhmert L, Lisber E, Block R, Cara-Carmona J, Pack LK, et al. Uptake and effects of orally ingested polystyrene microplastic particles in vitro and in vivo. Arch Toxicol. 2019;93(7):1817–1833. doi:10.1007/s00204-019-02478-7

45. Lu K, Lai KP, Stoeger T, Ji S, Lin Z, Lin X, et al. Detrimental effects of microplastic exposure on normal and asthmatic pulmonary physiology. J Hazard Mater. 2021;416:126069. doi:10.1016/j.jhazmat.2021.126069

46. Vianello A, Jensen RL, Liu L, Vollertsen J. Simulating human exposure to indoor airborne microplastics using a breathing thermal manikin. Sci Rep. 2019;9(1):8670. doi:10.1038/s41598-019-45054-w

47. Walczak AP, Hendriksen PJM, Woutersen RA, van der Zande M, Undas AK, Musber R, et al. Bioavailability and biodistribution of differently charged polystyrene nanoparticles upon oral exposure in rats. J Nanopart Res. 2015;17:231. doi:10.1007/s11051-015-3029-y

48. Suk JS, Xu Q, Kim N, Hanes J, Ensign LM. PEGylation as a strategy for improving nanoparticle-based drug and gene delivery. Adv Drug Deliv Rev. 2016;99(Pt A):28–51. doi:10.1016/j.addr.2015.09.012

49. Yee MS, Hii LW, Looi CK, Lim WM, Wong SF, Kok YY, et al. Impact of microplastics and nanoplastics on human health. Nanomaterials (Basel). 2021;11(2):496. doi:10.3390/nano11020496

50. Schwabl P, Köppel S, Königshofer P, Bucsics T, Trauner M, Reiberger T, et al. Detection of various microplastics in human stool: a prospective case series. Ann Intern Med. 2019;171(7):453–457. doi:10.7326/M19-0618

51. Choi HS, Liu W, Misra P, Tanaka E, Zimmer JP, Itty Ipe B, et al. Renal clearance of quantum dots. Nat Biotechnol. 2007;25(10):1165–1170. doi:10.1038/nbt1340

52. Kopatz V, Wen K, Kovács T, Keber AS, Gasber V, Vollnhofer N, et al. Micro- and nanoplastics breach the blood-brain barrier (BBB): biomolecular corona's role revealed. Nanomaterials (Basel). 2023;13(8):1404. doi:10.3390/nano13081404

53. Nihart AJ, Garcia MA, Liu R, Godoy Vitorino F, Garcia-Maldonado E, Orengo-Orengo L, et al. Bioaccumulation of microplastics in decedent human brains. Nat Med. 2025;31(4):1114–1119. doi:10.1038/s41591-024-03453-1

54. Jin H, Ma T, Sha X, Liu Z, Zhou Y, Meng X, et al. Polystyrene microplastics induced male reproductive toxicity in mice. J Hazard Mater. 2021;401:123430. doi:10.1016/j.jhazmat.2020.123430

55. An R, Wang X, Yang L, Zhang J, Wang N, Xu F, et al. Polystyrene microplastics cause granulosa cells apoptosis and fibrosis in ovary through oxidative stress in rats. Toxicology. 2021;449:152665. doi:10.1016/j.tox.2020.152665

56. Wick P, Malek A, Manser P, Meili D, Maeder-Althaus X, Diber L, et al. Barrier capacity of human placenta for nanosized materials. Environ Health Perspect. 2010;118(3):432–436. doi:10.1289/ehp.0901200

57. Ragusa A, Svelato A, Santacroce C, Catalano P, Notarstefano V, Carnevali O, et al. Plasticenta: first evidence of microplastics in human placenta. Environ Int. 2021;146:106274. doi:10.1016/j.envint.2020.106274

58. Turner JR. Intestinal mucosal barrier function in health and disease. Nat Rev Immunol. 2009;9(11):799–809. doi:10.1038/nri2653

59. Jani P, Halbert GW, Langridge J, Florence AT. Nanoparticle uptake by the rat gastrointestinal mucosa: quantitation and particle size dependency. J Pharm Pharmacol. 1990;42(12):821–826. doi:10.1111/j.2042-7158.1990.tb07033.x

60. Walczak AP, Kramer E, Hendriksen PJM, Tromp P, Helsper JPFG, van der Zande M, et al. Translocation of differently sized and charged polystyrene nanoparticles in in vitro intestinal cell models of increasing complexity. Nanotoxicology. 2015;9(4):453–461. doi:10.3109/17435390.2014.944599

61. des Rieux A, Ragnarsson EGE, Gullberg E, Préat V, Schneider YJ, Artursson P. Transport of nanoparticles across an in vitro model of the human intestinal follicle associated epithelium. Eur J Pharm Sci. 2005;25(4–5):455–465. doi:10.1016/j.ejps.2005.04.015

62. Luo T, Wang C, Pan Z, Jin C, Fu Z, Jin Y. Maternal polystyrene microplastic exposure during gestation and lactation altered metabolic homeostasis in the dams and their F1 and F2 offspring. Environ Sci Technol. 2019;53(18):10978–10992. doi:10.1021/acs.est.9b03191

63. Hesler M, Aengenheister L, Ellinger B, Drexel R, Straskraba S, Jost C, et al. Multi-endpoint toxicological assessment of polystyrene nano- and microparticles in different biological models in vitro. Toxicol In Vitro. 2019;61:104610. doi:10.1016/j.tiv.2019.104610

64. Lehner R, Weder C, Petri-Fink A, Rothen-Rutishauser B. Emergence of nanoplastic in the environment and possible impact on human health. Environ Sci Technol. 2019;53(4):1748–1765. doi:10.1021/acs.est.8b05512

65. Stock V, Fahrenson C, Thuenemann A, Dönmez MH, Böhmert L, Liese S, et al. Impact of artificial digestion on the sizes and shapes of microplastic particles. Food Chem Toxicol. 2020;135:111010. doi:10.1016/j.fct.2019.111010

66. Lai SK, Wang YY, Hanes J. Mucus-penetrating nanoparticles for drug and gene delivery to mucosal tissues. Adv Drug Deliv Rev. 2009;61(2):158–171. doi:10.1016/j.addr.2008.11.002

67. Weibel ER. Morphometry of the Human Lung. Berlin: Springer; 1963. doi:10.1007/978-3-642-87553-3

68. International Commission on Radiological Protection (ICRP). Human respiratory tract model for radiological protection. ICRP Publication 66. Ann ICRP. 1994;24(1–3):1–482.

69. Möller W, Felten K, Sommerer K, Scheuch G, Meyer G, Meyer P, et al. Deposition, retention, and translocation of ultrafine particles from the central airways and lung periphery. Am J Respir Crit Care Med. 2008;177(4):426–432. doi:10.1164/rccm.200602-301OC

70. Kreyling WG, Hirn S, Möller W, Schleh C, Wenk A, Celik G, et al. Air-blood barrier translocation of tracheally instilled gold nanoparticles inversely depends on particle size. ACS Nano. 2014;8(1):222–233. doi:10.1021/nn403256v

71. Oberdörster G, Oberdörster E, Oberdörster J. Nanotoxicology: an emerging discipline evolving from studies of ultrafine particles. Environ Health Perspect. 2005;113(7):823–839. doi:10.1289/ehp.7339

72. Donaldson K, Murphy FA, Duffin R, Poland CA. Asbestos, carbon nanotubes and the pleural mesothelium: a review of the hypothesis regarding the role of long fibre retention in the parietal pleura, inflammation and mesothelioma. Part Fibre Toxicol. 2010;7:5. doi:10.1186/1743-8977-7-5

73. Abbott NJ, Patabendige AAK, Dolman DEM, Yusof SR, Begley DJ. Structure and function of the blood-brain barrier. Neurobiol Dis. 2010;37(1):13–25. doi:10.1016/j.nbd.2009.07.030

74. Shan S, Zhang Y, Zhao H, Zeng T, Zhao X. Polystyrene nanoplastics penetrate across the blood-brain barrier and induce activation of microglia in the brain of mice. Chemosphere. 2022;298:134261. doi:10.1016/j.chemosphere.2022.134261

75. Teleanu DM, Chircov C, Grumezescu AM, Volceanov A, Teleanu RI. Blood-brain delivery methods using nanotechnology. Pharmaceutics. 2018;10(4):269. doi:10.3390/pharmaceutics10040269

76. Garcia MA, Liu R, Nihart AJ, El Hayek E, Castillo E, Wen M, et al. Quantitation and identification of microplastics accumulation in human placental specimens using pyrolysis gas chromatography-mass spectrometry. Toxicol Sci. 2024;199(2):235–245. PMID:38366932. doi:10.1093/toxsci/kfae021

77. Grafmueller S, Manser P, Diber L, Diber A, Jochum W, Krug HF, et al. Bidirectional transfer study of polystyrene nanoparticles across the placental barrier in an ex vivo human placental perfusion model. Environ Health Perspect. 2015;123(12):1280–1286. doi:10.1289/ehp.1409271

78. Halfar J, Čabanová K, Vávra K, Delongová P, Motyka O, Špaček R, et al. Microplastics and additives in patients with preterm birth: the first evidence of their presence in both human amniotic fluid and placenta. Chemosphere. 2023;343:140301. doi:10.1016/j.chemosphere.2023.140301

79. Braun T, Ehrlich L, Henrich W, Koeppel S, Lomako I, Schwabl P, et al. Detection of microplastic in human placenta and meconium in a clinical setting. Pharmaceutics. 2021;13(7):921. doi:10.3390/pharmaceutics13070921

80. Liu S, Liu X, Guo J, Yang R, Wang H, Sun Y, et al. The association between microplastics and microbiota in placentas and meconium: the first evidence in humans. Environ Sci Technol. 2023;57(46):17774–17785. doi:10.1021/acs.est.2c04706

81. Rejman J, Oberle V, Zuhorn IS, Hoekstra D. Size-dependent internalization of particles via the pathways of clathrin- and caveolae-mediated endocytosis. Biochem J. 2004;377(Pt 1):159–169. doi:10.1042/BJ20031253

82. Rejman J, Oberle V, Zuhorn IS, Hoekstra D. Size-dependent internalization of particles via the pathways of clathrin- and caveolae-mediated endocytosis. Biochem J. 2004;377(Pt 1):159–169. doi:10.1042/BJ20031253

83. Parton RG, del Pozo MA. Caveolae as plasma membrane sensors, protectors and organizers. Nat Rev Mol Cell Biol. 2013;14(2):98–112. doi:10.1038/nrm3512

84. Lim JP, Gleeson PA. Macropinocytosis: an endocytic pathway for internalising large gulps. Immunol Cell Biol. 2011;89(8):836–843. doi:10.1038/icb.2011.20

85. Gustafson HH, Holt-Casper D, Grainger DW, Ghandehari H. Nanoparticle uptake: the phagocyte problem. Nano Today. 2015;10(4):487–510. doi:10.1016/j.nantod.2015.06.006

86. Monopoli MP, Åberg C, Salvati A, Dawson KA. Biomolecular coronas provide the biological identity of nanosized materials. Nat Nanotechnol. 2012;7(12):779–786. doi:10.1038/nnano.2012.207

87. Tenzer S, Docter D, Kuharev J, Musyanovych A, Fetz V, Hecht R, et al. Rapid formation of plasma protein corona critically affects nanoparticle pathophysiology. Nat Nanotechnol. 2013;8(10):772–781. doi:10.1038/nnano.2013.181

88. Shen MC, Zhang YX, Zhu Y, Song B, Zeng GM, Hu DF, et al. Recent advances in toxicological research of nanoplastics in the environment: a review. Environ Pollut. 2019;252(Pt A):511–521. doi:10.1016/j.envpol.2019.05.102

89. Swanson KV, Deng M, Ting JP. The NLRP3 inflammasome: molecular activation and regulation to therapeutics. Nat Rev Immunol. 2019;19(8):477–489. doi:10.1038/s41577-019-0165-0

90. Li B, Ding Y, Cheng X, Sheng D, Xu Z, Rong Q, et al. Polyethylene microplastics affect the distribution of gut microbiota and inflammation development in mice. Chemosphere. 2020;244:125492. doi:10.1016/j.chemosphere.2019.125492

91. Rossi G, Barnoud J, Monticelli L. Polystyrene nanoparticles perturb lipid membranes. J Phys Chem Lett. 2014;5(1):241–246. doi:10.1021/jz402234c

92. Mabbott NA, Donaldson DS, Ohno H, Williams IR, Mahajan A. Microfold (M) cells: important immunosurveillance posts in the intestinal epithelium. Mucosal Immunol. 2013;6(4):666–677. doi:10.1038/mi.2013.30

93. Brayden DJ, Jepson MA, Baird AW. Keynote review: intestinal Peyer's patch M cells and oral vaccine targeting. Drug Discov Today. 2005;10(17):1145–1157. doi:10.1016/S1359-6446(05)03536-1

94. Florence AT. Nanoparticle uptake by the oral route: fulfilling its promise? Drug Discov Today Technol. 2005;2(1):75–81. doi:10.1016/j.ddtec.2005.05.019

95. des Rieux A, Fievez V, Garinot M, Schneider YJ, Préat V. Nanoparticles as potential oral delivery systems of proteins and vaccines: a mechanistic approach. J Control Release. 2006;116(1):1–27. doi:10.1016/j.jconrel.2006.08.013

96. Prata JC, da Costa JP, Lopes I, Duarte AC, Rocha-Santos T. Environmental exposure to microplastics: an overview on possible human health effects. Sci Total Environ. 2020;702:134455. doi:10.1016/j.scitotenv.2019.134455

97. Muhle H, Bellmann B, Creutzenberg O, Dasenbrock C, Ernst H, Kilpper R, et al. Pulmonary response to toner upon chronic inhalation exposure in rats. Fundam Appl Toxicol. 1991;17(2):280–299. doi:10.1016/0272-0590(91)90219-T

98. Geyer R, Jambeck JR, Law KL. Production, use, and fate of all plastics ever made. Sci Adv. 2017;3(7):e1700782. doi:10.1126/sciadv.1700782

99. Mohamed Nor NH, Kooi M, Diepens NJ, Koelmans AA. Lifetime accumulation of microplastic in children and adults. Environ Sci Technol. 2021;55(8):5084–5096. doi:10.1021/acs.est.0c07384

100. Boyer JL. Bile formation and secretion. Compr Physiol. 2013;3(3):1035–1078. doi:10.1002/cphy.c120027

101. Longmire M, Choyke PL, Kobayashi H. Clearance properties of nano-sized particles and molecules as imaging agents: considerations and caveats. Nanomedicine (Lond). 2008;3(5):703–717. doi:10.2217/17435889.3.5.703

102. Arslan B, Islamoğlu Y, Berçin AS, Akbulut S, Melikoğlu M. Detection and quantification of microplastics in cerumen. Turk J Med Sci. 2025;55(3):e6043. doi:10.55730/1300-0144.6043

103. Saraluck A, Techarang T, Bunyapipat P, Boonchuwong K, Pullaput Y, Mordmuang A. Detection of microplastics in human breast milk and its association with changes in human milk bacterial microbiota. J Clin Med. 2024;13(14):4029. doi:10.3390/jcm13144029

104. Tokiwa Y, Calabia BP, Ugwu CU, Aiba S. Biodegradability of plastics. Int J Mol Sci. 2009;10(9):3722–3742. doi:10.3390/ijms10093722

105. Poland CA, Duffin R, Kinloch I, Maynard A, Wallace WA, Seaton A, et al. Carbon nanotubes introduced into the abdominal cavity of mice show asbestos-like pathogenicity in a pilot study. Nat Nanotechnol. 2008;3(7):423–428. doi:10.1038/nnano.2008.111

106. Carson MJ, Doose JM, Melchior B, Schmid CD, Ploix CC. CNS immune privilege: hiding in plain sight. Immunol Rev. 2006;213:48–65. doi:10.1111/j.1600-065X.2006.00441.x

107. Shaw AC, Goldstein DR, Montgomery RR. Age-dependent dysregulation of innate immunity. Nat Rev Immunol. 2013;13(12):875–887. doi:10.1038/nri3547

108. Jones HM, Rowland-Yeo K. Basic concepts in physiologically based pharmacokinetic modeling in drug discovery and development. CPT Pharmacometrics Syst Pharmacol. 2013;2(8):e63. doi:10.1038/psp.2013.41

109. Nestorov I. Whole body pharmacokinetic models. Clin Pharmacokinet. 2003;42(10):883–908. doi:10.2165/00003088-200342100-00002

110. Li M, Al-Jamal KT, Kostarelos K, Reineke J. Physiologically based pharmacokinetic modeling of nanoparticles. ACS Nano. 2010;4(11):6303–6317. doi:10.1021/nn1018818

111. Choi HS, Liu W, Liu F, Nasr K, Misra P, Bawendi MG, et al. Design considerations for tumour-targeted nanoparticles. Nat Nanotechnol. 2010;5(1):42–47. doi:10.1038/nnano.2009.314

112. Sadauskas E, Wallin H, Stoltenberg M, Vogel U, Doering P, Larsen A, et al. Kupffer cells are central in the removal of nanoparticles from the organism. Part Fibre Toxicol. 2007;4:10. doi:10.1186/1743-8977-4-10

113. Lin Z, Monteiro-Riviere NA, Riviere JE. A physiologically based pharmacokinetic model for polyethylene glycol-coated gold nanoparticles of different sizes in adult mice. Nanotoxicology. 2016;10(2):162–172. doi:10.3109/17435390.2015.1027314

114. Bachler G, von Goetz N, Hungerbühler K. A physiologically based pharmacokinetic model for ionic silver and silver nanoparticles. Int J Nanomedicine. 2013;8:3365–3382. doi:10.2147/IJN.S46624

115. Lankveld DPK, Oomen AG, Krystek P, Neigh A, Troost-de Jong A, Noorlander CW, et al. The kinetics of the tissue distribution of silver nanoparticles of different sizes. Biomaterials. 2010;31(32):8350–8361. doi:10.1016/j.biomaterials.2010.07.045

116. Bachler G, von Goetz N, Hungerbühler K. Using physiologically based pharmacokinetic (PBPK) modeling for dietary risk assessment of titanium dioxide (TiO₂) nanoparticles. Nanotoxicology. 2015;9(3):373–380. doi:10.3109/17435390.2014.940404

117. Meta-analysis data generated by systematic review pipeline. Search conducted January 15, 2026, across PubMed, Europe PMC, and OpenAlex. Random-effects meta-analysis (DerSimonian-Laird method) of 11 studies, N = 474. [Original analysis for this chapter.]

118. DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7(3):177–188. doi:10.1016/0197-2456(86)90046-2

119. Freeman MF, Tukey JW. Transformations related to the angular and the square root. Ann Math Stat. 1950;21(4):607–611. doi:10.1214/aoms/1177729756

120. Leonard SVL, Liddle CR, Atherall CA, Chapman E, Watkins M, Calaminus SDJ, et al. Microplastics in human blood: polymer types, concentrations and characterisation using μFTIR. Environ Int. 2024;189:108751. PMID:38761430. doi:10.1016/j.envint.2024.108751

121. Leslie HA, van Velzen MJM, Brandsma SH, Kroesbergen J, Allen JG, Lamoree MH. Discovery and quantification of plastic particle pollution in human blood. Environ Int. 2022;163:107199. PMID:35367073. doi:10.1016/j.envint.2022.107199

122. Padarya SK, Asati AA, Singh A, Saad T, Mishra S, Tripathi L. Bronchoalveolar lavage as a diagnostic window into human exposure to microplastics and associated lung changes. J Pharm Bioallied Sci. 2025. PMID:41523016. doi:10.4103/jpbs.jpbs_1095_25

123. Özgen Alpaydin A, Uçan ES, Köktürk M, Atamanalp M, Kalyoncu Ç, Yiğit S, et al. Microplastics, as a risk factor in the development of interstitial lung disease. Environ Pollut. 2024;361:125054. PMID:39368622. doi:10.1016/j.envpol.2024.125054

124. Zakynthinos GE, Dimeas IE, Salmas C, Pagonis AD, Papanikolaou ID, Oikonomou E, et al. Detection of microplastics in human bronchoalveolar lavage fluid: preliminary evidence. Cureus. 2025;17(11):e97632. PMID:41450414. doi:10.7759/cureus.97632

125. Horvatits T, Tamminga M, Liu B, Schwarz M, Tzavaras S, Tiesmeyer M, et al. Microplastics detected in cirrhotic liver tissue. EBioMedicine. 2022;82:104147. doi:10.1016/j.ebiom.2022.104147

126. Ke D, Zheng J, Liu X, Xu X, Zhao L, Gu Y, et al. Occurrence of microplastics and disturbance of gut microbiota: a pilot study of preschool children in Xiamen, China. EBioMedicine. 2023;97:104828. PMID:37837933. doi:10.1016/j.ebiom.2023.104828

127. Santini S, Exposito N, Sierra J, Cincinelli A, Rovira J. Set up and validation of a method to analyse microplastics in stool and small intestine samples. MethodsX. 2024;12:102777. PMID:38883588. doi:10.1016/j.mex.2024.102777

128. Demirelli E, Tepe Y, Oğuz U, Aydın H, Kodat M, Tok DS, et al. The first reported values of microplastics in prostate. BMC Urol. 2024;24(1):106. PMID:38745203. doi:10.1186/s12894-024-01495-8

129. Qu J, Zeng J, Mou L, Wu X, Ha M, Liu C. Plastic tableware use, microplastic accumulation, and sperm quality: from epidemiological evidence to FOXA1/p38 mechanistic insights. J Nanobiotechnology. 2025;23:287. PMID:41068921. doi:10.1186/s12951-025-03747-7

130. Nair AB, Jacob S. A simple practice guide for dose conversion between animals and human. J Basic Clin Pharm. 2016;7(2):27–31. doi:10.4103/0976-0105.177703

131. Lusher AL, Hollman PCH, Mendoza-Hill JJ. Microplastics in fisheries and aquaculture: status of knowledge on their occurrence and implications for aquatic organisms and food safety. FAO Fisheries and Aquaculture Technical Paper. No. 615. Rome: FAO; 2017.

132. Cox KD, Covernton GA, Davies HL, Dower JF, Juanes F, Dudas SE. Human consumption of microplastics. Environ Sci Technol. 2019;53(12):7068–7074. PMID:31184127. doi:10.1021/acs.est.9b01517

133. Riviere JE. Pharmacokinetics of nanomaterials: an overview of carbon nanotubes, fullerenes and quantum dots. Wiley Interdiscip Rev Nanomed Nanobiotechnol. 2009;1(1):26–34. doi:10.1002/wnan.24

134. Dris R, Gasperi J, Mirande C, Mandin C, Guerrouache M, Langlois V, et al. A first overview of textile fibers, including microplastics, in indoor and outdoor environments. Environ Pollut. 2017;221:453–458. doi:10.1016/j.envpol.2016.12.013

135. Schneider M, Stracke F, Hansen S, Schaefer UF. Nanoparticles and their interactions with the dermal barrier. Dermatoendocrinol. 2009;1(4):197–206. doi:10.4161/derm.1.4.9501

136. Fournier SB, D'Errico JN, Adler DS, Kollontzi S, Goedken MJ, Fabris L, et al. Nanopolystyrene translocation and fetal deposition after acute lung exposure during late-stage pregnancy. Part Fibre Toxicol. 2020;17(1):55. doi:10.1186/s12989-020-00385-9

---

## Figure Legends

**Figure 8.1.** Schematic illustration of iMNP translocation mechanisms across the intestinal epithelial barrier. Three primary pathways are depicted: (A) Transcellular transport via clathrin-mediated, caveolae-mediated, and macropinocytic endocytosis in enterocytes; (B) Paracellular transport through tight junction disruption induced by inflammatory signaling and oxidative stress; (C) M-cell mediated transcytosis in Peyer's patch follicle-associated epithelium, with basolateral delivery to subepithelial dendritic cells and entry into lymphatic circulation.

**Figure 8.2.** Schematic of a physiologically-based pharmacokinetic (PBPK) model for iMNP biodistribution. Organ compartments (lungs, liver, spleen, kidneys, gut, brain, placenta) are connected by arterial and venous blood flows. Key model parameters include: size-dependent barrier-crossing rate constants, reticuloendothelial system (RES) saturable uptake kinetics, and tissue-specific partition/sequestration coefficients. Clearance pathways (fecal, biliary, renal, mucociliary) are indicated.

**Figure 8.3.** Forest plot of random-effects meta-analysis showing detection prevalence of iMNPs across human tissues. Individual study estimates are displayed with 95% confidence intervals (blue squares; square size proportional to study weight). Tissue-specific subgroup pooled estimates are shown as red diamonds. The overall pooled prevalence was 68.5% (95% CI: 56.3–78.6%), with an I² of 77%. [Adapted from systematic review meta-analysis, see Section 8.9.]

**Figure 8.4.** Funnel plot for assessment of publication bias in the meta-analysis. Individual studies are plotted as logit-transformed prevalence (x-axis) versus standard error (y-axis, inverted). The dashed vertical line represents the pooled estimate. Pseudo 95% confidence limits are shown as dashed lines. Egger's regression test indicated no significant asymmetry (intercept = 1.72, p = 0.909). [Adapted from systematic review meta-analysis, see Section 8.9.]

**Figure 8.5.** Anatomical diagram showing iMNP detection prevalence across human tissues based on meta-analysis results. Detection prevalence is indicated for each tissue type: blood/serum (90.0%), placenta (90.0%), stool/gut (85.3%), cerumen/ear (83.3%), lung (66.1%), semen/testis (55.2%), and breast milk (39.0%). Exposure routes (oral ingestion, inhalation, dermal contact) and systemic distribution pathways are illustrated. Overall pooled prevalence across all tissues: 68.5% (95% CI: 56.3–78.6%, N = 474).

---

*Chapter word count: approximately 10,500 words (excluding references and tables)*
*Figures: 5 (Figures 8.1–8.5)*
*Tables: 4 (Tables 8.1–8.4)*
