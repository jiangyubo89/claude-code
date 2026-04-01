"""
Step 2: Full-Text Data Extraction
=================================
Extracts quantitative data from full text and abstracts:
- Prevalence / detection rate
- Sample size
- Tissue type
- Polymer type
- Detection method
- Concentration data
- Study design classification
- Quality assessment (Modified Newcastle-Ottawa Scale)
"""
import pandas as pd
import numpy as np
import re, json, os
from xml.etree import ElementTree as ET
from pathlib import Path
from collections import Counter

OUT = Path(r"D:\桌面\meta_analysis_output\optimized_v3")

# ══════════════════════════════════════════════════════════════════════════
# TISSUE CLASSIFICATION
# ══════════════════════════════════════════════════════════════════════════
TISSUE_PATTERNS = {
    'Blood/Serum/Plasma': [
        r'\b(?:blood|serum|plasma|whole\s+blood|peripheral\s+blood|venous\s+blood)\b',
        r'\b(?:red\s+blood\s+cell|erythrocyte|leukocyte|white\s+blood\s+cell)\b',
        r'\bblood\s+sample',
    ],
    'Lung': [
        r'\b(?:lung|pulmonary|bronch|alveol|bronchoalveolar|BALF|bronchial)\b',
        r'\blung\s+tissue',
        r'\brespiratory\s+(?:tract|system|tissue)',
    ],
    'Placenta': [
        r'\b(?:placenta|placental|decidua|meconium|amniotic\s+fluid|cord\s+blood|umbilical)\b',
        r'\bfetal\s+(?:membrane|tissue)',
    ],
    'Stool/Gut': [
        r'\b(?:stool|feces|fecal|faecal|gut|intestin|colon|colonic|rectal|gastrointestin|GI\s+tract)\b',
        r'\b(?:bowel|duoden|jejun|ileum|cecum|appendix)\b',
    ],
    'Liver': [
        r'\b(?:liver|hepat|hepatocyte|bile|biliary|gallbladder)\b',
    ],
    'Breast milk': [
        r'\b(?:breast\s+milk|human\s+milk|lactation|colostrum|breastmilk)\b',
        r'\bmaternal\s+milk',
    ],
    'Semen/Testis': [
        r'\b(?:semen|sperm|testis|testicular|seminal|ejaculat|spermatozoa)\b',
    ],
    'Brain/Nerve': [
        r'\b(?:brain|cerebr|cortex|hippocamp|cerebrospinal|CSF|neural|neuron|nervous\s+system)\b',
        r'\b(?:olfactory|temporal\s+lobe|frontal\s+lobe)\b',
    ],
    'Kidney': [
        r'\b(?:kidney|renal|nephro|urine|urinary|urothelial)\b',
    ],
    'Adipose/Fat': [
        r'\b(?:adipose|fat\s+tissue|subcutaneous\s+fat|visceral\s+fat|lipoma)\b',
    ],
    'Thyroid': [
        r'\b(?:thyroid|thyroidectomy|thyroidal)\b',
    ],
    'Heart/Artery': [
        r'\b(?:heart|cardiac|myocard|atrial|ventricul|arteri|arterial|carotid|aort|cardiovascular)\b',
        r'\b(?:atheroscler|thrombus|atheromatous|vascular\s+tissue)\b',
    ],
    'Cervix/Uterus': [
        r'\b(?:cervi[xc]|uterus|uterine|endometri|ovary|ovarian|fallopian)\b',
    ],
    'Skin': [
        r'\b(?:skin|dermal|cutaneous|epiderm|subcutaneous)\b',
    ],
    'Bone Marrow': [
        r'\b(?:bone\s+marrow|hematopoietic)\b',
    ],
    'Salivary/Oral': [
        r'\b(?:saliva|salivary|oral\s+cavity|buccal|gingival)\b',
    ],
}

# ══════════════════════════════════════════════════════════════════════════
# POLYMER CLASSIFICATION
# ══════════════════════════════════════════════════════════════════════════
POLYMER_PATTERNS = {
    'Polyethylene (PE)': r'\b(?:polyethylene|PE\b|HDPE|LDPE|LLDPE)',
    'Polypropylene (PP)': r'\b(?:polypropylene|PP\b)',
    'Polystyrene (PS)': r'\b(?:polystyrene|PS\b|EPS\b|expanded\s+polystyrene)',
    'PET': r'\b(?:polyethylene\s+terephthalate|PET\b|polyester)',
    'PVC': r'\b(?:polyvinyl\s+chloride|PVC\b|vinyl\s+chloride)',
    'Nylon/PA': r'\b(?:nylon|polyamide|PA6|PA66|PA\s*6)',
    'Polycarbonate (PC)': r'\b(?:polycarbonate|PC\b|bisphenol)',
    'PTFE/Teflon': r'\b(?:polytetrafluoroethylene|PTFE|Teflon)',
    'Acrylic/PMMA': r'\b(?:poly\(?methyl\s*methacrylate\)?|PMMA|acrylic)',
    'Polyurethane (PU)': r'\b(?:polyurethane|PU\b)',
    'Tire rubber': r'\b(?:tire|tyre|rubber|styrene.?butadiene|SBR)',
    'Cellulose (modified)': r'\b(?:rayon|viscose|cellophane|cellulose\s+acetate)',
}

# ══════════════════════════════════════════════════════════════════════════
# DETECTION METHOD CLASSIFICATION
# ══════════════════════════════════════════════════════════════════════════
METHOD_PATTERNS = {
    'μFTIR/FTIR': r'\b(?:FTIR|μFTIR|micro.?FTIR|Fourier.?transform|ATR.?FTIR|FT.?IR)\b',
    'Raman/μRaman': r'\b(?:Raman|μRaman|micro.?Raman|confocal\s+Raman)\b',
    'Pyrolysis-GC/MS': r'\b(?:Py.?GC.?MS|pyrolysis.?GC|thermal\s+desorption.?GC|TED.?GC)\b',
    'SEM/TEM/AFM': r'\b(?:SEM|TEM|AFM|scanning\s+electron|transmission\s+electron|atomic\s+force)\b',
    'LC-MS/GC-MS': r'\b(?:LC.?MS|GC.?MS|HPLC.?MS|liquid\s+chromatograph|gas\s+chromatograph)\b',
    'Fluorescence': r'\b(?:fluorescen|Nile\s+Red|fluorescent\s+(?:stain|dye|microscop))\b',
    'Flow Cytometry': r'\b(?:flow\s+cytom|FACS|fluorescence.?activated)\b',
    'ICP-MS/ICP-OES': r'\b(?:ICP.?MS|ICP.?OES|inductively\s+coupled)\b',
    'LDIR': r'\b(?:LDIR|laser\s+direct\s+infrared)\b',
    'TOF-SIMS': r'\b(?:TOF.?SIMS|time.?of.?flight\s+secondary\s+ion)\b',
}

# ══════════════════════════════════════════════════════════════════════════
# STUDY DESIGN CLASSIFICATION
# ══════════════════════════════════════════════════════════════════════════
STUDY_TYPE_PATTERNS = {
    'Cross-sectional': r'\bcross.?section',
    'Cohort': r'\bcohort\b',
    'Case-control': r'\bcase.?control',
    'Systematic review': r'\bsystematic\s+review',
    'Meta-analysis': r'\bmeta.?analy',
    'Randomized controlled trial': r'\brandomized?\s+(?:controlled?\s+)?trial|RCT\b',
    'In vitro': r'\bin\s+vitro\b',
    'In vivo (animal)': r'\b(?:mouse|mice|rat|animal\s+model|murine|zebrafish|in\s+vivo)\b',
    'Exposure assessment': r'\bexposure\s+(?:assessment|study|evaluation)',
    'Biomonitoring': r'\bbiomonitor',
    'Autopsy/Surgical': r'\b(?:autopsy|surgical|post.?mortem|biopsy|resect)',
    'Environmental': r'\b(?:environmental\s+monitor|ambient\s+air|drinking\s+water|food\s+contamin)',
}


# ══════════════════════════════════════════════════════════════════════════
# QUANTITATIVE DATA EXTRACTION
# ══════════════════════════════════════════════════════════════════════════
def extract_prevalence(text):
    """
    Extract detection prevalence from text.
    Returns list of dicts: [{value, numerator, denominator, tissue, context}]
    """
    if not text or len(text) < 20:
        return []

    text_lower = text.lower()
    results = []

    # Pattern 1: "X out of Y (Z%) samples"
    for m in re.finditer(
        r'(\d+)\s+(?:out\s+of|of)\s+(\d+)\s*(?:\((\d+(?:\.\d+)?)\s*%\))?\s*'
        r'(?:(?:human\s+)?(?:samples?|participants?|subjects?|patients?|donors?|individuals?|specimens?|cases?))',
        text_lower
    ):
        n, total = int(m.group(1)), int(m.group(2))
        if 0 < n <= total and total >= 3:
            pct = round(100 * n / total, 1)
            results.append({
                'prevalence_pct': pct,
                'numerator': n,
                'denominator': total,
                'context': text_lower[max(0, m.start()-50):m.end()+50],
                'pattern': 'n_of_N',
            })

    # Pattern 2: "detected in X% of [N] samples"
    for m in re.finditer(
        r'(?:detected?|found|identified|present|observed)\s+in\s+'
        r'(\d+(?:\.\d+)?)\s*%\s*'
        r'(?:of\s+(?:the\s+)?(?:(\d+)\s+)?'
        r'(?:(?:human\s+)?(?:samples?|participants?|subjects?|patients?|donors?|individuals?)))?',
        text_lower
    ):
        pct = float(m.group(1))
        n = int(m.group(2)) if m.group(2) else None
        if 0 < pct <= 100:
            results.append({
                'prevalence_pct': pct,
                'numerator': round(pct * n / 100) if n else None,
                'denominator': n,
                'context': text_lower[max(0, m.start()-50):m.end()+50],
                'pattern': 'pct_of_N',
            })

    # Pattern 3: "X% detection rate / prevalence / positive rate"
    for m in re.finditer(
        r'(\d+(?:\.\d+)?)\s*%\s*'
        r'(?:detection\s+rate|prevalence|positive\s+rate|positivity\s+rate|incidence)',
        text_lower
    ):
        pct = float(m.group(1))
        if 0 < pct <= 100:
            results.append({
                'prevalence_pct': pct,
                'numerator': None,
                'denominator': None,
                'context': text_lower[max(0, m.start()-50):m.end()+50],
                'pattern': 'pct_rate',
            })

    # Pattern 4: "X% of samples/blood/etc contained/showed/had"
    for m in re.finditer(
        r'(\d+(?:\.\d+)?)\s*%\s*of\s*(?:the\s+)?(?:(\d+)\s+)?'
        r'(?:(?:human\s+)?(?:samples?|blood|plasma|serum|stool|fecal|lung|placenta|'
        r'tissue|urine|semen|milk|brain|liver|kidney))'
        r'\s*(?:samples?\s+)?'
        r'(?:contained?|showed?|had|were\s+positive|tested\s+positive|'
        r'were\s+found|were\s+detected|were\s+identified)',
        text_lower
    ):
        pct = float(m.group(1))
        n = int(m.group(2)) if m.group(2) else None
        if 0 < pct <= 100:
            results.append({
                'prevalence_pct': pct,
                'numerator': round(pct * n / 100) if n else None,
                'denominator': n,
                'context': text_lower[max(0, m.start()-50):m.end()+50],
                'pattern': 'pct_of_samples',
            })

    # Pattern 5: "microplastics were detected/found/identified in X%"
    for m in re.finditer(
        r'(?:micro\s*plastics?|nano\s*plastics?|MPs?|NPs?|plastic\s+particles?)\s+'
        r'(?:were\s+)?(?:detected?|found|identified|observed|present)\s+'
        r'in\s+(\d+(?:\.\d+)?)\s*%',
        text_lower
    ):
        pct = float(m.group(1))
        if 0 < pct <= 100:
            results.append({
                'prevalence_pct': pct,
                'numerator': None,
                'denominator': None,
                'context': text_lower[max(0, m.start()-50):m.end()+50],
                'pattern': 'mp_detected_pct',
            })

    # Pattern 6: "all N samples" (= 100%)
    for m in re.finditer(
        r'(?:all|every(?:\s+single)?)\s+(\d+)\s+'
        r'(?:(?:human\s+)?(?:samples?|participants?|subjects?|patients?|specimens?))\s+'
        r'(?:contained?|showed?|had|were\s+positive|tested\s+positive)',
        text_lower
    ):
        n = int(m.group(1))
        if n >= 3:
            results.append({
                'prevalence_pct': 100.0,
                'numerator': n,
                'denominator': n,
                'context': text_lower[max(0, m.start()-50):m.end()+50],
                'pattern': 'all_N',
            })

    # Pattern 7: "none of the N samples" (= 0%)
    for m in re.finditer(
        r'(?:none|no)\s+(?:of\s+(?:the\s+)?)?(\d+)\s+'
        r'(?:samples?|participants?|subjects?|patients?)',
        text_lower
    ):
        n = int(m.group(1))
        if n >= 3:
            results.append({
                'prevalence_pct': 0.0,
                'numerator': 0,
                'denominator': n,
                'context': text_lower[max(0, m.start()-50):m.end()+50],
                'pattern': 'none_N',
            })

    return results


def extract_sample_size(text):
    """Extract sample size with validation."""
    if not text or len(text) < 20:
        return None

    text_lower = text.lower()
    candidates = []

    patterns = [
        # "n = 42 participants/samples"
        (r'n\s*=\s*(\d+)\s*(?:human|adult|healthy|male|female|participants?|subjects?|'
         r'patients?|donors?|volunteers?|individuals?|samples?|specimens?)', 'n_equals'),
        # "enrolled/recruited/included N subjects"
        (r'(?:enrolled?|recruited?|included?|analyzed?|examined?|collected\s+from)\s+'
         r'(\d+)\s*(?:human|adult|healthy|male|female)?\s*(?:participants?|subjects?|'
         r'patients?|donors?|volunteers?|individuals?)', 'verb_N'),
        # "a total of N samples"
        (r'(?:a\s+)?total\s+of\s+(\d+)\s+(?:samples?|subjects?|patients?|participants?|'
         r'donors?|individuals?|specimens?)', 'total_of'),
        # "N blood/stool/lung/etc samples"
        (r'(\d+)\s+(?:human\s+)?(?:blood|serum|plasma|stool|fecal|lung|placenta|'
         r'tissue|brain|liver|kidney|urine|semen|milk|breast\s+milk|saliva)\s+'
         r'samples?', 'tissue_N'),
        # "from N participants"
        (r'(?:from|of|with)\s+(\d+)\s+(?:healthy\s+)?(?:participants?|subjects?|'
         r'patients?|donors?|volunteers?|individuals?)', 'from_N'),
        # "sample size of N"
        (r'sample\s+size\s+(?:of\s+)?(?:n\s*=\s*)?(\d+)', 'sample_size'),
    ]

    for pat, label in patterns:
        for m in re.finditer(pat, text_lower):
            try:
                n = int(m.group(1))
                if 3 <= n <= 100000:
                    candidates.append((n, label, m.start()))
            except:
                pass

    if not candidates:
        return None

    # Prefer "n =" patterns, then "total of", then others
    priority = {'n_equals': 0, 'total_of': 1, 'sample_size': 1,
                'verb_N': 2, 'tissue_N': 2, 'from_N': 3}
    candidates.sort(key=lambda x: (priority.get(x[1], 9), x[2]))

    return candidates[0][0]


def extract_concentration(text):
    """Extract MP/NP concentration data."""
    if not text:
        return []

    text_lower = text.lower()
    results = []

    # Pattern: "X.X ± Y.Y particles/unit"
    for m in re.finditer(
        r'(\d+(?:\.\d+)?)\s*[±+/-]+\s*(\d+(?:\.\d+)?)\s*'
        r'(particles?|items?|pieces?|counts?|fibers?|fragments?)\s*'
        r'(?:per\s+|/\s*)'
        r'(g|mg|kg|ml|mL|L|l|μL|individual|person|sample)',
        text_lower
    ):
        results.append({
            'mean': float(m.group(1)),
            'sd': float(m.group(2)),
            'unit_numerator': m.group(3),
            'unit_denominator': m.group(4),
            'context': text_lower[max(0, m.start()-30):m.end()+30],
        })

    # Pattern: "median X (IQR: Y-Z) particles/unit"
    for m in re.finditer(
        r'median\s+(?:of\s+)?(\d+(?:\.\d+)?)\s*'
        r'(?:\((?:IQR|interquartile)?\s*:?\s*(\d+(?:\.\d+)?)\s*[-–]\s*(\d+(?:\.\d+)?)\))?',
        text_lower
    ):
        results.append({
            'median': float(m.group(1)),
            'iqr_lo': float(m.group(2)) if m.group(2) else None,
            'iqr_hi': float(m.group(3)) if m.group(3) else None,
            'context': text_lower[max(0, m.start()-30):m.end()+30],
        })

    return results


# ══════════════════════════════════════════════════════════════════════════
# FULL TEXT PARSING
# ══════════════════════════════════════════════════════════════════════════
def parse_fulltext_xml(xml_text):
    """Extract structured text sections from PMC XML."""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return ''

    sections = []

    # Body text
    for p in root.findall('.//body//p'):
        text = ''.join(p.itertext()).strip()
        if text:
            sections.append(text)

    # Table captions and cells
    for table in root.findall('.//table-wrap'):
        caption = table.find('.//caption')
        if caption is not None:
            sections.append(''.join(caption.itertext()).strip())
        for cell in table.findall('.//td'):
            text = ''.join(cell.itertext()).strip()
            if text:
                sections.append(text)

    # Figure captions
    for fig in root.findall('.//fig'):
        caption = fig.find('.//caption')
        if caption is not None:
            sections.append(''.join(caption.itertext()).strip())

    return ' '.join(sections)


# ══════════════════════════════════════════════════════════════════════════
# QUALITY ASSESSMENT (Modified Newcastle-Ottawa Scale for Cross-Sectional)
# ══════════════════════════════════════════════════════════════════════════
def assess_quality(text, has_fulltext=False):
    """
    Automated quality assessment based on Modified NOS for cross-sectional studies.
    Returns score (0-9) and individual criteria.
    Maximum: 9 stars (Selection: 4, Comparability: 2, Outcome: 3)
    """
    text_lower = text.lower() if text else ''
    score = {}

    # Selection (max 4 stars)
    # S1: Representativeness of sample
    s1 = 0
    if re.search(r'\b(?:representative|random(?:ly)?|stratified|population.?based|national|multi.?cent[re])\b', text_lower):
        s1 = 1
    score['S1_representativeness'] = s1

    # S2: Sample size adequate
    s2 = 0
    if re.search(r'\b(?:power\s+(?:analysis|calculation)|sample\s+size\s+(?:was\s+)?calculat|adequate\s+sample)\b', text_lower):
        s2 = 1
    elif re.search(r'(?:n\s*=\s*\d{2,}|total\s+of\s+\d{2,})', text_lower):
        s2 = 1  # >= 10 samples is considered adequate for MP studies
    score['S2_sample_size'] = s2

    # S3: Non-respondents described
    s3 = 0
    if re.search(r'\b(?:non.?respond|exclusion\s+criteria|excluded?\s+\d|dropout|attrition|loss\s+to\s+follow)\b', text_lower):
        s3 = 1
    score['S3_nonrespondents'] = s3

    # S4: Ascertainment of exposure/detection
    s4 = 0
    if re.search(r'\b(?:validated|standard(?:ized)?|calibrat|quality\s+(?:control|assurance)|reference\s+(?:material|standard)|SOP|ISO)\b', text_lower):
        s4 = 1
    score['S4_ascertainment'] = s4

    # Comparability (max 2 stars)
    # C1: Controls for confounders
    c1 = 0
    if re.search(r'\b(?:adjust|confound|covariate|multivariat|regression|stratif|control(?:led)?\s+for)\b', text_lower):
        c1 = 1
    score['C1_confounders'] = c1

    # C2: Contamination controls
    c2 = 0
    if re.search(r'\b(?:blank|negative\s+control|procedural\s+blank|field\s+blank|contamination\s+control|'
                 r'clean\s+room|laminar\s+flow|background\s+contamination)\b', text_lower):
        c2 = 1
    score['C2_contamination_ctrl'] = c2

    # Outcome (max 3 stars)
    # O1: Assessment method validated
    o1 = 0
    if re.search(r'\b(?:FTIR|Raman|Py.?GC.?MS|SEM|TEM|validated\s+method|gold\s+standard)\b', text_lower):
        o1 = 1
    score['O1_method_validated'] = o1

    # O2: Statistical analysis appropriate
    o2 = 0
    if re.search(r'\b(?:mean|median|standard\s+deviation|confidence\s+interval|p\s*[<>=]|statistical|'
                 r'ANOVA|t.?test|chi.?square|Mann.?Whitney|Kruskal|regression)\b', text_lower):
        o2 = 1
    score['O2_statistics'] = o2

    # O3: Response rate / recovery rate reported
    o3 = 0
    if re.search(r'\b(?:recovery\s+rate|spike\s+recovery|recov\w+\s+\d+\s*%|limit\s+of\s+detect|LOD|LOQ)\b', text_lower):
        o3 = 1
    score['O3_recovery'] = o3

    total = sum(score.values())
    quality_grade = 'High' if total >= 7 else ('Moderate' if total >= 4 else 'Low')

    return {
        'quality_total': total,
        'quality_grade': quality_grade,
        **score,
    }


# ══════════════════════════════════════════════════════════════════════════
# COUNTRY EXTRACTION
# ══════════════════════════════════════════════════════════════════════════
COUNTRY_KEYWORDS = {
    'China': ['china', 'chinese', 'beijing', 'shanghai', 'guangzhou', 'nanjing', 'wuhan',
              'hangzhou', 'shenzhen', 'chengdu', 'tianjin', 'jinan', 'xiamen',
              'dalian', 'qingdao', 'harbin', 'shenyang', 'chinese academy'],
    'USA': ['united states', 'usa', 'u.s.a', 'u.s.', 'american', 'new york',
            'california', 'boston', 'chicago', 'texas', 'michigan', 'harvard',
            'stanford', 'columbia university', 'johns hopkins', 'yale', 'mit',
            'national institutes of health', 'NIH'],
    'Italy': ['italy', 'italian', 'rome', 'milan', 'napl', 'turin', 'florence',
              'bologna', 'padova', 'catania', 'università'],
    'Germany': ['germany', 'german', 'berlin', 'munich', 'hamburg', 'frankfurt',
                'heidelberg', 'universität'],
    'UK': ['united kingdom', 'england', 'london', 'oxford', 'cambridge', 'edinburgh',
           'manchester', 'bristol', 'uk', 'imperial college', 'ucl'],
    'South Korea': ['south korea', 'korea', 'korean', 'seoul', 'busan', 'incheon'],
    'India': ['india', 'indian', 'delhi', 'mumbai', 'bangalore', 'chennai', 'hyderabad',
              'kolkata', 'pune'],
    'Spain': ['spain', 'spanish', 'madrid', 'barcelona', 'valencia', 'seville'],
    'France': ['france', 'french', 'paris', 'lyon', 'marseille', 'toulouse'],
    'Netherlands': ['netherlands', 'dutch', 'amsterdam', 'utrecht', 'leiden',
                    'wageningen', 'rotterdam'],
    'Iran': ['iran', 'iranian', 'tehran', 'isfahan', 'tabriz', 'shiraz'],
    'Brazil': ['brazil', 'brazilian', 'são paulo', 'rio de janeiro'],
    'Australia': ['australia', 'australian', 'sydney', 'melbourne', 'brisbane'],
    'Japan': ['japan', 'japanese', 'tokyo', 'osaka', 'kyoto', 'nagoya'],
    'Canada': ['canada', 'canadian', 'toronto', 'vancouver', 'montreal', 'ottawa'],
    'Turkey': ['turkey', 'turkish', 'istanbul', 'ankara', 'izmir'],
    'Austria': ['austria', 'austrian', 'vienna', 'graz', 'innsbruck', 'medizinische universität wien'],
    'Poland': ['poland', 'polish', 'warsaw', 'krakow', 'wroclaw'],
    'Portugal': ['portugal', 'portuguese', 'lisbon', 'porto'],
    'Denmark': ['denmark', 'danish', 'copenhagen', 'aarhus'],
    'Sweden': ['sweden', 'swedish', 'stockholm', 'gothenburg', 'lund'],
    'Switzerland': ['switzerland', 'swiss', 'zurich', 'zürich', 'bern', 'geneva', 'ETH'],
    'Belgium': ['belgium', 'belgian', 'brussels', 'ghent', 'antwerp'],
    'Mexico': ['mexico', 'mexican'],
    'Saudi Arabia': ['saudi arabia', 'saudi', 'riyadh', 'jeddah'],
    'Egypt': ['egypt', 'egyptian', 'cairo'],
    'Greece': ['greece', 'greek', 'athens'],
    'Czech Republic': ['czech', 'prague', 'brno'],
    'Norway': ['norway', 'norwegian', 'oslo', 'bergen'],
    'Finland': ['finland', 'finnish', 'helsinki'],
    'Thailand': ['thailand', 'thai', 'bangkok'],
    'Indonesia': ['indonesia', 'indonesian', 'jakarta'],
    'Malaysia': ['malaysia', 'malaysian', 'kuala lumpur'],
    'Singapore': ['singapore'],
    'Pakistan': ['pakistan', 'pakistani', 'islamabad', 'karachi', 'lahore'],
    'Vietnam': ['vietnam', 'vietnamese', 'hanoi'],
    'Nigeria': ['nigeria', 'nigerian', 'lagos'],
    'South Africa': ['south africa', 'johannesburg', 'cape town'],
    'Argentina': ['argentina', 'argentine', 'buenos aires'],
    'Chile': ['chile', 'chilean', 'santiago'],
    'Colombia': ['colombia', 'colombian', 'bogota'],
    'New Zealand': ['new zealand', 'auckland', 'wellington'],
    'Israel': ['israel', 'israeli', 'tel aviv', 'jerusalem'],
    'Taiwan': ['taiwan', 'taiwanese', 'taipei'],
    'Ireland': ['ireland', 'irish', 'dublin'],
    'Romania': ['romania', 'romanian', 'bucharest'],
    'Croatia': ['croatia', 'croatian', 'zagreb'],
    'Slovakia': ['slovakia', 'slovak', 'bratislava'],
    'Slovenia': ['slovenia', 'slovenian', 'ljubljana'],
    'Hungary': ['hungary', 'hungarian', 'budapest'],
    'Serbia': ['serbia', 'serbian', 'belgrade'],
}


def extract_country(affiliation):
    """Extract country from affiliation string."""
    if not affiliation:
        return 'Unknown'
    aff_lower = affiliation.lower()
    for country, keywords in COUNTRY_KEYWORDS.items():
        for kw in keywords:
            if kw in aff_lower:
                return country
    return 'Other'


# ══════════════════════════════════════════════════════════════════════════
# MAIN EXTRACTION PIPELINE
# ══════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 60)
    print("STEP 2: Full-Text Data Extraction")
    print("=" * 60)

    # Load raw data
    df = pd.read_csv(OUT / 'all_papers_raw.csv', encoding='utf-8-sig')
    print(f"Loaded {len(df)} records")

    # Load full text cache
    fulltext_cache = {}
    ft_path = OUT / 'fulltext_cache.json'
    if ft_path.exists():
        with open(ft_path, 'r', encoding='utf-8') as f:
            fulltext_cache = json.load(f)
        print(f"Full text cache: {len(fulltext_cache)} articles")

    # Fill missing fields
    for col in ['abstract', 'title', 'keywords', 'affiliation', 'journal', 'mesh',
                'pmc_id', 'doi', 'pmid', 'pub_type']:
        if col not in df.columns:
            df[col] = ''
        df[col] = df[col].fillna('')

    df['year'] = pd.to_numeric(df['year'], errors='coerce').astype('Int64')

    # ── Build analysis text ──────────────────────────────────────────
    print("\nExtracting full text content...")
    full_texts = []
    for _, row in df.iterrows():
        # Try full text first
        ft_key = row.get('pmc_id', '') or row.get('pmid', '') or row.get('doi', '')
        ft_xml = fulltext_cache.get(ft_key, '')
        if ft_xml:
            parsed = parse_fulltext_xml(ft_xml)
            full_texts.append(parsed)
        else:
            # Fall back to abstract
            full_texts.append(str(row.get('abstract', '')))

    df['analysis_text'] = full_texts
    df['has_fulltext'] = [bool(fulltext_cache.get(
        row.get('pmc_id', '') or row.get('pmid', '') or row.get('doi', ''), ''
    )) for _, row in df.iterrows()]
    df['text_combined'] = (df['title'].str.lower() + ' ' + df['analysis_text'].str.lower())

    # ── Country extraction ───────────────────────────────────────────
    print("Extracting countries...")
    df['country'] = df['affiliation'].apply(extract_country)

    # ── Tissue classification ────────────────────────────────────────
    print("Classifying tissues...")
    for tissue, patterns in TISSUE_PATTERNS.items():
        col_name = f'tissue_{tissue}'
        combined_pat = '|'.join(patterns)
        df[col_name] = df['text_combined'].apply(
            lambda t: 1 if re.search(combined_pat, str(t), re.IGNORECASE) else 0
        )
    tissue_cols_actual = [f'tissue_{t}' for t in TISSUE_PATTERNS.keys()]
    df['tissue_any'] = df[tissue_cols_actual].max(axis=1)

    # ── Polymer classification ───────────────────────────────────────
    print("Classifying polymers...")
    for poly, pat in POLYMER_PATTERNS.items():
        col_name = f'polymer_{poly}'
        df[col_name] = df['text_combined'].apply(
            lambda t: 1 if re.search(pat, str(t), re.IGNORECASE) else 0
        )

    # ── Method classification ────────────────────────────────────────
    print("Classifying detection methods...")
    for method, pat in METHOD_PATTERNS.items():
        col_name = f'method_{method}'
        df[col_name] = df['text_combined'].apply(
            lambda t: 1 if re.search(pat, str(t), re.IGNORECASE) else 0
        )

    # ── Study design classification ──────────────────────────────────
    print("Classifying study designs...")
    def classify_study_type(text):
        text_lower = str(text).lower()
        for stype, pat in STUDY_TYPE_PATTERNS.items():
            if re.search(pat, text_lower):
                return stype
        # Default classification
        if re.search(r'\bhuman\b.*\b(?:detect|found|identif|quantif)\b', text_lower):
            return 'Observational'
        if re.search(r'\b(?:review|overview|survey)\b', text_lower):
            return 'Review'
        return 'Other'

    df['study_type'] = df['text_combined'].apply(classify_study_type)

    # ── Prevalence extraction ────────────────────────────────────────
    print("Extracting prevalence data...")
    prev_results = []
    sample_sizes = []
    for _, row in df.iterrows():
        text = row['analysis_text']
        prev_list = extract_prevalence(text)
        if prev_list:
            # Take the first valid result (usually the main finding)
            best = prev_list[0]
            prev_results.append(best['prevalence_pct'])
            if best.get('denominator'):
                sample_sizes.append(best['denominator'])
            else:
                ss = extract_sample_size(text)
                sample_sizes.append(ss)
        else:
            prev_results.append(None)
            ss = extract_sample_size(text)
            sample_sizes.append(ss)

    df['prevalence_pct'] = prev_results
    df['sample_size'] = sample_sizes

    # ── Concentration extraction ─────────────────────────────────────
    print("Extracting concentration data...")
    conc_data = []
    for _, row in df.iterrows():
        concs = extract_concentration(row['analysis_text'])
        if concs:
            conc_data.append(json.dumps(concs[0]))
        else:
            conc_data.append('')
    df['concentration_data'] = conc_data

    # ── Quality assessment ───────────────────────────────────────────
    print("Running quality assessment...")
    quality_results = []
    for _, row in df.iterrows():
        qa = assess_quality(row['analysis_text'], row.get('has_fulltext', False))
        quality_results.append(qa)

    qa_df = pd.DataFrame(quality_results)
    df = pd.concat([df, qa_df], axis=1)

    # ── Primary tissue assignment ────────────────────────────────────
    def get_primary_tissue(row):
        for t in TISSUE_PATTERNS.keys():
            col = f'tissue_{t}'
            if col in row.index and row[col] == 1:
                return t
        return 'Other/Unspecified'

    df['primary_tissue'] = df.apply(get_primary_tissue, axis=1)

    # ── Screening: Include only original research on human samples ───
    print("\nApplying inclusion/exclusion criteria...")

    # --- Exclusion 1: Reviews ---
    # Catch reviews by study_type AND by keyword patterns in title/abstract
    review_patterns = (
        r'\b(?:review|overview|meta.?analy|systematic\s+review|narrative\s+review|'
        r'scoping\s+review|critical\s+review|mini.?review|state.?of.?the.?art\s+review|'
        r'comprehensive\s+review|literature\s+review|current\s+(?:status|knowledge|understanding)|'
        r'recent\s+advances|recent\s+progress)\b'
    )
    title_lower = df['title'].str.lower().fillna('')
    abstract_lower = df['abstract'].str.lower().fillna('')
    df['exclude_review'] = (
        df['study_type'].isin(['Systematic review', 'Meta-analysis', 'Review']) |
        title_lower.str.contains(review_patterns, regex=True, na=False) |
        # Catch "this review" / "in this review" / "we review" in abstract
        abstract_lower.str.contains(
            r'\b(?:this\s+review|we\s+review|this\s+paper\s+reviews?|'
            r'aim(?:s|ed)?\s+(?:to|of)\s+(?:this\s+)?review|'
            r'we\s+summarize|we\s+provide\s+(?:a|an)\s+(?:overview|summary)|'
            r'this\s+(?:article|paper|study)\s+(?:reviews?|summarizes?|provides?\s+an?\s+overview))\b',
            regex=True, na=False
        )
    ).astype(int)

    # --- Exclusion 2: Non-human subjects ---
    # Detect studies primarily about animals, fish, marine organisms, birds
    non_human_title_patterns = (
        r'\b(?:mice|mouse|rat[s]?\b|murine|zebrafish|daphnia|'
        r'seabass|fish(?:es)?|bivalve|mussel|oyster|clam|mollusk|mollusc|'
        r'shrimp|crab|lobster|invertebrate|insect|bee[s]?\b|ant[s]?\b|'
        r'bird[s]?\b|avian|nestling|tit[s]?\b|sparrow|pigeon|'
        r'sponge[s]?\b|coral|algae|algal|plankton|copepod|'
        r'dog[s]?\b|cat[s]?\b|pig[s]?\b|porcine|bovine|cattle|'
        r'worm[s]?\b|nematode|earthworm|'
        r'marine\s+(?:ecosystem|organism|species|environment)|'
        r'aquaculture|aquatic\s+(?:organism|species|ecosystem)|'
        r'terrestrial\s+vertebrate|wildlife|'
        r'soil\s+(?:organism|fauna|ecosystem)|sediment\s+(?:sample|core))\b'
    )
    df['exclude_nonhuman'] = (
        title_lower.str.contains(non_human_title_patterns, regex=True, na=False)
    ).astype(int)

    # --- Exclusion 3: Animal studies ---
    df['exclude_animal'] = (
        df['study_type'] == 'In vivo (animal)'
    ).astype(int)

    # --- Exclusion 4: In vitro studies (without human tissue data) ---
    # More aggressive: exclude if title says "in vitro" and no human tissue
    invitro_patterns = r'\b(?:in\s+vitro|cell\s+(?:line|culture)|HeLa|HepG2|A549|Caco|MCF|cytotoxicity\s+assay)\b'
    df['exclude_invitro'] = (
        (df['study_type'] == 'In vitro') |
        (title_lower.str.contains(invitro_patterns, regex=True, na=False) &
         ~abstract_lower.str.contains(r'\bhuman\s+(?:blood|tissue|sample|specimen|patient|subject|donor|volunteer)\b', regex=True, na=False))
    ).astype(int)

    # --- Exclusion 5: Environmental/food/water studies without human data ---
    env_food_patterns = (
        r'\b(?:seafood|drinking\s+water|tap\s+water|bottled\s+water|'
        r'wastewater|sewage|agricultural|farmland|'
        r'food\s+(?:product|packaging|contamination|safety)|'
        r'air\s+(?:quality|pollution|sample|particulate|filter)|'
        r'marine\s+pollution|ocean|beach|river|lake|estuary|'
        r'sediment\s+(?:contamination|pollution|analysis))\b'
    )
    df['exclude_env'] = (
        (df['study_type'] == 'Environmental') |
        (title_lower.str.contains(env_food_patterns, regex=True, na=False) &
         ~title_lower.str.contains(r'\bhuman\b', regex=True, na=False))
    ).astype(int)

    # --- Exclusion 6: Non-microplastics studies ---
    # Studies that mention MP keywords incidentally but are about other topics
    not_mp_patterns = (
        r'\b(?:bladder\s+cancer|intraoperative\s+cancer|tumor\s+(?:assessment|detection)|'
        r'harmful\s+algal\s+bloom|pesticide(?:s)?\s+exposure|herbicide|fungicide)\b'
    )
    df['exclude_not_mp'] = (
        title_lower.str.contains(not_mp_patterns, regex=True, na=False) &
        ~title_lower.str.contains(r'\b(?:microplastic|nanoplastic|plastic\s+particle)\b', regex=True, na=False)
    ).astype(int)

    # --- Exclusion 7: No abstract ---
    df['exclude_no_abstract'] = (df['abstract'].str.len() < 50).astype(int)

    df['excluded'] = df[['exclude_review', 'exclude_nonhuman', 'exclude_animal',
                         'exclude_invitro', 'exclude_env', 'exclude_not_mp',
                         'exclude_no_abstract']].max(axis=1)

    n_total = len(df)
    n_excluded = df['excluded'].sum()
    n_included = n_total - n_excluded

    print(f"  Total records: {n_total}")
    print(f"  Excluded total: {n_excluded}")
    print(f"    Reviews: {int(df['exclude_review'].sum())}")
    print(f"    Non-human subjects: {int(df['exclude_nonhuman'].sum())}")
    print(f"    Animal studies: {int(df['exclude_animal'].sum())}")
    print(f"    In vitro: {int(df['exclude_invitro'].sum())}")
    print(f"    Environmental/food: {int(df['exclude_env'].sum())}")
    print(f"    Not about MP/NP: {int(df['exclude_not_mp'].sum())}")
    print(f"    No abstract: {int(df['exclude_no_abstract'].sum())}")
    print(f"  Included: {n_included}")

    # ── Save PRISMA screening counts ─────────────────────────────────
    prisma_path = OUT / 'prisma_counts.json'
    if prisma_path.exists():
        with open(prisma_path, 'r') as f:
            prisma = json.load(f)
    else:
        prisma = {}

    prisma.update({
        'screening_total': n_total,
        'excluded_reviews': int(df['exclude_review'].sum()),
        'excluded_nonhuman': int(df['exclude_nonhuman'].sum()),
        'excluded_animal': int(df['exclude_animal'].sum()),
        'excluded_invitro': int(df['exclude_invitro'].sum()),
        'excluded_env': int(df['exclude_env'].sum()),
        'excluded_not_mp': int(df['exclude_not_mp'].sum()),
        'excluded_no_abstract': int(df['exclude_no_abstract'].sum()),
        'total_excluded': int(n_excluded),
        'included_studies': int(n_included),
        'with_prevalence': int(df.loc[df['excluded']==0, 'prevalence_pct'].notna().sum()),
        'with_sample_size': int(df.loc[df['excluded']==0, 'sample_size'].notna().sum()),
        'with_both': int(((df['excluded']==0) & df['prevalence_pct'].notna() & df['sample_size'].notna()).sum()),
        'with_fulltext': int(df.loc[df['excluded']==0, 'has_fulltext'].sum()),
    })

    with open(prisma_path, 'w') as f:
        json.dump(prisma, f, indent=2)

    # ── Summary stats ────────────────────────────────────────────────
    inc = df[df['excluded'] == 0]
    print(f"\n{'='*60}")
    print("EXTRACTION SUMMARY (Included studies only)")
    print(f"{'='*60}")
    print(f"Total included: {len(inc)}")
    print(f"With full text: {inc['has_fulltext'].sum()}")
    print(f"With prevalence data: {inc['prevalence_pct'].notna().sum()}")
    print(f"With sample size: {inc['sample_size'].notna().sum()}")
    print(f"With BOTH (for meta-analysis): {((inc['prevalence_pct'].notna()) & (inc['sample_size'].notna())).sum()}")
    print(f"\nTissue distribution:")
    for t in TISSUE_PATTERNS.keys():
        col = f'tissue_{t}'
        n = inc[col].sum()
        if n > 0:
            print(f"  {t:<25}: {n:4d} papers")
    print(f"\nStudy types:")
    print(inc['study_type'].value_counts().to_string())
    print(f"\nCountry distribution (top 15):")
    print(inc['country'].value_counts().head(15).to_string())
    print(f"\nQuality grades:")
    print(inc['quality_grade'].value_counts().to_string())

    # ── Save enhanced dataset ────────────────────────────────────────
    # Drop large columns for CSV
    save_cols = [c for c in df.columns if c not in ['analysis_text', 'text_combined', 'fulltext_xml']]
    df[save_cols].to_csv(OUT / 'papers_extracted.csv', index=False, encoding='utf-8-sig')
    print(f"\nSaved extracted data to {OUT / 'papers_extracted.csv'}")

    return df


if __name__ == '__main__':
    main()
