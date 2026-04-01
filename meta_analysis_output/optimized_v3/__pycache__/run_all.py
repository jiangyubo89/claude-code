"""
Master Runner — Optimized Meta-Analysis Pipeline v3
====================================================
Runs all 3 steps sequentially:
  1. Multi-database search
  2. Full-text data extraction
  3. Meta-analysis + figure generation
"""
import sys, time

print("=" * 70)
print("  OPTIMIZED META-ANALYSIS PIPELINE v3")
print("  Microplastics/Nanoplastics in Human Tissues")
print("  Multi-Database | Full-Text | PRISMA 2020 | Publication-Quality")
print("=" * 70)

start = time.time()

print("\n\n" + "▓" * 70)
print("  STEP 1/3: Multi-Database Literature Search")
print("▓" * 70 + "\n")
from importlib import import_module
step1 = import_module('01_multi_db_search')
step1.main()

print("\n\n" + "▓" * 70)
print("  STEP 2/3: Full-Text Data Extraction")
print("▓" * 70 + "\n")
step2 = import_module('02_fulltext_extract')
step2.main()

print("\n\n" + "▓" * 70)
print("  STEP 3/3: Meta-Analysis & Figure Generation")
print("▓" * 70 + "\n")
step3 = import_module('03_meta_analysis')
step3.main()

elapsed = time.time() - start
print(f"\n\n{'=' * 70}")
print(f"  ALL STEPS COMPLETE in {elapsed/60:.1f} minutes")
print(f"{'=' * 70}")
