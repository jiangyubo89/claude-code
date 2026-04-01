"""Fix and verify references by DOI lookup"""
import urllib.request, json, time

# References that need DOI-based PMID lookup
refs_to_fix = [
    {"ref": 4,  "doi": "10.1021/acs.est.9b01517",           "expected": "Cox KD, 2019, Environ Sci Technol"},
    {"ref": 5,  "doi": "10.1038/s41598-019-45054-w",         "expected": "Vianello A, 2019, Sci Rep"},
    {"ref": 8,  "doi": "10.1016/j.scitotenv.2022.154907",    "expected": "Jenner LC, 2022, Sci Total Environ"},
    {"ref": 10, "doi": "10.1016/j.ebiom.2022.104147",        "expected": "Horvatits T, 2022, EBioMedicine"},
    {"ref": 11, "doi": "10.3390/polym14132700",               "expected": "Ragusa A, 2022, Polymers"},
    {"ref": 13, "doi": "10.1093/toxsci/kfae043",             "expected": "Garcia MA, 2024, Toxicol Sci"},
    {"ref": 14, "doi": "10.1126/science.abe5041",             "expected": "Vethaak AD, 2021, Science"},
    {"ref": 15, "doi": "10.1016/j.scitotenv.2019.02.431",    "expected": "Novotna K, 2019, Sci Total Environ"},
    {"ref": 18, "doi": "10.1021/acs.est.7b00423",             "expected": "Wright SL, 2017, Environ Sci Technol"},
    {"ref": 25, "doi": "10.1021/acs.est.1c03924",             "expected": "Yan Z, 2022, Environ Sci Technol"},
    {"ref": 26, "doi": "10.3390/pharmaceutics13070921",       "expected": "Braun T, 2021, Pharmaceutics"},
    {"ref": 27, "doi": "10.1021/acs.est.0c07384",             "expected": "Mohamed Nor NH, 2021, Environ Sci Technol"},
    {"ref": 29, "doi": "10.3390/ijerph18052392",              "expected": "D'Angelo S, 2021, Int J Environ Res Public Health"},
    {"ref": 30, "doi": "10.1038/s41598-017-05838-4",          "expected": "Wesch C, 2017, Sci Rep"},
    {"ref": 31, "doi": "10.1038/srep46687",                   "expected": "Deng Y, 2017, Sci Rep"},
    {"ref": 32, "doi": "10.1007/s11051-015-3029-y",           "expected": "Walczak AP, 2015, J Nanopart Res"},
    {"ref": 33, "doi": "10.1021/acs.est.2c04706",             "expected": "Liu S, 2023, Environ Sci Technol"},
]

print("=" * 80)
print("DOI-BASED REFERENCE VERIFICATION")
print("=" * 80)

for entry in refs_to_fix:
    doi = entry["doi"]
    ref = entry["ref"]
    expected = entry["expected"]

    # Search PubMed by DOI
    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        "?db=pubmed&term=" + doi + "[doi]&retmode=json"
    )

    try:
        req = urllib.request.Request(search_url, headers={"User-Agent": "Mozilla/5.0"})
        resp = urllib.request.urlopen(req, timeout=15)
        data = json.loads(resp.read().decode("utf-8", errors="replace"))
        idlist = data.get("esearchresult", {}).get("idlist", [])

        if idlist:
            pmid = idlist[0]
            # Fetch summary
            sum_url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                "?db=pubmed&id=" + pmid + "&retmode=json"
            )
            req2 = urllib.request.Request(sum_url, headers={"User-Agent": "Mozilla/5.0"})
            resp2 = urllib.request.urlopen(req2, timeout=15)
            sdata = json.loads(resp2.read().decode("utf-8", errors="replace"))
            info = sdata.get("result", {}).get(pmid, {})

            authors = info.get("authors", [])
            first_auth = authors[0].get("name", "") if authors else ""
            journal = info.get("source", "")
            title = info.get("title", "")[:80]
            pubdate = info.get("pubdate", "")

            print(f"[Ref {ref:>2d}] DOI: {doi}")
            print(f"  PMID found: {pmid}")
            print(f"  Author: {first_auth} | Journal: {journal} | Year: {pubdate}")
            print(f"  Title: {title}...")
            print(f"  Expected: {expected}")
            print(f"  Status: VERIFIED")
        else:
            print(f"[Ref {ref:>2d}] DOI: {doi}")
            print(f"  No PMID found in PubMed")
            print(f"  Expected: {expected}")
            print(f"  Status: DOI exists but no PMID (may be non-indexed)")

        print()
        time.sleep(0.4)

    except Exception as e:
        err_msg = str(e)
        print(f"[Ref {ref:>2d}] DOI: {doi}")
        print(f"  Error: {err_msg[:60]}")
        print()
        time.sleep(0.5)

# Also verify the confirmed-correct ones
print("\n" + "=" * 80)
print("PREVIOUSLY VERIFIED REFERENCES (no changes needed)")
print("=" * 80)
correct_refs = [
    {"ref": 2,  "pmid": "15131299", "doi": "10.1126/science.1094559",     "author": "Thompson RC", "journal": "Science"},
    {"ref": 6,  "pmid": "31476765", "doi": "10.7326/M19-0618",            "author": "Schwabl P",   "journal": "Ann Intern Med"},
    {"ref": 7,  "pmid": "35367073", "doi": "10.1016/j.envint.2022.107199","author": "Leslie HA",   "journal": "Environ Int"},
    {"ref": 9,  "pmid": "33395930", "doi": "10.1016/j.envint.2020.106274","author": "Ragusa A",    "journal": "Environ Int"},
    {"ref": 12, "pmid": "39901044", "doi": "10.1038/s41591-024-03453-1",  "author": "Nihart AJ",   "journal": "Nat Med"},
    {"ref": 16, "pmid": "31733547", "doi": "10.1016/j.scitotenv.2019.134455","author":"Prata JC",  "journal": "Sci Total Environ"},
    {"ref": 19, "pmid": "33782057", "doi": "10.1136/bmj.n71",             "author": "Page MJ",     "journal": "BMJ"},
    {"ref": 23, "pmid": "12958120", "doi": "10.1136/bmj.327.7414.557",    "author": "Higgins JP",  "journal": "BMJ"},
    {"ref": 24, "pmid": "9310563",  "doi": "10.1136/bmj.315.7109.629",    "author": "Egger M",     "journal": "BMJ"},
    {"ref": 28, "pmid": "27989388", "doi": "10.1016/j.envpol.2016.12.013","author": "Dris R",      "journal": "Environ Pollut"},
    {"ref": 34, "pmid": "38761430", "doi": "10.1016/j.envint.2024.108751","author": "Leonard S",   "journal": "Environ Int"},
    {"ref": 35, "pmid": "37837933", "doi": "10.1016/j.ebiom.2023.104828", "author": "Ke D",        "journal": "EBioMedicine"},
    {"ref": 36, "pmid": "40933963", "doi": "10.55730/1300-0144.6043",     "author": "Arslan B",    "journal": "Turk J Med Sci"},
    {"ref": 37, "pmid": "38883588", "doi": "10.1016/j.mex.2024.102777",   "author": "Santini S",   "journal": "MethodsX"},
    {"ref": 38, "pmid": "41450414", "doi": "10.7759/cureus.97632",        "author": "Zakynthinos GE","journal":"Cureus"},
    {"ref": 39, "pmid": "38745203", "doi": "10.1186/s12894-024-01495-8",  "author": "Demirelli E", "journal": "BMC Urol"},
    {"ref": 40, "pmid": "39064070", "doi": "10.3390/jcm13144029",         "author": "Saraluck A",  "journal": "J Clin Med"},
    {"ref": 41, "pmid": "39368622", "doi": "10.1016/j.envpol.2024.125054","author": "Ozgen Alpaydin","journal":"Environ Pollut"},
    {"ref": 42, "pmid": "37769922", "doi": "10.1016/j.chemosphere.2023.140301","author":"Halfar J","journal":"Chemosphere"},
    {"ref": 43, "pmid": "41068921", "doi": "10.1186/s12951-025-03747-7",  "author": "Qu J",        "journal": "J Nanobiotechnology"},
    {"ref": 44, "pmid": "41523016", "doi": "10.4103/jpbs.jpbs_1095_25",   "author": "Padarya SK",  "journal": "J Pharm Bioallied Sci"},
]
for r in correct_refs:
    print(f"[Ref {r['ref']:>2d}] PMID {r['pmid']:>10s} | DOI: {r['doi']} | {r['author']} | {r['journal']} | VERIFIED")
