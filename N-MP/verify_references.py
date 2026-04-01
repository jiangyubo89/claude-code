"""Verify all references via PubMed NCBI E-utilities"""
import urllib.request, json, time

pmids_to_verify = {
    '15131299': {'ref':2, 'expected_doi':'10.1126/science.1094559'},
    '31108247': {'ref':4, 'expected_doi':'10.1021/acs.est.9b01517'},
    '31197215': {'ref':5, 'expected_doi':'10.1038/s41598-019-45054-w'},
    '31476765': {'ref':6, 'expected_doi':'10.7326/M19-0618'},
    '35367073': {'ref':7, 'expected_doi':'10.1016/j.envint.2022.107199'},
    '35395230': {'ref':8, 'expected_doi':'10.1016/j.scitotenv.2022.154907'},
    '33395930': {'ref':9, 'expected_doi':'10.1016/j.envint.2020.106274'},
    '35777307': {'ref':10, 'expected_doi':'10.1016/j.ebiom.2022.104147'},
    '35808793': {'ref':11, 'expected_doi':'10.3390/polym14132700'},
    '39901044': {'ref':12, 'expected_doi':'10.1038/s41591-024-03453-1'},
    '38650467': {'ref':13, 'expected_doi':'10.1093/toxsci/kfae043'},
    '33542124': {'ref':14, 'expected_doi':'10.1126/science.abe5041'},
    '31733547': {'ref':16, 'expected_doi':'10.1016/j.scitotenv.2019.134455'},
    '28548990': {'ref':18, 'expected_doi':'10.1021/acs.est.7b00423'},
    '33782057': {'ref':19, 'expected_doi':'10.1136/bmj.n71'},
    '12958120': {'ref':23, 'expected_doi':'10.1136/bmj.327.7414.557'},
    '9310563':  {'ref':24, 'expected_doi':'10.1136/bmj.315.7109.629'},
    '34931829': {'ref':25, 'expected_doi':'10.1021/acs.est.1c03924'},
    '33641332': {'ref':27, 'expected_doi':'10.1021/acs.est.0c07384'},
    '27989388': {'ref':28, 'expected_doi':'10.1016/j.envpol.2016.12.013'},
    '28710363': {'ref':30, 'expected_doi':'10.1038/s41598-017-05838-4'},
    '28426019': {'ref':31, 'expected_doi':'10.1038/srep46687'},
    '37963230': {'ref':33, 'expected_doi':'10.1021/acs.est.2c04706'},
    '38761430': {'ref':34, 'expected_doi':'10.1016/j.envint.2024.108751'},
    '37837933': {'ref':35, 'expected_doi':'10.1016/j.ebiom.2023.104828'},
    '40933963': {'ref':36, 'expected_doi':'10.55730/1300-0144.6043'},
    '38883588': {'ref':37, 'expected_doi':'10.1016/j.mex.2024.102777'},
    '41450414': {'ref':38, 'expected_doi':'10.7759/cureus.97632'},
    '38745203': {'ref':39, 'expected_doi':'10.1186/s12894-024-01495-8'},
    '39064070': {'ref':40, 'expected_doi':'10.3390/jcm13144029'},
    '39368622': {'ref':41, 'expected_doi':'10.1016/j.envpol.2024.125054'},
    '37769922': {'ref':42, 'expected_doi':'10.1016/j.chemosphere.2023.140301'},
    '41068921': {'ref':43, 'expected_doi':'10.1186/s12951-025-03747-7'},
    '41523016': {'ref':44, 'expected_doi':'10.4103/jpbs.jpbs_1095_25'},
}

all_pmids = list(pmids_to_verify.keys())
url = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    "?db=pubmed&id=" + ",".join(all_pmids) + "&retmode=json"
)

print("Fetching from NCBI E-utilities...")
req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
resp = urllib.request.urlopen(req, timeout=30)
raw = resp.read().decode("utf-8", errors="replace")
data = json.loads(raw)
result_data = data.get("result", {})

print("=" * 80)
print("REFERENCE VERIFICATION REPORT")
print("=" * 80)

verified = 0
failed = 0
doi_list = []

for pmid in sorted(all_pmids, key=lambda x: pmids_to_verify[x]["ref"]):
    info = result_data.get(pmid, {})
    expected = pmids_to_verify[pmid]
    ref_num = expected["ref"]
    exp_doi = expected["expected_doi"]

    doi_found = ""
    for aid in info.get("articleids", []):
        if aid.get("idtype") == "doi":
            doi_found = aid.get("value", "")
            break

    first_auth = ""
    authors = info.get("authors", [])
    if authors:
        first_auth = authors[0].get("name", "")

    title = info.get("title", "")[:90]
    journal = info.get("source", "")
    exists = bool(title)

    status = "VERIFIED" if exists else "NOT FOUND"
    if exists:
        verified += 1
    else:
        failed += 1

    final_doi = doi_found if doi_found else exp_doi
    doi_list.append((ref_num, pmid, final_doi, first_auth, journal))

    print(f"[Ref {ref_num:>2d}] PMID {pmid}: {status}")
    print(f"  Author: {first_auth} | Journal: {journal}")
    print(f"  DOI: {doi_found}")
    if doi_found and exp_doi.lower() not in doi_found.lower():
        print(f"  NOTE: Expected DOI was {exp_doi}")
    print()

print(f"TOTAL: {verified} verified, {failed} failed / {len(all_pmids)}")

# Write DOI compilation file
doi_list.sort(key=lambda x: x[0])
output_path = r"D:\桌面\claude\N-MP\manuscript\reference_DOIs.csv"
with open(output_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["Ref_Number", "PMID", "DOI", "First_Author", "Journal"])
    for row in doi_list:
        writer.writerow(row)
print(f"\nDOI list written to: {output_path}")

import csv
