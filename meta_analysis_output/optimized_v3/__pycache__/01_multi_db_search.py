"""
Step 1: Multi-Database Literature Search
========================================
Searches PubMed, Europe PMC, OpenAlex, and Semantic Scholar
for microplastics/nanoplastics in human tissues.
Deduplicates by DOI/PMID and retrieves full text from PMC.
"""
import requests, json, time, os, re, csv, sys
import pandas as pd
from xml.etree import ElementTree as ET
from pathlib import Path

OUT = Path(r"D:\桌面\meta_analysis_output\optimized_v3")
OUT.mkdir(parents=True, exist_ok=True)

HEADERS = {'User-Agent': 'MetaAnalysis/1.0 (academic research)'}

# ── 1. PubMed Search ─────────────────────────────────────────────────────
def search_pubmed(query, max_results=5000):
    """Search PubMed and retrieve full metadata."""
    print("=" * 60)
    print("[1/4] Searching PubMed...")

    # Step 1: ESearch
    esearch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    params = {
        'db': 'pubmed',
        'term': query,
        'retmax': max_results,
        'retmode': 'json',
        'usehistory': 'y',
        'sort': 'relevance',
    }
    r = requests.get(esearch_url, params=params, timeout=60)
    data = r.json()
    ids = data['esearchresult']['idlist']
    webenv = data['esearchresult'].get('webenv', '')
    qkey = data['esearchresult'].get('querykey', '')
    total_count = int(data['esearchresult']['count'])
    print(f"  PubMed total hits: {total_count}, retrieving: {len(ids)}")

    # Step 2: EFetch in batches
    efetch_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    records = []
    BATCH = 200

    for start in range(0, len(ids), BATCH):
        batch_ids = ids[start:start + BATCH]
        params = {
            'db': 'pubmed',
            'id': ','.join(batch_ids),
            'rettype': 'xml',
            'retmode': 'xml',
        }
        for attempt in range(3):
            try:
                r = requests.get(efetch_url, params=params, timeout=90)
                if r.status_code == 200:
                    break
            except Exception as e:
                print(f"  Retry {attempt+1}: {e}")
                time.sleep(3)

        root = ET.fromstring(r.content)
        for article in root.findall('.//PubmedArticle'):
            rec = parse_pubmed_article(article)
            rec['source_db'] = 'PubMed'
            records.append(rec)

        print(f"  Fetched {min(start + BATCH, len(ids))}/{len(ids)}", end='\r')
        time.sleep(0.35)

    print(f"\n  PubMed records: {len(records)}")
    return records


def parse_pubmed_article(article):
    """Parse a single PubMed XML article element."""
    rec = {}

    # PMID
    pmid_el = article.find('.//PMID')
    rec['pmid'] = pmid_el.text if pmid_el is not None else ''

    # Title
    title_el = article.find('.//ArticleTitle')
    rec['title'] = ''.join(title_el.itertext()) if title_el is not None else ''

    # Abstract (structured)
    abs_parts = []
    for at in article.findall('.//AbstractText'):
        label = at.get('Label', '')
        text = ''.join(at.itertext())
        if label:
            abs_parts.append(f"{label}: {text}")
        else:
            abs_parts.append(text)
    rec['abstract'] = ' '.join(abs_parts)

    # Journal
    journal_el = article.find('.//Journal/Title')
    rec['journal'] = journal_el.text if journal_el is not None else ''
    jabr_el = article.find('.//Journal/ISOAbbreviation')
    rec['journal_abbr'] = jabr_el.text if jabr_el is not None else ''

    # Year
    pub_year = article.find('.//PubDate/Year')
    pub_med = article.find('.//PubDate/MedlineDate')
    if pub_year is not None:
        rec['year'] = pub_year.text
    elif pub_med is not None:
        rec['year'] = pub_med.text[:4]
    else:
        rec['year'] = ''

    # Authors
    authors = []
    for auth in article.findall('.//Author'):
        ln = auth.find('LastName')
        fn = auth.find('ForeName')
        if ln is not None:
            name = ln.text
            if fn is not None:
                name += ' ' + fn.text
            authors.append(name)
    rec['authors'] = '; '.join(authors)
    rec['first_author'] = authors[0] if authors else ''
    rec['n_authors'] = len(authors)

    # Affiliations (all)
    affs = []
    for aff_el in article.findall('.//AffiliationInfo/Affiliation'):
        if aff_el.text:
            affs.append(aff_el.text)
    rec['affiliation'] = ' | '.join(affs)

    # Keywords
    kws = [kw.text.strip() for kw in article.findall('.//Keyword') if kw.text]
    rec['keywords'] = '; '.join(kws)

    # MeSH
    mesh = [mh.text for mh in article.findall('.//MeshHeading/DescriptorName') if mh.text]
    rec['mesh'] = '; '.join(mesh)

    # DOI
    doi_el = article.find('.//ELocationID[@EIdType="doi"]')
    rec['doi'] = doi_el.text if doi_el is not None else ''
    # Also check ArticleId for DOI
    if not rec['doi']:
        for aid in article.findall('.//ArticleId'):
            if aid.get('IdType') == 'doi' and aid.text:
                rec['doi'] = aid.text
                break

    # PMC ID
    for aid in article.findall('.//ArticleId'):
        if aid.get('IdType') == 'pmc' and aid.text:
            rec['pmc_id'] = aid.text
            break
    else:
        rec['pmc_id'] = ''

    # Publication type
    pub_types = [pt.text for pt in article.findall('.//PublicationType') if pt.text]
    rec['pub_type'] = '; '.join(pub_types)

    return rec


# ── 2. Europe PMC Search ─────────────────────────────────────────────────
def search_europe_pmc(query, max_results=2000):
    """Search Europe PMC (includes PubMed, PMC, preprints)."""
    print("\n" + "=" * 60)
    print("[2/4] Searching Europe PMC...")

    records = []
    page_size = 1000
    cursor = '*'

    for page in range(0, max_results, page_size):
        url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
        params = {
            'query': query,
            'format': 'json',
            'pageSize': min(page_size, max_results - page),
            'cursorMark': cursor,
            'resultType': 'core',  # Full metadata
        }
        try:
            r = requests.get(url, params=params, headers=HEADERS, timeout=60)
            data = r.json()
        except Exception as e:
            print(f"  Europe PMC error: {e}")
            break

        results = data.get('resultList', {}).get('result', [])
        if not results:
            break

        for item in results:
            rec = {
                'pmid': str(item.get('pmid', '')),
                'pmc_id': item.get('pmcid', ''),
                'doi': item.get('doi', ''),
                'title': item.get('title', ''),
                'abstract': item.get('abstractText', ''),
                'journal': item.get('journalTitle', ''),
                'journal_abbr': item.get('journalIssn', ''),
                'year': str(item.get('pubYear', '')),
                'authors': '; '.join(
                    [f"{a.get('lastName', '')} {a.get('firstName', '')}"
                     for a in item.get('authorList', {}).get('author', [])]
                ),
                'first_author': '',
                'n_authors': len(item.get('authorList', {}).get('author', [])),
                'affiliation': item.get('affiliation', ''),
                'keywords': '; '.join(item.get('keywordList', {}).get('keyword', [])),
                'mesh': '',
                'pub_type': item.get('pubType', ''),
                'source_db': 'Europe_PMC',
                'cited_by': item.get('citedByCount', 0),
                'has_fulltext': item.get('isOpenAccess', 'N') == 'Y',
            }
            auth_list = item.get('authorList', {}).get('author', [])
            if auth_list:
                a0 = auth_list[0]
                rec['first_author'] = f"{a0.get('lastName', '')} {a0.get('firstName', '')}"
            records.append(rec)

        cursor = data.get('nextCursorMark', '')
        if not cursor or cursor == '*':
            break

        print(f"  Fetched {len(records)} records...", end='\r')
        time.sleep(0.5)

    print(f"\n  Europe PMC records: {len(records)}")
    return records


# ── 3. OpenAlex Search ───────────────────────────────────────────────────
def search_openalex(query_terms, max_results=2000):
    """Search OpenAlex for broader coverage."""
    print("\n" + "=" * 60)
    print("[3/4] Searching OpenAlex...")

    records = []
    page_size = 200
    base_url = 'https://api.openalex.org/works'

    # OpenAlex uses a different query format
    search_query = query_terms

    for page in range(1, (max_results // page_size) + 2):
        params = {
            'search': search_query,
            'per_page': page_size,
            'page': page,
            'filter': 'type:article,language:en',
            'select': 'id,doi,title,publication_year,authorships,primary_location,abstract_inverted_index,cited_by_count,concepts,open_access,type',
        }
        try:
            r = requests.get(base_url, params=params, headers={
                **HEADERS,
                'Accept': 'application/json',
            }, timeout=60)
            if r.status_code != 200:
                print(f"  OpenAlex HTTP {r.status_code}")
                break
            data = r.json()
        except Exception as e:
            print(f"  OpenAlex error: {e}")
            break

        results = data.get('results', [])
        if not results:
            break

        for item in results:
            # Reconstruct abstract from inverted index
            abstract = ''
            inv_idx = item.get('abstract_inverted_index')
            if inv_idx:
                positions = {}
                for word, pos_list in inv_idx.items():
                    for pos in pos_list:
                        positions[pos] = word
                abstract = ' '.join(positions[k] for k in sorted(positions.keys()))

            # Extract authors
            authorships = item.get('authorships', [])
            authors = []
            first_aff = ''
            for aship in authorships:
                auth = aship.get('author', {})
                name = auth.get('display_name', '')
                if name:
                    authors.append(name)
                if not first_aff and aship.get('institutions'):
                    inst = aship['institutions'][0]
                    first_aff = f"{inst.get('display_name', '')}, {inst.get('country_code', '')}"

            doi_raw = item.get('doi', '') or ''
            doi = doi_raw.replace('https://doi.org/', '') if doi_raw else ''

            rec = {
                'pmid': '',
                'pmc_id': '',
                'doi': doi,
                'title': item.get('title', '') or '',
                'abstract': abstract,
                'journal': '',
                'journal_abbr': '',
                'year': str(item.get('publication_year', '')),
                'authors': '; '.join(authors),
                'first_author': authors[0] if authors else '',
                'n_authors': len(authors),
                'affiliation': first_aff,
                'keywords': '',
                'mesh': '',
                'pub_type': item.get('type', ''),
                'source_db': 'OpenAlex',
                'cited_by': item.get('cited_by_count', 0),
            }

            # Journal from primary_location
            loc = item.get('primary_location', {}) or {}
            source = loc.get('source', {}) or {}
            rec['journal'] = source.get('display_name', '')

            records.append(rec)

        if len(records) >= max_results:
            break

        print(f"  Fetched {len(records)} records...", end='\r')
        time.sleep(0.2)

    print(f"\n  OpenAlex records: {len(records)}")
    return records


# ── 4. Semantic Scholar Search ───────────────────────────────────────────
def search_semantic_scholar(query, max_results=1000):
    """Search Semantic Scholar for citation-enriched data."""
    print("\n" + "=" * 60)
    print("[4/4] Searching Semantic Scholar...")

    records = []
    base_url = 'https://api.semanticscholar.org/graph/v1/paper/search'
    page_size = 100

    for offset in range(0, max_results, page_size):
        params = {
            'query': query,
            'offset': offset,
            'limit': page_size,
            'fields': 'title,abstract,year,authors,externalIds,journal,citationCount,publicationTypes,openAccessPdf',
        }
        try:
            r = requests.get(base_url, params=params, headers=HEADERS, timeout=60)
            if r.status_code == 429:
                print("  Rate limited, waiting 30s...")
                time.sleep(30)
                r = requests.get(base_url, params=params, headers=HEADERS, timeout=60)
            if r.status_code != 200:
                print(f"  S2 HTTP {r.status_code}")
                break
            data = r.json()
        except Exception as e:
            print(f"  Semantic Scholar error: {e}")
            break

        results = data.get('data', [])
        if not results:
            break

        for item in results:
            ext_ids = item.get('externalIds', {}) or {}
            authors_list = item.get('authors', []) or []
            journal_info = item.get('journal', {}) or {}

            rec = {
                'pmid': str(ext_ids.get('PubMed', '') or ''),
                'pmc_id': str(ext_ids.get('PubMedCentral', '') or ''),
                'doi': ext_ids.get('DOI', '') or '',
                'title': item.get('title', '') or '',
                'abstract': item.get('abstract', '') or '',
                'journal': journal_info.get('name', '') or '',
                'journal_abbr': '',
                'year': str(item.get('year', '') or ''),
                'authors': '; '.join([a.get('name', '') for a in authors_list]),
                'first_author': authors_list[0].get('name', '') if authors_list else '',
                'n_authors': len(authors_list),
                'affiliation': '',
                'keywords': '',
                'mesh': '',
                'pub_type': '; '.join(item.get('publicationTypes', []) or []),
                'source_db': 'SemanticScholar',
                'cited_by': item.get('citationCount', 0),
            }
            records.append(rec)

        if data.get('next') is None and offset + page_size >= data.get('total', 0):
            break

        print(f"  Fetched {len(records)} records...", end='\r')
        time.sleep(1.0)  # S2 rate limit is strict

    print(f"\n  Semantic Scholar records: {len(records)}")
    return records


# ── 5. Deduplication ─────────────────────────────────────────────────────
def deduplicate(all_records):
    """Deduplicate by DOI, then by PMID, then by title similarity."""
    print("\n" + "=" * 60)
    print("Deduplicating records...")

    seen_doi = {}
    seen_pmid = {}
    seen_title = {}
    unique = []
    dup_count = 0

    # Priority: PubMed > Europe PMC > OpenAlex > Semantic Scholar
    priority = {'PubMed': 0, 'Europe_PMC': 1, 'OpenAlex': 2, 'SemanticScholar': 3}
    all_records.sort(key=lambda x: priority.get(x.get('source_db', ''), 9))

    for rec in all_records:
        doi = rec.get('doi', '').strip().lower()
        pmid = rec.get('pmid', '').strip()
        title_norm = re.sub(r'[^a-z0-9]', '', rec.get('title', '').lower())

        is_dup = False

        # Check DOI
        if doi and len(doi) > 5:
            if doi in seen_doi:
                is_dup = True
                # Merge citation counts
                idx = seen_doi[doi]
                if rec.get('cited_by', 0) > unique[idx].get('cited_by', 0):
                    unique[idx]['cited_by'] = rec['cited_by']
                # Fill in missing fields
                for field in ['pmid', 'pmc_id', 'abstract', 'keywords', 'mesh', 'affiliation']:
                    if not unique[idx].get(field) and rec.get(field):
                        unique[idx][field] = rec[field]
            else:
                seen_doi[doi] = len(unique)

        # Check PMID
        if not is_dup and pmid and pmid != 'None' and pmid != '':
            if pmid in seen_pmid:
                is_dup = True
                idx = seen_pmid[pmid]
                if rec.get('cited_by', 0) > unique[idx].get('cited_by', 0):
                    unique[idx]['cited_by'] = rec['cited_by']
            else:
                seen_pmid[pmid] = len(unique)

        # Check title (fuzzy)
        if not is_dup and title_norm and len(title_norm) > 30:
            if title_norm in seen_title:
                is_dup = True
            else:
                seen_title[title_norm] = len(unique)

        if is_dup:
            dup_count += 1
        else:
            # Track which databases found this record
            rec['found_in_dbs'] = rec.get('source_db', '')
            unique.append(rec)

    print(f"  Total input: {len(all_records)}")
    print(f"  Duplicates removed: {dup_count}")
    print(f"  Unique records: {len(unique)}")

    return unique


# ── 6. Fetch Full Text from Europe PMC ───────────────────────────────────
def fetch_fulltext_pmc(records, max_fetch=500):
    """Fetch full text XML from Europe PMC for open-access articles."""
    print("\n" + "=" * 60)
    print("Fetching full text from Europe PMC...")

    pmc_records = [r for r in records if r.get('pmc_id') and r['pmc_id'].startswith('PMC')]
    print(f"  Records with PMC ID: {len(pmc_records)}")

    fetched = 0
    for rec in pmc_records[:max_fetch]:
        pmc_id = rec['pmc_id']
        url = f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc_id}/fullTextXML'
        try:
            r = requests.get(url, headers=HEADERS, timeout=30)
            if r.status_code == 200:
                rec['fulltext_xml'] = r.text
                rec['has_fulltext'] = True
                fetched += 1
            else:
                rec['fulltext_xml'] = ''
                rec['has_fulltext'] = False
        except:
            rec['fulltext_xml'] = ''
            rec['has_fulltext'] = False

        if fetched % 20 == 0:
            print(f"  Full text fetched: {fetched}/{len(pmc_records[:max_fetch])}", end='\r')
        time.sleep(0.3)

    print(f"\n  Full text retrieved: {fetched}")
    return records


# ── MAIN ─────────────────────────────────────────────────────────────────
def main():
    # Search queries
    pubmed_query = (
        '(microplastic*[TIAB] OR nanoplastic*[TIAB] OR "micro-plastic"[TIAB] OR "nano-plastic"[TIAB] '
        'OR "plastic particle"[TIAB] OR "plastic debris"[TIAB]) '
        'AND (human[TIAB] OR patient*[TIAB] OR "biological sample"[TIAB] OR biomonitoring[TIAB]) '
        'AND (detect*[TIAB] OR identif*[TIAB] OR quantif*[TIAB] OR measur*[TIAB] OR analyz*[TIAB] '
        'OR found[TIAB] OR present[TIAB] OR concentrat*[TIAB])'
    )

    epmc_query = (
        '(microplastic OR nanoplastic OR "micro-plastic" OR "nano-plastic" OR "plastic particle") '
        'AND (human OR patient OR "biological sample" OR biomonitoring) '
        'AND (detect OR identif OR quantif OR measur OR found OR present OR concentrat)'
    )

    openalex_query = 'microplastics nanoplastics human tissue detection identification blood lung placenta'

    s2_query = 'microplastics nanoplastics human tissue detection blood lung placenta gut'

    # Search all databases
    all_records = []

    pm_records = search_pubmed(pubmed_query, max_results=3000)
    all_records.extend(pm_records)

    epmc_records = search_europe_pmc(epmc_query, max_results=2000)
    all_records.extend(epmc_records)

    oa_records = search_openalex(openalex_query, max_results=2000)
    all_records.extend(oa_records)

    s2_records = search_semantic_scholar(s2_query, max_results=1000)
    all_records.extend(s2_records)

    # Deduplicate
    unique_records = deduplicate(all_records)

    # Fetch full text
    unique_records = fetch_fulltext_pmc(unique_records, max_fetch=500)

    # Save full text separately (large)
    fulltext_data = {}
    for rec in unique_records:
        ft = rec.pop('fulltext_xml', '')
        if ft:
            key = rec.get('pmc_id') or rec.get('pmid') or rec.get('doi', '')
            fulltext_data[key] = ft

    # Save
    df = pd.DataFrame(unique_records)
    df.to_csv(OUT / 'all_papers_raw.csv', index=False, encoding='utf-8-sig')

    # Save full text
    with open(OUT / 'fulltext_cache.json', 'w', encoding='utf-8') as f:
        json.dump(fulltext_data, f, ensure_ascii=False)

    # Save PRISMA counts
    prisma = {
        'pubmed_hits': len(pm_records),
        'europe_pmc_hits': len(epmc_records),
        'openalex_hits': len(oa_records),
        'semantic_scholar_hits': len(s2_records),
        'total_retrieved': len(all_records),
        'after_dedup': len(unique_records),
        'with_fulltext': sum(1 for r in unique_records if r.get('has_fulltext')),
        'with_abstract': sum(1 for r in unique_records if len(r.get('abstract', '')) > 50),
    }
    with open(OUT / 'prisma_counts.json', 'w') as f:
        json.dump(prisma, f, indent=2)

    print("\n" + "=" * 60)
    print("SEARCH SUMMARY")
    print("=" * 60)
    for k, v in prisma.items():
        print(f"  {k:<30}: {v}")
    print(f"\nSaved to {OUT}")

    return df, fulltext_data


if __name__ == '__main__':
    main()
