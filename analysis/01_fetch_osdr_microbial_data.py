#!/usr/bin/env python3
"""
01_fetch_osdr_microbial_data.py

Phase 1 of the OSD-528 re-analysis: harvest NASA OSDR metadata.

1. Raw study metadata for OSD-528 (Mycobacterium marinum; 3D clinostat vs RPM 2.0
   vs 1 g), kept for provenance. Sample groups live in config/samples.tsv.
2. A pan-microbial OSDR compendium (Supplementary Table S16), rebuilt from the
   OSDR API rather than from a hand-written list.

Compendium method (every decision is rule-based and written to a log):
  a. Enumerate studies.  The OSDR search index (type=cgene) is read in full, and
     https://osdr.nasa.gov/osdr/data/osd/meta/<n> is also probed for every
     accession number from 1 to (largest indexed number + SCAN_MARGIN), because
     some studies are served by the metadata API but are absent from the search
     index.
  b. Organisms come from each study's ISA "Organism" characteristic (sources
     and samples), together with the NCBI Taxonomy or MeSH term that OSDR
     attaches; the search-index "organism" field is the fallback.  "Host
     Organism" values are kept separately.
  c. Each organism is classified from its NCBI Taxonomy lineage (E-utilities
     efetch): Bacteria, Archaea, Viruses, Fungi, or another eukaryote outside
     Metazoa and Viridiplantae ("protist").  The MeSH "Microbiota" term, NCBI
     metagenome taxa, and the root taxon "cellular organisms" when the assay is
     a 16S/ITS amplicon or shotgun metagenome, count as a microbial community.
     Metazoa and Viridiplantae (including unicellular green algae) are not
     counted as microbial.  Names without an NCBI taxon term are looked up in
     NCBI Taxonomy (exact name, then synonyms, then genus); anything still
     unresolved is logged, never guessed.
  d. A study with at least one microbial organism is kept only if it addresses
     spaceflight or simulated microgravity: its OSDR Project Type is
     "Spaceflight Study"/"Spaceflight Project", or one of its factors (name,
     type or level) names spaceflight, microgravity, weightlessness,
     clinorotation, a random positioning machine, a rotating-wall / HARV vessel
     or low-shear modelled microgravity (RELEVANCE_* below).  Aircraft,
     high-altitude and parabolic-flight studies are therefore excluded unless a
     factor qualifies.  The basis is recorded per study, and every microbial
     study screened (kept or not) is written to the screening log.
  e. Title, organisms, assays and factors are recorded as the API returns them.
     A DOI is recorded only if OSDR's own "DOI" comment provides it and DataCite
     resolves it to a landing URL for that same OSD (or listed GLDS) accession;
     otherwise the DOI cell is left empty.

Raw API responses are cached (gzip) under data/raw/microbial_osdr_api_cache/
for the microbial studies, together with a compact index of every accession
scanned (data/raw/microbial_osdr_all_studies_index.json).  Re-running uses the
cache; --refresh re-downloads everything.

Usage:
  python analysis/01_fetch_osdr_microbial_data.py                    # OSD-528 metadata + Table S16
  python analysis/01_fetch_osdr_microbial_data.py --compendium-only  # Table S16 only
  python analysis/01_fetch_osdr_microbial_data.py --compendium-only --refresh
"""

import argparse
import csv
import datetime
import glob
import gzip
import html
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_RAW = os.path.join(PROJECT_DIR, 'data', 'raw')
DATA_PROCESSED = os.path.join(PROJECT_DIR, 'data', 'processed')
RESULTS_DIR = os.path.join(PROJECT_DIR, 'results_tables')
CACHE_DIR = os.path.join(DATA_RAW, 'microbial_osdr_api_cache')
INDEX_FILE = os.path.join(DATA_RAW, 'microbial_osdr_all_studies_index.json')
CATALOG_FILE = os.path.join(DATA_RAW, 'microbial_osdr_catalog.json')
SCREEN_LOG = os.path.join(DATA_PROCESSED, 'microbial_osdr_compendium_screening_log.tsv')
TABLE_PATTERN = 'Table_S16_pan_microbial_osdr_compendium_{n}studies.tsv'

# Plain User-Agent; no e-mail address is ever sent to any API.
HEADERS = {'User-Agent': 'osd528-reanalysis/1.0'}

OSDR_SEARCH = 'https://osdr.nasa.gov/osdr/data/search?type=cgene&size={size}&from={start}'
OSDR_META = 'https://osdr.nasa.gov/osdr/data/osd/meta/{n}'
NCBI_EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={ids}'
NCBI_ESEARCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&retmode=json&term={term}'
DATACITE = 'https://api.datacite.org/dois/{doi}'

SCAN_MARGIN = 150          # probe this many accession numbers past the largest indexed one
META_WORKERS = 3           # light parallelism for the OSDR metadata API

# NCBI Taxonomy lineage anchors
TAX_BACTERIA, TAX_ARCHAEA, TAX_VIRUSES, TAX_EUKARYOTA = 2, 2157, 10239, 2759
TAX_FUNGI, TAX_METAZOA, TAX_VIRIDIPLANTAE = 4751, 33208, 33090
TAX_METAGENOMES = 408169
TAX_CELLULAR_ORGANISMS = 131567
COMMUNITY_ASSAY = re.compile(r'amplicon|16S|\bITS\b|metagenom', re.I)
MESH_MICROBIOTA = 'D064307'
NOT_AN_ORGANISM = {'', 'not applicable', 'n/a', 'na', 'none'}
MICROBIAL_CATEGORIES = ('bacteria', 'archaea', 'fungus', 'protist', 'virus', 'microbial community')

RELEVANCE_PROJECT_TYPE = re.compile(r'^space ?flight( study| project)?$', re.I)
RELEVANCE_FACTOR = re.compile(
    r'space ?flight|micro-?gravity|weightless|clino|random positioning|\bRPM\b|'
    r'rotating[- ]wall|\bRWV\b|\bHARV\b|high aspect ratio vessel|low[- ]shear|\bLSMMG\b|'
    r'rotary cell culture', re.I)


def log(msg):
    print(msg, flush=True)


# --------------------------------------------------------------------------- HTTP

def http_get(url, retries=5, timeout=60):
    """GET with retries and exponential backoff.  Returns (status, bytes)."""
    delay = 2.0
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.status, resp.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return 404, b''
            last = e
        except Exception as e:  # network errors, timeouts
            last = e
        time.sleep(delay)
        delay *= 2
    raise RuntimeError(f'GET failed after {retries} attempts: {url}: {last}')


def fetch_json(url):
    """Legacy helper used by fetch_osd_metadata()."""
    try:
        status, body = http_get(url)
        return json.loads(body.decode('utf-8')) if status == 200 else None
    except Exception as e:
        print(f"Warning: Failed to fetch {url}: {e}")
        return None


def write_gz(path, data_bytes):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, 'wb') as f:
        f.write(data_bytes)


def read_gz(path):
    with gzip.open(path, 'rb') as f:
        return f.read()


def cached_get(url, cache_path, refresh=False):
    """Fetch url, caching the raw body (gzip) at cache_path."""
    if not refresh and os.path.exists(cache_path):
        return 200, read_gz(cache_path)
    status, body = http_get(url)
    if status == 200:
        write_gz(cache_path, body)
    return status, body


def fetch_osd_metadata(study_id):
    print(f"Fetching metadata for OSD-{study_id} from OSDR REST API...")
    url = f"https://osdr.nasa.gov/osdr/data/osd/meta/{study_id}"
    data = fetch_json(url)
    if data:
        out_file = os.path.join(DATA_RAW, f"OSD-{study_id}_metadata.json")
        with open(out_file, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)
        print(f"Saved: {out_file} ({os.path.getsize(out_file):,} bytes)")
        return data
    return None


# ------------------------------------------------------------ OSDR study parsing

def _annot(v):
    if isinstance(v, dict):
        return (v.get('annotationValue') or '').strip(), (v.get('termAccession') or '').strip()
    return (str(v).strip() if v is not None else ''), ''


def _uniq(seq):
    out = []
    for x in seq:
        if x not in out:
            out.append(x)
    return out


def extract_study(n, meta):
    """Compact, API-faithful record for one OSD study from its metadata JSON."""
    acc = f'OSD-{n}'
    st = meta['study'][acc]
    s = st['studies'][0]
    comments = {c.get('name', ''): (c.get('value') or '').strip() for c in s.get('comments', [])}

    cat_names = {c['@id']: (c.get('characteristicType', {}).get('annotationValue') or '').strip()
                 for c in s.get('characteristicCategories', [])}
    for a in s.get('assays', []):
        for c in a.get('characteristicCategories', []):
            cat_names.setdefault(c['@id'], (c.get('characteristicType', {}).get('annotationValue') or '').strip())
    organisms, hosts = [], []
    materials = s.get('materials', {})
    for group in ('sources', 'samples'):
        for m in materials.get(group, []):
            for ch in m.get('characteristics', []):
                cname = cat_names.get(ch.get('category', {}).get('@id', ''), '').lower()
                name, term = _annot(ch.get('value'))
                if name.lower() in NOT_AN_ORGANISM:
                    continue
                if cname == 'organism':
                    organisms.append((name, term))
                elif cname == 'host organism':
                    hosts.append((name, term))

    factor_names = {}
    factors = []
    for f in s.get('factors', []):
        fname = (f.get('factorName') or '').strip()
        ftype = (f.get('factorType', {}).get('annotationValue') or '').strip()
        factor_names[f.get('@id', '')] = fname
        factors.append({'name': fname, 'type': ftype, 'levels': []})
    by_name = {f['name']: f for f in factors}
    for m in materials.get('samples', []):
        for fv in m.get('factorValues', []):
            fname = factor_names.get(fv.get('category', {}).get('@id', ''))
            val, _ = _annot(fv.get('value'))
            if fname in by_name and val and val not in by_name[fname]['levels']:
                by_name[fname]['levels'].append(val)

    # Some responses leave the ISA factor values empty but carry the sample
    # table in additionalInformation; read the "Factor Value: <name>" columns.
    ai = st.get('additionalInformation') or {}
    for tbl in (ai.get('samples') or {}).values():
        for h in tbl.get('header', []):
            title = (h.get('title') or '') if isinstance(h, dict) else ''
            if not title.lower().startswith('factor value:'):
                continue
            fname = title.split(':', 1)[1].strip().lower()
            target = next((f for f in factors if f['name'].lower() == fname and not f['levels']), None)
            if target is None:
                continue
            for row in tbl.get('table', []):
                val = html.unescape(re.sub(r'<[^>]+>', '', str(row.get(h.get('field'), '')))).strip()
                if val and val not in target['levels']:
                    target['levels'].append(val)

    assays = []
    for a in s.get('assays', []):
        mt, _ = _annot(a.get('measurementType'))
        tt, _ = _annot(a.get('technologyType'))
        tp = a.get('technologyPlatform')
        tp = (tp.get('annotationValue') if isinstance(tp, dict) else tp) or ''
        label = ' / '.join(x for x in (mt, tt) if x)
        if tp:
            label += f' ({tp.strip()})'
        assays.append(label)

    return {
        'accession': acc,
        'title': (s.get('title') or '').strip(),
        'organisms': [{'name': nme, 'term': trm} for nme, trm in _uniq(organisms)],
        'host_organisms': [{'name': nme, 'term': trm} for nme, trm in _uniq(hosts)],
        'assays': _uniq([a for a in assays if a]),
        'factors': factors,
        'project_type': comments.get('Project Type', ''),
        'flight_program': comments.get('Flight Program', ''),
        'mission_name': comments.get('Mission Name', ''),
        'experiment_platform': comments.get('Experiment Platform', ''),
        'identifiers': comments.get('Identifiers', ''),
        'osdr_doi': comments.get('DOI', ''),
        'osdr_version': st.get('version'),
    }


def clean(value):
    """Collapse whitespace (incl. tabs/newlines) so every TSV cell is one line."""
    return re.sub(r'\s+', ' ', str(value if value is not None else '')).strip()


def split_search_organisms(value):
    parts = re.split(r'\s{2,}|;', str(value or ''))
    return [p.strip() for p in parts if p.strip() and p.strip().lower() not in NOT_AN_ORGANISM]


# ---------------------------------------------------------------- enumeration

def load_search_index(refresh):
    path = os.path.join(CACHE_DIR, 'search_cgene_all.json.gz')
    size, start, hits, total = 1000, 0, [], None
    if not refresh and os.path.exists(path):
        data = json.loads(read_gz(path))
        return data['hits']['hits'], data['hits']['total']
    while True:
        status, body = http_get(OSDR_SEARCH.format(size=size, start=start))
        page = json.loads(body)
        total = page['hits']['total']
        total = total['value'] if isinstance(total, dict) else total
        hits.extend(page['hits']['hits'])
        start += size
        if start >= total or not page['hits']['hits']:
            break
    write_gz(path, json.dumps({'hits': {'total': total, 'hits': hits}}).encode())
    return hits, total


def scan_metadata(numbers):
    """Fetch metadata for each accession number. Returns {n: raw_bytes or None}."""
    results, errors = {}, {}

    def one(n):
        try:
            status, body = http_get(OSDR_META.format(n=n))
            if status != 200:
                return n, None, None
            d = json.loads(body)
            return n, (body if d.get('hits', 0) and d.get('study') else None), None
        except Exception as e:
            return n, None, str(e)

    with ThreadPoolExecutor(max_workers=META_WORKERS) as pool:
        for i, (n, body, err) in enumerate(pool.map(one, numbers), 1):
            if err:
                errors[n] = err
            else:
                results[n] = body
            if i % 100 == 0:
                log(f'  scanned {i}/{len(numbers)} accession numbers')
    return results, errors


# ------------------------------------------------------------- NCBI taxonomy

def _taxid_from_term(term):
    m = re.search(r'NCBITAXON/(\d+)$', term or '')
    return int(m.group(1)) if m else None


def load_json(path, default):
    if os.path.exists(path):
        with open(path, encoding='utf-8') as f:
            return json.load(f)
    return default


def save_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(obj, f, indent=1, ensure_ascii=False, sort_keys=True)
        f.write('\n')


def _esearch_ids(term):
    _, body = http_get(NCBI_ESEARCH.format(term=urllib.parse.quote(term)))
    time.sleep(0.4)   # NCBI: <= 3 requests/s without an API key
    return json.loads(body).get('esearchresult', {}).get('idlist', [])


def resolve_names_to_taxids(names, refresh):
    """NCBI Taxonomy lookup for organism names that OSDR gave no NCBI taxon term.

    Tried in order, accepting only a single unambiguous hit:
      1. exact scientific name;
      2. any taxonomy name field (catches synonyms, e.g. renamed species);
      3. the genus (first word) as an exact scientific name of rank genus.
    The method used is cached with the result; failures stay unresolved."""
    path = os.path.join(CACHE_DIR, 'ncbi_taxonomy_name_lookup.json')
    cache = {} if refresh else load_json(path, {})
    for name in sorted(set(names)):
        if name in cache:
            continue
        result = {'taxid': None, 'method': 'unresolved'}
        if re.search(r'[A-Za-z]{3,}', name):
            tries = [('exact scientific name', f'"{name}"[Scientific Name]'),
                     ('taxonomy name search (synonyms)', f'"{name}"')]
            genus = name.split()[0]
            if len(name.split()) > 1 and re.fullmatch(r'[A-Z][a-z]+', genus):
                tries.append(('genus of listed name', f'"{genus}"[Scientific Name] AND genus[Rank]'))
            for method, term in tries:
                ids = _esearch_ids(term)
                if len(ids) == 1:
                    result = {'taxid': int(ids[0]), 'method': method}
                    break
        cache[name] = result
    save_json(path, cache)
    return {k: v['taxid'] for k, v in cache.items()}, cache


def fetch_lineages(taxids, refresh):
    """taxid -> {name, rank, lineage_ids, lineage}.  Raw efetch XML is cached."""
    path = os.path.join(CACHE_DIR, 'ncbi_taxonomy_lineages.json')
    cache = {} if refresh else {int(k): v for k, v in load_json(path, {}).items()}
    todo = sorted({t for t in taxids if t is not None and t not in cache})
    for i in range(0, len(todo), 100):
        batch = todo[i:i + 100]
        _, body = http_get(NCBI_EFETCH.format(ids=','.join(map(str, batch))))
        write_gz(os.path.join(CACHE_DIR, 'ncbi_taxonomy_efetch', f'efetch_{batch[0]}_{batch[-1]}.xml.gz'), body)
        root = ET.fromstring(body)
        for tx in root.findall('Taxon'):
            tid = int(tx.findtext('TaxId'))
            rec = {
                'name': tx.findtext('ScientificName'),
                'rank': tx.findtext('Rank'),
                'lineage_ids': [int(t.findtext('TaxId')) for t in tx.findall('LineageEx/Taxon')] + [tid],
                'lineage': tx.findtext('Lineage'),
            }
            cache[tid] = rec
            for aka in tx.findall('AkaTaxIds/TaxId'):
                cache[int(aka.text)] = dict(rec, merged_into=tid)
        time.sleep(0.4)
    save_json(path, {str(k): v for k, v in cache.items()})
    return cache


def classify_taxon(rec):
    if rec is None:
        return 'unresolved'
    ids = set(rec['lineage_ids'])
    if TAX_BACTERIA in ids:
        return 'bacteria'
    if TAX_ARCHAEA in ids:
        return 'archaea'
    if TAX_VIRUSES in ids:
        return 'virus'
    if TAX_METAGENOMES in ids:
        return 'microbial community'
    if TAX_EUKARYOTA in ids:
        if TAX_FUNGI in ids:
            return 'fungus'
        if TAX_METAZOA in ids:
            return 'animal'
        if TAX_VIRIDIPLANTAE in ids:
            return 'plant'
        return 'protist'
    return 'unresolved'


def classify_organism(org, name_lookup, lineages, assays=()):
    term = org.get('term', '')
    if term.endswith(MESH_MICROBIOTA) or org['name'].strip().lower() == 'microbiota':
        return 'microbial community', None
    tid = _taxid_from_term(term)
    if tid is None:
        tid = name_lookup.get(org['name'])
    if tid == TAX_CELLULAR_ORGANISMS:
        # NCBI root-level "cellular organisms": a community sample, which only
        # counts as microbial when the assay profiles microbes (16S/ITS
        # amplicons or shotgun metagenomes).
        if any(COMMUNITY_ASSAY.search(a) for a in assays):
            return 'microbial community', tid
        return 'unresolved', tid
    return (classify_taxon(lineages.get(tid)) if tid is not None else 'unresolved'), tid


# ----------------------------------------------------------------- DataCite

def check_doi(doi, study, refresh):
    """Return (ok, datacite_url, datacite_title, reason)."""
    if not doi:
        return False, '', '', 'no DOI in OSDR metadata'
    safe = re.sub(r'[^A-Za-z0-9._-]', '_', doi)
    path = os.path.join(CACHE_DIR, 'datacite', f'{safe}.json.gz')
    status, body = cached_get(DATACITE.format(doi=urllib.parse.quote(doi, safe='/')), path, refresh)
    if status != 200:
        return False, '', '', f'DataCite HTTP {status}'
    attr = json.loads(body)['data']['attributes']
    url = attr.get('url') or ''
    title = '; '.join(t.get('title', '') for t in attr.get('titles', []))
    n = study['accession'].split('-')[1]
    accepted = {f'OSD-{n}'} | set(re.findall(r'GLDS-\d+', study.get('identifiers', '')))
    m = re.search(r'((?:OSD|GLDS)-\d+)/?$', url)
    if m and m.group(1) in accepted:
        return True, url, title, 'resolves to this study'
    return False, url, title, f'DataCite URL does not point to {study["accession"]}'


# ----------------------------------------------------------------- compendium

def relevance(study):
    reasons = []
    pt = study.get('project_type', '').strip()
    if RELEVANCE_PROJECT_TYPE.match(pt):
        reasons.append(f'Project Type = {pt}')
    for f in study['factors']:
        hits = [x for x in [f['name'], f['type']] + f['levels'] if x and RELEVANCE_FACTOR.search(x)]
        if hits:
            reasons.append(f"factor {f['name']}: " + ', '.join(_uniq(hits)))
    return reasons


def fmt_factors(study):
    out = []
    for f in study['factors']:
        label = f['name']
        if f['type'] and f['type'].lower() != f['name'].lower():
            label += f" [{f['type']}]"
        if f['levels']:
            label += ' (' + ', '.join(f['levels']) + ')'
        out.append(label)
    return '; '.join(out)


def build_compendium(refresh=False):
    log('=== Pan-microbial OSDR compendium (Table S16) ===')
    retrieved = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
    index = None if refresh else load_json(INDEX_FILE, None)
    raw_meta = {}

    if index is None:
        hits, total = load_search_index(refresh=True)
        search = {}
        for h in hits:
            src = h.get('_source', {})
            acc = src.get('Accession') or ''
            if acc.startswith('OSD-') and acc[4:].isdigit():
                search[acc] = src
        log(f'OSDR search index: {total} hits, {len(search)} unique OSD accessions')
        max_n = max(int(a[4:]) for a in search)
        numbers = list(range(1, max_n + SCAN_MARGIN + 1))
        log(f'Probing metadata API for OSD-1 .. OSD-{numbers[-1]} ...')
        raw, errors = scan_metadata(numbers)
        if errors:
            for n, e in sorted(errors.items()):
                log(f'  ERROR OSD-{n}: {e}')
            sys.exit(f'{len(errors)} metadata requests failed; not writing an incomplete table. Re-run.')
        studies, empty = {}, []
        for n in numbers:
            body = raw.get(n)
            acc = f'OSD-{n}'
            if body is None:
                empty.append(n)
                if acc in search:
                    log(f'  note: {acc} is in the search index but the metadata API returned no study')
                continue
            rec = extract_study(n, json.loads(body))
            rec['in_search_index'] = acc in search
            if acc in search and not rec['identifiers']:
                rec['identifiers'] = str(search[acc].get('Identifiers') or '').strip()
            if not rec['organisms'] and acc in search:
                rec['organisms'] = [{'name': x, 'term': ''} for x in split_search_organisms(search[acc].get('organism'))]
                rec['organism_source'] = 'search index'
            else:
                rec['organism_source'] = 'ISA Organism characteristic'
            studies[acc] = rec
            raw_meta[acc] = body
        index = {
            'retrieved_utc': retrieved,
            'search_index_hits': total,
            'search_index_unique_osd': len(search),
            'scanned_numbers': [numbers[0], numbers[-1]],
            'empty_numbers': empty,
            'studies': studies,
        }
    else:
        log(f"Using cached study index from {index['retrieved_utc']} ({len(index['studies'])} studies); --refresh to re-download")

    studies = index['studies']
    # taxonomy
    need_names = [o['name'] for s in studies.values() for o in s['organisms'] + s['host_organisms']
                  if _taxid_from_term(o.get('term')) is None and not o.get('term', '').endswith(MESH_MICROBIOTA)
                  and o['name'].strip().lower() != 'microbiota']
    name_lookup, name_methods = resolve_names_to_taxids(need_names, refresh)
    taxids = [_taxid_from_term(o.get('term')) for s in studies.values() for o in s['organisms'] + s['host_organisms']]
    taxids += [t for t in name_lookup.values() if t]
    lineages = fetch_lineages(taxids, refresh)

    for s in studies.values():
        for o in s['organisms'] + s['host_organisms']:
            o['category'], o['taxid'] = classify_organism(o, name_lookup, lineages, s['assays'])
            if o['taxid'] is not None and _taxid_from_term(o.get('term')) is None:
                o['taxid_lookup'] = name_methods.get(o['name'], {}).get('method')
        s['microbial_categories'] = _uniq([o['category'] for o in s['organisms'] if o['category'] in MICROBIAL_CATEGORIES])
    save_json(INDEX_FILE, index)

    microbial = {a: s for a, s in studies.items() if s['microbial_categories']}
    unresolved = sorted({(a, o['name']) for a, s in studies.items() for o in s['organisms'] if o['category'] == 'unresolved'})
    log(f"{len(studies)} studies scanned; {len(microbial)} list at least one microbial organism; "
        f"{len(unresolved)} organism entries unresolved")

    # cache raw metadata for every microbial study (provenance for each table cell)
    meta_dir = os.path.join(CACHE_DIR, 'meta')
    if raw_meta:   # fresh scan: drop files left over from earlier rule sets
        for stale in glob.glob(os.path.join(meta_dir, 'OSD-*.json.gz')):
            if os.path.basename(stale)[:-8] not in microbial:
                os.remove(stale)
    for acc in microbial:
        path = os.path.join(meta_dir, f'{acc}.json.gz')
        if acc in raw_meta:
            write_gz(path, raw_meta[acc])
        elif not os.path.exists(path):
            cached_get(OSDR_META.format(n=acc[4:]), path, refresh=True)

    rows, log_rows = [], []
    for acc in sorted(microbial, key=lambda a: int(a[4:])):
        s = microbial[acc]
        basis = relevance(s)
        ok, dc_url, dc_title, doi_reason = check_doi(s['osdr_doi'], s, refresh)
        decision = 'included' if basis else 'excluded: no spaceflight / simulated-microgravity project type or factor'
        log_rows.append({
            'accession': acc, 'decision': decision, 'inclusion_basis': ' | '.join(basis),
            'title': s['title'], 'organisms': '; '.join(f"{o['name']} [{o['category']}]" for o in s['organisms']),
            'project_type': s['project_type'], 'factors': fmt_factors(s),
            'osdr_doi': s['osdr_doi'], 'doi_check': doi_reason, 'datacite_url': dc_url, 'datacite_title': dc_title,
            'in_search_index': s.get('in_search_index'),
        })
        if not basis:
            continue
        rows.append({
            'accession': acc,
            'title': s['title'],
            'organisms': '; '.join(o['name'] for o in s['organisms']),
            'microbial_categories': '; '.join(s['microbial_categories']),
            'other_listed_organisms': '; '.join(f"{o['name']} [{o['category']}]" for o in s['organisms']
                                                if o['category'] not in MICROBIAL_CATEGORIES),
            'host_organisms': '; '.join(o['name'] for o in s['host_organisms']),
            'assay_types': '; '.join(s['assays']),
            'factors': fmt_factors(s),
            'project_type': s['project_type'],
            'flight_program': s['flight_program'],
            'mission': s['mission_name'],
            'inclusion_basis': ' | '.join(basis),
            'doi': s['osdr_doi'] if ok else '',
        })

    rows = [{k: clean(v) for k, v in r.items()} for r in rows]
    log_rows = [{k: clean(v) for k, v in r.items()} for r in log_rows]
    with open(SCREEN_LOG, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(log_rows[0].keys()), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(log_rows)
        for acc, name in unresolved:
            f.write(f'# organism entry not resolved in NCBI Taxonomy (any study): {acc}\t{name}\n')

    catalog = {
        'source': 'NASA OSDR API (osdr.nasa.gov/osdr/data/search, /osdr/data/osd/meta/<n>); '
                  'built by analysis/01_fetch_osdr_microbial_data.py',
        'retrieved_utc': index['retrieved_utc'],
        'n_studies': len(rows),
        'studies': {r['accession']: dict(microbial[r['accession']], inclusion_basis=r['inclusion_basis'],
                                          doi_verified=r['doi']) for r in rows},
    }
    save_json(CATALOG_FILE, catalog)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    table = os.path.join(RESULTS_DIR, TABLE_PATTERN.format(n=len(rows)))
    for old in glob.glob(os.path.join(RESULTS_DIR, TABLE_PATTERN.format(n='*'))):
        if os.path.abspath(old) != os.path.abspath(table):
            os.remove(old)
            log(f'Removed stale {os.path.basename(old)}')
    with open(table, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    n_doi = sum(1 for r in rows if r['doi'])
    log(f'Wrote {table}: {len(rows)} studies, {n_doi} DataCite-verified DOIs')
    log(f'Wrote {SCREEN_LOG} ({len(log_rows)} microbial studies screened)')
    log(f'Wrote {CATALOG_FILE}')
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--compendium-only', action='store_true',
                    help='only rebuild the pan-microbial compendium (Table S16)')
    ap.add_argument('--refresh', action='store_true', help='ignore the API cache and re-download')
    args = ap.parse_args()
    os.makedirs(DATA_RAW, exist_ok=True)
    os.makedirs(DATA_PROCESSED, exist_ok=True)
    print("=== Phase 1: Fetching NASA OSDR Microbial Microgravity Data ===")
    if not args.compendium_only:
        fetch_osd_metadata(528)  # raw OSD-528 study metadata, kept for provenance
    build_compendium(refresh=args.refresh)
    print("Phase 1 completed successfully.")


if __name__ == '__main__':
    main()
