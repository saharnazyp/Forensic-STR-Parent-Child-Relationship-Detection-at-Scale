# src/codechallenge2025/participant_solution.py
"""
Easy Participant Template for #codechallenge2025

You ONLY need to implement the function: match_single

The find_matches function is provided for you — no need to change it!
"""

import pandas as pd
from typing import List, Dict, Any


def match_single(
    query_profile: Dict[str, Any], database_df: pd.DataFrame
) -> List[Dict]:
    """
    Find the top 10 candidate matches for a SINGLE query profile.

    Args:
        query_profile: dict with 'PersonID' and locus columns (e.g. {'PersonID': 'Q001', 'TH01': '9,9.3', ...})
        database_df: Full database as pandas DataFrame (500k rows)

    Returns:
        List of up to 10 candidate dicts, sorted by strength (best first):
        [
            {
                "person_id": "P000123",
                "clr": 1e15,                    # Combined Likelihood Ratio
                "posterior": 0.99999,           # Optional: posterior probability
                "consistent_loci": 20,
                "mutated_loci": 1,
                "inconclusive_loci": 0
            },
            ...
        ]
    """
    # === Helper: parse alleles safely ===
    def parse_alleles(val):
        if pd.isna(val):
            return None
        s = str(val).strip()
        if s in ('-', '', 'nan', 'None'):
            return None
        try:
            if ',' in s:
                return frozenset(float(x.strip()) for x in s.split(','))
            return frozenset([float(s)])
        except (ValueError, TypeError):
            return None

    # === Build/retrieve cached index (using function attribute) ===
    db_id = id(database_df)
    if not hasattr(match_single, '_cache') or match_single._cache.get('db_id') != db_id:
        # Build allele index and frequency table
        loci = [c for c in database_df.columns if c != 'PersonID']
        allele_index = {}  # (locus, allele) -> set of person_ids
        allele_counts = {}  # (locus, allele) -> count
        profiles = {}  # person_id -> {locus: frozenset of alleles}

        for _, row in database_df.iterrows():
            pid = row['PersonID']
            profile = {}
            for locus in loci:
                alleles = parse_alleles(row[locus])
                profile[locus] = alleles
                if alleles:
                    for a in alleles:
                        key = (locus, a)
                        if key not in allele_index:
                            allele_index[key] = set()
                            allele_counts[key] = 0
                        allele_index[key].add(pid)
                        allele_counts[key] += 1
            profiles[pid] = profile

        # Compute frequencies
        total_profiles = len(database_df)
        allele_freqs = {k: v / total_profiles for k, v in allele_counts.items()}

        match_single._cache = {
            'db_id': db_id,
            'loci': loci,
            'allele_index': allele_index,
            'allele_freqs': allele_freqs,
            'profiles': profiles
        }

    cache = match_single._cache
    loci = cache['loci']
    allele_index = cache['allele_index']
    allele_freqs = cache['allele_freqs']
    profiles = cache['profiles']

    query_id = query_profile['PersonID']

    # Parse query profile
    query_parsed = {}
    for locus in loci:
        query_parsed[locus] = parse_alleles(query_profile.get(locus))

    # === Pre-filter candidates using inverted index ===
    candidate_match_count = {}
    for locus in loci:
        q_alleles = query_parsed[locus]
        if not q_alleles:
            continue
        for allele in q_alleles:
            key = (locus, allele)
            if key in allele_index:
                for pid in allele_index[key]:
                    if pid != query_id:
                        candidate_match_count[pid] = candidate_match_count.get(pid, 0) + 1

    # Only consider candidates that share alleles at multiple loci
    promising = [pid for pid, cnt in candidate_match_count.items() if cnt >= 8]

    # === Score promising candidates ===
    results = []

    for pid in promising:
        cand_profile = profiles[pid]

        clr = 1.0
        consistent = 0
        mutated = 0
        inconclusive = 0
        exclusions = 0
        identity_count = 0
        compared = 0

        for locus in loci:
            q_alleles = query_parsed[locus]
            c_alleles = cand_profile.get(locus)

            if q_alleles is None or c_alleles is None:
                inconclusive += 1
                continue

            compared += 1
            if q_alleles == c_alleles:
                identity_count += 1

            shared = q_alleles & c_alleles

            if shared:
                consistent += 1
                # Simple LR: transmission_prob / frequency
                trans = 1.0 if len(c_alleles) == 1 else 0.5
                clr *= trans / 0.15  # Average frequency
            else:
                # Check mutation
                is_mutation = any(0 < abs(a - b) <= 1.0 for a in q_alleles for b in c_alleles)
                if is_mutation:
                    mutated += 1
                    clr *= 0.01  # Mutation penalty
                else:
                    exclusions += 1
                    clr *= 0.001  # Exclusion penalty

        # Filter criteria
        # True parent-child should have 0 exclusions (rarely 1 due to mutation)
        if exclusions > 1:
            continue
        if consistent < 8:
            continue
        if compared > 0 and identity_count / compared > 0.85:
            continue

        posterior = clr / (clr + 1.0) if clr > 0 else 0.0

        results.append({
            "person_id": pid,
            "clr": clr,
            "posterior": posterior,
            "consistent_loci": consistent,
            "mutated_loci": mutated,
            "inconclusive_loci": inconclusive
        })

    results.sort(key=lambda x: -x['clr'])
    return results[:10]


# ============================================================
# DO NOT MODIFY BELOW THIS LINE — This runs your function!
# ============================================================


def find_matches(database_path: str, queries_path: str) -> List[Dict]:
    """
    Main entry point — automatically tested by CI.
    Loads data and calls your match_single for each query.
    """
    print("Loading database and queries...")
    database_df = pd.read_csv(database_path)
    queries_df = pd.read_csv(queries_path)

    results = []

    print(f"Processing {len(queries_df)} queries...")
    for _, query_row in queries_df.iterrows():
        query_id = query_row["PersonID"]
        query_profile = query_row.to_dict()

        print(f"  Matching query {query_id}...")
        top_candidates = match_single(query_profile, database_df)

        results.append(
            {
                "query_id": query_id,
                "top_candidates": top_candidates[:10],  # Ensure max 10
            }
        )

    print("All queries processed.")
    return results
