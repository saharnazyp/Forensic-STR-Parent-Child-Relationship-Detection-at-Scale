# src/codechallenge2025/participant_solution.py
"""
Easy Participant Template for #codechallenge2025

You ONLY need to implement the function: match_single

The find_matches function is provided for you — no need to change it!
"""

import pandas as pd
from typing import List, Dict, Any


# Allele frequencies (from dataset generator)
ALLELE_FREQS = {
    "D3S1358": {14: 0.15, 15: 0.25, 16: 0.22, 17: 0.20, 18: 0.13, 19: 0.05},
    "vWA": {14: 0.10, 15: 0.12, 16: 0.20, 17: 0.25, 18: 0.20, 19: 0.10, 20: 0.03},
    "FGA": {19: 0.05, 20: 0.10, 21: 0.15, 22: 0.20, 23: 0.18, 24: 0.15, 25: 0.10, 26: 0.07},
    "D8S1179": {10: 0.05, 11: 0.08, 12: 0.10, 13: 0.30, 14: 0.25, 15: 0.15, 16: 0.07},
    "D21S11": {27: 0.05, 28: 0.15, 29: 0.20, 30: 0.25, 31: 0.15, 32: 0.10, 30.2: 0.08, 31.2: 0.02},
    "D18S51": {12: 0.08, 13: 0.15, 14: 0.20, 15: 0.18, 16: 0.12, 17: 0.10, 18: 0.08, 19: 0.06, 20: 0.03},
    "D5S818": {9: 0.05, 10: 0.08, 11: 0.25, 12: 0.30, 13: 0.20, 14: 0.10, 15: 0.02},
    "D13S317": {8: 0.05, 9: 0.08, 10: 0.10, 11: 0.25, 12: 0.20, 13: 0.18, 14: 0.12, 15: 0.02},
    "D7S820": {8: 0.10, 9: 0.12, 10: 0.25, 11: 0.28, 12: 0.15, 13: 0.08, 14: 0.02},
    "D16S539": {8: 0.05, 9: 0.20, 10: 0.15, 11: 0.25, 12: 0.20, 13: 0.10, 14: 0.05},
    "TH01": {6: 0.20, 7: 0.15, 8: 0.18, 9: 0.22, 9.3: 0.15, 10: 0.08, 11: 0.02},
    "TPOX": {8: 0.40, 9: 0.10, 10: 0.12, 11: 0.25, 12: 0.10, 13: 0.03},
    "CSF1PO": {9: 0.05, 10: 0.20, 11: 0.25, 12: 0.30, 13: 0.12, 14: 0.08},
    "D2S1338": {17: 0.08, 18: 0.05, 19: 0.10, 20: 0.15, 21: 0.08, 22: 0.07, 23: 0.12, 24: 0.15, 25: 0.15},
    "D19S433": {13: 0.15, 14: 0.30, 14.2: 0.05, 15: 0.20, 15.2: 0.05, 16: 0.15, 17: 0.10},
    "D22S1045": {11: 0.10, 14: 0.08, 15: 0.30, 16: 0.35, 17: 0.12, 18: 0.05},
    "D10S1248": {11: 0.05, 12: 0.08, 13: 0.25, 14: 0.30, 15: 0.20, 16: 0.10, 17: 0.02},
    "D1S1656": {12: 0.10, 13: 0.08, 14: 0.05, 15: 0.12, 16: 0.15, 17: 0.20, 17.3: 0.10, 18: 0.10, 18.3: 0.05},
    "D12S391": {17: 0.05, 18: 0.15, 19: 0.12, 20: 0.20, 21: 0.18, 22: 0.15, 23: 0.10, 24: 0.05},
    "D2S441": {10: 0.10, 11: 0.20, 11.3: 0.05, 12: 0.08, 13: 0.10, 14: 0.25, 15: 0.15, 16: 0.07},
    "SE33": {19: 0.05, 20: 0.08, 21: 0.10, 22: 0.12, 23: 0.10, 24: 0.08, 25: 0.12, 26: 0.10, 27: 0.10, 28: 0.08, 29: 0.07},
}

DEFAULT_FREQ = 0.01


def get_freq(locus: str, allele: float) -> float:
    """Get allele frequency."""
    if locus in ALLELE_FREQS:
        return ALLELE_FREQS[locus].get(allele, DEFAULT_FREQ)
    return DEFAULT_FREQ


def normalize_str(value) -> str:
    """Normalize allele string for comparison."""
    if pd.isna(value) or str(value).strip() in ("-", ""):
        return "-"
    return str(value).strip()


def parse_alleles(value) -> set:
    """Parse allele string into a set of float values."""
    if pd.isna(value) or str(value).strip() in ("-", ""):
        return set()
    
    parts = str(value).strip().split(",")
    alleles = set()
    for p in parts:
        p = p.strip()
        if p and p != "-":
            try:
                alleles.add(float(p))
            except ValueError:
                pass
    return alleles


def analyze_locus(query_str, candidate_str, locus: str) -> dict:
    """
    Analyze one locus.
    
    Returns:
    - type: "identical", "partial", "mutation", "mismatch", "missing"
    - score: contribution to CLR
    """
    q_norm = normalize_str(query_str)
    c_norm = normalize_str(candidate_str)
    
    # Missing data
    if q_norm == "-" or c_norm == "-":
        return {"type": "missing", "score": 1.0}
    
    q_alleles = parse_alleles(query_str)
    c_alleles = parse_alleles(candidate_str)
    
    # Check if they share at least one allele
    shared = q_alleles & c_alleles
    
    if not shared:
        # Check for mutation (±1 step)
        for q in q_alleles:
            for c in c_alleles:
                if 0 < abs(q - c) <= 1.0:
                    return {"type": "mutation", "score": 0.1}
        # Complete mismatch
        return {"type": "mismatch", "score": 0.0001}
    
    # They share allele(s)
    # Calculate LR based on rarest shared allele
    min_freq = min(get_freq(locus, a) for a in shared)
    lr = 1.0 / min_freq
    
    # Check if strings are identical or different
    if q_norm == c_norm:
        return {"type": "identical", "score": lr}
    else:
        # Different strings but shared allele = strong parent-child signal
        return {"type": "partial", "score": lr}


def match_single(
    query_profile: Dict[str, Any], database_df: pd.DataFrame
) -> List[Dict]:
    """
    Find the top 10 candidate matches for a SINGLE query profile.
    
    Version 7: Better same-person vs parent-child distinction.
    """
    query_id = query_profile["PersonID"]
    loci = [col for col in database_df.columns if col != "PersonID"]
    
    candidates = []
    
    for idx, row in database_df.iterrows():
        candidate_id = row["PersonID"]
        
        if candidate_id == query_id:
            continue
        
        total_lr = 1.0
        identical_count = 0
        partial_count = 0
        mismatch_count = 0
        mutation_count = 0
        missing_count = 0
        
        for locus in loci:
            q_val = query_profile.get(locus, "-")
            c_val = row[locus]
            
            result = analyze_locus(q_val, c_val, locus)
            total_lr *= result["score"]
            
            if result["type"] == "identical":
                identical_count += 1
            elif result["type"] == "partial":
                partial_count += 1
            elif result["type"] == "mismatch":
                mismatch_count += 1
            elif result["type"] == "mutation":
                mutation_count += 1
            else:
                missing_count += 1
        
        # Skip if too many mismatches (not related)
        if mismatch_count > 2:
            continue
        
        # Key insight: 
        # - Same person (C000xxx): almost all loci are "identical", very few "partial"
        # - Parent (P000xxx): many loci are "partial" (different strings, shared allele)
        
        # Require minimum partial matches to be considered parent-child
        # With 21 loci, ~5% dropout, ~8% single allele, expect ~15-18 comparable loci
        # Parent-child should have many partial matches (child's other allele differs)
        
        if partial_count < 3:
            # Too few partial matches - likely same person, not parent-child
            continue
        
        # Score: base LR + bonus for partial matches
        # Partial matches are the KEY indicator of parent-child relationship
        score = total_lr * (1 + partial_count * 0.5)
        
        candidates.append({
            "person_id": candidate_id,
            "clr": score,
            "posterior": partial_count / 21.0,
            "consistent_loci": identical_count + partial_count,
            "mutated_loci": mutation_count,
            "inconclusive_loci": missing_count,
        })
    
    # Sort by score (highest first)
    candidates.sort(key=lambda x: x["clr"], reverse=True)
    return candidates[:10]


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
