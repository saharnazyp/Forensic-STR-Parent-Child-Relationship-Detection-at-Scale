# src/codechallenge2025/dataset_generator.py
"""
Synthetic STR Dataset Generator for #codechallenge2025
Generates realistic forensic DNA profiles with hidden parent-child relationships.
"""

import random
import pandas as pd
import os

# -------------------------------
# Configuration
# -------------------------------

NUM_DB_PROFILES = 500000  # Total mixed profiles in database
NUM_QUERIES = 40  # Number of query profiles
NUM_TRUE_PAIRS = 35  # Number of queries with a true match in DB

# 21 common forensic loci (CODIS + expanded)
LOCI = [
    "D3S1358",
    "vWA",
    "FGA",
    "D8S1179",
    "D21S11",
    "D18S51",
    "D5S818",
    "D13S317",
    "D7S820",
    "D16S539",
    "TH01",
    "TPOX",
    "CSF1PO",
    "D2S1338",
    "D19S433",
    "D22S1045",
    "D10S1248",
    "D1S1656",
    "D12S391",
    "D2S441",
    "SE33",
]

# Dropout & mutation rates
DROPOUT_RATE = 0.05  # Complete locus missing
SINGLE_ALLELE_RATE = 0.08  # Only one allele observed
MUTATION_RATE = 0.002  # Per locus per generation

# -------------------------------
# Allele frequencies (based on real population data, includes microvariants)
# -------------------------------

ALLELE_FREQS = {
    "D3S1358": {14: 0.15, 15: 0.25, 16: 0.22, 17: 0.20, 18: 0.13, 19: 0.05},
    "vWA": {14: 0.10, 15: 0.12, 16: 0.20, 17: 0.25, 18: 0.20, 19: 0.10, 20: 0.03},
    "FGA": {
        19: 0.05,
        20: 0.10,
        21: 0.15,
        22: 0.20,
        23: 0.18,
        24: 0.15,
        25: 0.10,
        26: 0.07,
    },
    "D8S1179": {10: 0.05, 11: 0.08, 12: 0.10, 13: 0.30, 14: 0.25, 15: 0.15, 16: 0.07},
    "D21S11": {
        27: 0.05,
        28: 0.15,
        29: 0.20,
        30: 0.25,
        31: 0.15,
        32: 0.10,
        30.2: 0.08,
        31.2: 0.02,
    },
    "D18S51": {
        12: 0.08,
        13: 0.15,
        14: 0.20,
        15: 0.18,
        16: 0.12,
        17: 0.10,
        18: 0.08,
        19: 0.06,
        20: 0.03,
    },
    "D5S818": {9: 0.05, 10: 0.08, 11: 0.25, 12: 0.30, 13: 0.20, 14: 0.10, 15: 0.02},
    "D13S317": {
        8: 0.05,
        9: 0.08,
        10: 0.10,
        11: 0.25,
        12: 0.20,
        13: 0.18,
        14: 0.12,
        15: 0.02,
    },
    "D7S820": {8: 0.10, 9: 0.12, 10: 0.25, 11: 0.28, 12: 0.15, 13: 0.08, 14: 0.02},
    "D16S539": {8: 0.05, 9: 0.20, 10: 0.15, 11: 0.25, 12: 0.20, 13: 0.10, 14: 0.05},
    "TH01": {6: 0.20, 7: 0.15, 8: 0.18, 9: 0.22, 9.3: 0.15, 10: 0.08, 11: 0.02},
    "TPOX": {8: 0.40, 9: 0.10, 10: 0.12, 11: 0.25, 12: 0.10, 13: 0.03},
    "CSF1PO": {9: 0.05, 10: 0.20, 11: 0.25, 12: 0.30, 13: 0.12, 14: 0.08},
    "D2S1338": {
        17: 0.08,
        18: 0.05,
        19: 0.10,
        20: 0.15,
        21: 0.08,
        22: 0.07,
        23: 0.12,
        24: 0.15,
        25: 0.15,
    },
    "D19S433": {
        13: 0.15,
        14: 0.30,
        14.2: 0.05,
        15: 0.20,
        15.2: 0.05,
        16: 0.15,
        17: 0.10,
    },
    "D22S1045": {11: 0.10, 14: 0.08, 15: 0.30, 16: 0.35, 17: 0.12, 18: 0.05},
    "D10S1248": {11: 0.05, 12: 0.08, 13: 0.25, 14: 0.30, 15: 0.20, 16: 0.10, 17: 0.02},
    "D1S1656": {
        12: 0.10,
        13: 0.08,
        14: 0.05,
        15: 0.12,
        16: 0.15,
        17: 0.20,
        17.3: 0.10,
        18: 0.10,
        18.3: 0.05,
    },
    "D12S391": {
        17: 0.05,
        18: 0.15,
        19: 0.12,
        20: 0.20,
        21: 0.18,
        22: 0.15,
        23: 0.10,
        24: 0.05,
    },
    "D2S441": {
        10: 0.10,
        11: 0.20,
        11.3: 0.05,
        12: 0.08,
        13: 0.10,
        14: 0.25,
        15: 0.15,
        16: 0.07,
    },
    "SE33": {
        19: 0.05,
        20: 0.08,
        21: 0.10,
        22: 0.12,
        23: 0.10,
        24: 0.08,
        25: 0.12,
        26: 0.10,
        27: 0.10,
        28: 0.08,
        29: 0.07,
    },
}

# Normalize frequencies
for locus in LOCI:
    total = sum(ALLELE_FREQS[locus].values())
    for allele in ALLELE_FREQS[locus]:
        ALLELE_FREQS[locus][allele] /= total

# Precompute weighted choices
WEIGHTED_ALLELES = {}
for locus in LOCI:
    alleles = list(ALLELE_FREQS[locus].keys())
    weights = list(ALLELE_FREQS[locus].values())
    WEIGHTED_ALLELES[locus] = (alleles, weights)

# -------------------------------
# Helper functions
# -------------------------------


def sample_allele(locus):
    alleles, weights = WEIGHTED_ALLELES[locus]
    return random.choices(alleles, weights=weights, k=1)[0]


def format_alleles(a1, a2):
    """Format alleles according to challenge spec (e.g., '13', '13,14', '9.3')"""
    alleles = sorted([a1, a2], key=lambda x: float(x))
    if alleles[0] == alleles[1]:
        return f"{alleles[0]:.1f}".rstrip("0").rstrip(".")
    else:
        a1_str = f"{alleles[0]:.1f}".rstrip("0").rstrip(".")
        a2_str = f"{alleles[1]:.1f}".rstrip("0").rstrip(".")
        return f"{a1_str},{a2_str}"


def generate_profile(person_id):
    profile = {"PersonID": person_id}
    for locus in LOCI:
        a1 = sample_allele(locus)
        a2 = sample_allele(locus)

        # Simulate dropout
        if random.random() < DROPOUT_RATE:
            profile[locus] = "-"
            continue
        if random.random() < SINGLE_ALLELE_RATE:
            allele = random.choice([a1, a2])
            profile[locus] = f"{allele:.1f}".rstrip("0").rstrip(".")
            continue

        profile[locus] = format_alleles(a1, a2)
    return profile


def mutate_allele(allele):
    """±1 step mutation, preserves microvariant if present"""
    allele_float = float(allele)
    step = random.choice([-1, 1])
    new_val = allele_float + step
    if "." in str(allele):
        base = int(allele_float)
        micro = round((allele_float - base) * 10)
        return f"{base + step}.{micro}"
    return str(int(new_val))


def generate_child_profile(parent_profile, person_id):
    profile = {"PersonID": person_id}
    for locus in LOCI:
        parent_str = parent_profile.get(locus, "-")
        if parent_str in ("-", "") or parent_str == "-":
            profile[locus] = "-"
            continue

        # Parse parent alleles
        if "," in parent_str:
            p_alleles = [float(x) for x in parent_str.split(",")]
        else:
            p_alleles = [float(parent_str)] * 2

        # Transmit one allele (with possible mutation)
        transmitted = random.choice(p_alleles)
        if random.random() < MUTATION_RATE:
            transmitted = float(mutate_allele(str(transmitted)))

        # Other allele from population
        other = float(sample_allele(locus))

        # Apply dropout
        if random.random() < DROPOUT_RATE:
            profile[locus] = "-"
            continue
        if random.random() < SINGLE_ALLELE_RATE:
            allele = random.choice([transmitted, other])
            profile[locus] = f"{allele:.1f}".rstrip("0").rstrip(".")
            continue

        profile[locus] = format_alleles(transmitted, other)
    return profile


# -------------------------------
# Dataset generation
# -------------------------------

if __name__ == "__main__":
    print("Generating synthetic STR dataset for #codechallenge2025...")

    os.makedirs("data", exist_ok=True)

    profiles = []

    # Generate founder parents
    parent_ids = [f"P{i:06d}" for i in range(NUM_TRUE_PAIRS)]
    parent_profiles = []
    for pid in parent_ids:
        prof = generate_profile(pid)
        profiles.append(prof)
        parent_profiles.append(prof)

    # Generate their biological children
    child_profiles = []
    for i, pid in enumerate(parent_ids):
        child_id = f"C{i:06d}"
        child_prof = generate_child_profile(parent_profiles[i], child_id)
        profiles.append(child_prof)
        child_profiles.append(child_prof)

    # Fill remaining database with unrelated individuals
    remaining = NUM_DB_PROFILES - len(profiles)
    for i in range(remaining):
        pid = f"U{i + 1:06d}"
        profiles.append(generate_profile(pid))

    # Shuffle database
    random.shuffle(profiles)

    # Save database
    db_df = pd.DataFrame(profiles)
    db_df = db_df[["PersonID"] + LOCI]
    db_df.to_csv("data/str_database.csv", index=False)
    print(f"Database saved: data/str_database.csv ({len(db_df):,} profiles)")

    # Generate query profiles (35 children + 5 unrelated)
    query_profiles = []
    for i in range(NUM_TRUE_PAIRS):
        query_prof = child_profiles[i].copy()
        query_prof["PersonID"] = f"Q{i + 1:03d}"
        query_profiles.append(query_prof)

    for i in range(NUM_TRUE_PAIRS, NUM_QUERIES):
        query_profiles.append(generate_profile(f"Q{i + 1:03d}"))

    random.shuffle(query_profiles)
    query_df = pd.DataFrame(query_profiles)
    query_df = query_df[["PersonID"] + LOCI]
    query_df.to_csv("data/str_queries.csv", index=False)
    print(f"Queries saved: data/str_queries.csv ({len(query_df)} profiles)")

    # Ground truth (for validation only)
    ground_truth = []
    for i in range(NUM_TRUE_PAIRS):
        ground_truth.append(
            {"QueryID": f"Q{i + 1:03d}", "TrueCounterpartID": parent_ids[i]}
        )
    gt_df = pd.DataFrame(ground_truth)
    gt_df.to_csv("data/ground_truth.csv", index=False)
    print("Ground truth saved: data/ground_truth.csv")

    print(
        "\nDataset generation complete! Ready for the challenge on PYDay Iran, 2025 🧬"
    )
