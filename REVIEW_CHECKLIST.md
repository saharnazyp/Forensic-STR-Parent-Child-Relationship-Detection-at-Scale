# 🧪 Code Review Checklist – #codechallenge2025  
**Total Review Score: 0–20 points**  
*Assigned by organizer after inspecting participant’s `match_single()` implementation.*

**Challenge**  
Develop a highly scalable Python program that identifies likely single parent-child relationships (father or mother only) within a large mixed dataset of DNA profiles.

The dataset consists of approximately 500,000 known STR profiles (parents and children mixed together, roles unknown) plus approximately 40 unknown query profiles. For each query profile, the program must search the large database and return a ranked list of the top candidate matches that could be either the parent of the query individual or the child of the query individual (bidirectional relationship detection).

**Key Biological Concepts (Integrated)**  
- **STR (Short Tandem Repeat)**: Highly variable regions in human DNA where a short sequence (typically 2–6 bases) is repeated a variable number of times. These are the primary markers used in forensic and paternity DNA testing.  
- **Locus (plural: loci)**: A specific physical location on a chromosome containing an STR marker. Standard forensic panels use 13–24 loci (e.g., D8S1179, TH01, vWA, FGA, D21S11, etc.).  
- **Allele**: The number of repeats observed at a locus (e.g., 13, 14, or 14.2 for microvariants). Each person has two alleles per locus – one inherited from each biological parent.  
- **Homozygous**: Both alleles at a locus are identical (e.g., 13,13). **Heterozygous**: The two alleles differ (e.g., 13,14).  
- **Genotype / STR Profile**: The complete set of allele pairs across all tested loci for an individual.  
- **Mendelian Inheritance**: A child must receive exactly one allele at each locus from the biological mother and one from the biological father. In single-parent testing, at least one of the child’s alleles must match one of the candidate parent’s alleles at every locus (unless explained by mutation).  
- **Mutation**: A rare change in repeat count during meiosis (typical rate 0.001–0.004 per locus per generation), usually differing by ±1 repeat. Mutations allow small discrepancies without excluding a true relationship.  
- **Missing Data / Allele Dropout / Partial Profiles**: Common in degraded or low-quantity DNA samples; some loci may have no result ("-") or only one allele may amplify.  
- **Likelihood Ratio (LR) / Paternity Index**: The statistical measure used to evaluate evidence. For each locus, LR = (probability of observed alleles if parent-child) / (probability of observed alleles if unrelated). The Combined Likelihood Ratio (CLR or CPI) is the product of individual locus LRs. Typical interpretation: CLR > 10,000 indicates strong support; > 100,000 is considered virtually proof.  
- **Allele Frequencies**: Population-specific frequencies required to calculate the “unrelated” probability. Rare alleles shared between individuals dramatically increase the LR.

**Input Format**  
CSV file with columns:  
- `PersonID` (unique identifier)  
- One column per locus (e.g., `D8S1179`, `TH01`, `vWA`, etc.)  
Allele values formatted as:  
- `"13,14"` (heterozygous)  
- `"13"` or `"13,13"` (homozygous)  
- `"-"` or blank (missing locus/data)

**Requirements**  
- Compute accurate per-locus likelihood ratios handling exact matches, ±1 step mutations, allele dropout, and complete mismatches.  
- Calculate the overall Combined Likelihood Ratio (CLR) for each candidate pair.  
- Support bidirectional matching (same logic works for “query is child” or “query is parent”).  
- Handle partial profiles and missing data robustly.  
- Scale to ~500,000 database profiles × ~40 queries (millions of comparisons) – naive pairwise comparison is too slow. Implement efficient candidate pre-filtering/indexing (e.g., using shared rare alleles, multi-locus signatures, or hashing) to achieve practical runtime.  
- For each query, output the top 10 candidates ranked by CLR (descending), including: PersonID, CLR value, estimated posterior probability (using a reasonable prior, e.g., 50%), number of consistent loci, number of loci requiring mutation, and number of missing/inconclusive loci.  

**Bonus features:**  
- Support microvariants (e.g., 9.3)  
- Tri-allelic loci, null alleles  
- Visualize matches or export detailed reports

This challenge mirrors real-world forensic bioinformatics and paternity testing systems used in laboratories and courts worldwide. Good luck uncovering the hidden family connections! 🧬

> ⚠️ **Goal**: Assess whether the solution **correctly and robustly** handles real-world forensic STR challenges—not just passes the current dataset.


---

## 📋 Evaluation Criteria

#| # | Feature / Requirement | Max Points | How to Check |
#|---|------------------------|------------|--------------|
#| 1 | **Microvariant support**<br>(e.g., `9.3`, `30.2`, `17.3`) | 4 | Does code parse alleles as **floats** (not strings)?<br>Does it handle comparisons like `9.3 == 9.3` and mutation `9.3 → 10.3`? |
#| 2 | **Mutation tolerance**<br>(±1 repeat per locus) | 4 | Does matching allow **1-step difference** when evaluating parent-child consistency?<br>Is it applied **per locus**, not globally? |
#| 3 | **Missing data handling**<br>(dropout `"-"`, single-allele) | 4 | Does it skip or probabilistically handle:<br>– `"-"` (no data)<br>– `"13"` (only one allele)?<br>Does it avoid crashing or false mismatches? |
#| 4 | **Bidirectional matching**<br>(query = parent **or** child) | 3 | Is the logic **symmetric**?<br>Does it work whether the DB profile is parent or child of the query? |
#| 5 | **Statistical scoring (LR/CLR)**<br>Uses allele frequencies | 3 | Does it compute **Likelihood Ratio** per locus using population frequencies?<br>Or does it just count matching loci (❌ insufficient)? |
#| 6 | **Efficiency / Scalability**<br>(~500k × 40 ≈ 20M comparisons) | 2 | Does it use **indexing, filtering, or hashing**?<br>Or is it brute-force O(N) per query (❌ too slow at scale)? |


> ✅ **Full points** require **correct implementation**, not just presence of keywords.  
> ❌ **0 points** for unhandled edge cases that break correctness.


---

## 📝 Scoring Guide

- **18–20**: Excellent — handles all cases robustly, clean code, efficient, proper LR.
- **14–17**: Good — minor gaps (e.g., no microvariants, or hard-coded loci).
- **10–13**: Fair — works on basic data but fails on edge cases (mutations/dropout).
- **0–9**: Poor — crashes, ignores key biology, or uses naive matching.


> 💡 Tip: Run the solution mentally (or with debug prints) on a profile with `"TH01": "9.3"` and a mutation like `13 → 14`.


---

## 📤 Output Fields (for `all_results.csv`)

After review, fill these columns:

- `review_score`: Numeric (0–20)
- `features`: Short summary, e.g.,  
  `"✅ Micro, ✅ Mut, ❌ LR"` or `"✅ All"`
- `comments`: Optional notes, e.g.,  
  `"Used string comparison — fails on 9.3"`


---

> ℹ️ **Note**: The synthetic dataset **includes microvariants, mutations, dropout, and partial profiles** — a robust solution **must** handle them to be forensically valid.
