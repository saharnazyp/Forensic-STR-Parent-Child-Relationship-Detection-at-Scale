# How I Solved the PyDay 2025 DNA Matching Challenge

## The Challenge

I was given a Python coding challenge that simulates real forensic DNA analysis. The task: find parent-child relationships in a database of 5,000 DNA profiles. For each of 40 mystery people, I needed to identify their biological parent from the database.

The tricky part? The database was a mix of parents, children, and random unrelated people — all shuffled together with no labels telling me who was who.

**Note on Methodology**

To address this challenge, I utilized Claude (developed by Anthropic) as an assistive AI tool for ideation, problem-solving, and refining approaches. The complete history of my interactions with Claude is documented in the attached file `codechallenge2025/data-2025-12-18-22-30-52-batch-0000.zip` for transparency, research and reproducibility.

## Understanding the Data

Before writing any code, I spent time understanding what DNA profiles actually look like. Each person's profile is essentially a form with 21 questions (called "loci" in biology). At each question, every person has two answers (called "alleles") — one inherited from their mother, one from their father.

Here's what the data looked like:

```
PersonID  | TH01    | D8S1179 | vWA     | ...
----------|---------|---------|---------|----
P000001   | 8,9     | 13,16   | 14,16   | ...
C000001   | 9,14    | 13,15   | 16,17   | ...
```

The fundamental rule I learned: **a child must share exactly one allele with their parent at each locus**. So if a child has `9,14` at TH01, one of those numbers came from mom and one from dad.

## My First Attempts (And Why They Failed)

My initial approach was simple: count how many loci share at least one allele. More shared alleles = more likely to be related. This gave me only 28% accuracy.

I then tried using statistical likelihood ratios, weighting rare alleles more heavily. This crashed to 0% accuracy — a complete disaster.

The problem? I was missing something fundamental about how the test data was structured.

## The Breakthrough

After analyzing the dataset generator code, I discovered the hidden structure:

- **P000xxx**: Parents (35 of them)
- **C000xxx**: Children of those parents (also in the database!)
- **Q001-Q035**: Query profiles that are **copies** of the children
- **U000xxx**: Random unrelated people

This meant when I searched for Q001's parent, I was finding C000001 first (because it's literally the same person with a different ID), not P000001 (the actual parent).

## The Key Insight

I realized there's a critical difference between "same person" and "parent-child":

| Scenario | What Happens |
|----------|--------------|
| Same person | The allele strings are **identical** (e.g., both have `"9,14"`) |
| Parent-child | They share ONE allele but strings are **different** (e.g., `"9,14"` vs `"8,9"`) |

I called these "identical" matches versus "partial" matches. A same-person comparison would have almost all identical matches. A parent-child comparison would have many partial matches — because while they share one allele per locus, the child's other allele comes from the unknown second parent.

## The Final Solution

My algorithm works in four steps:

**Step 1: Compare each locus**

For every question in the DNA profile, I classify the comparison into one of five categories:
- **Identical**: Same string exactly
- **Partial**: Different strings but share at least one allele
- **Mismatch**: No shared alleles at all
- **Mutation**: No exact match but off by ±1 (rare genetic change)
- **Missing**: One or both profiles have no data

**Step 2: Filter out non-candidates**

I remove anyone with more than 2 mismatches — unrelated people typically have 5-10 mismatches, so this filter eliminates most of the database quickly.

**Step 3: Distinguish parent from same-person**

Here's the critical filter: I require at least 3 "partial" matches. Same-person comparisons have almost zero partial matches (everything is identical). Parent-child comparisons have 10-15 partial matches. This single rule solved the core problem.

**Step 4: Score and rank**

For remaining candidates, I calculate a score based on:
- The rarity of shared alleles (sharing a rare allele like `TH01=11` which only 2% of people have is stronger evidence than sharing `TPOX=8` which 40% have)
- A bonus for the number of partial matches

The candidate with the highest score becomes my #1 prediction.

## The Results

After implementing this approach, I've tested the function over 16 runs and achieved the following results:

| Metric             | Highest Score          | On Average             | Lowest Score           |
|--------------------|------------------------|------------------------|------------------------|
| **Accuracy**       | 100.0%                 | ~90.71%                | 80.0%                  |
| **Execution time** | 65.71 seconds per run  | ~61.52 seconds per run | 58.54 seconds per run  |
| **Final score**    | 120.0 out of 120       | ~110.71 out of 120     | 100.0 out of 120       |

The missed cases are primarily edge cases where random generation produced ambiguous data. Given the use of randomly generated datasets, some variance in results is unavoidable.

## What I Learned

This challenge taught me that understanding the problem structure matters more than fancy algorithms. My statistical likelihood ratio approach was mathematically sound but failed completely because I hadn't understood the test data structure.

The winning insight wasn't about biology or statistics — it was recognizing that "identical strings" versus "partial matches" could distinguish same-person from parent-child relationships. Once I saw that pattern, the solution became straightforward.

The allele frequency data I used isn't cheating — it's publicly available in the challenge repository and mirrors how real forensic labs operate. They use published population statistics to assess how meaningful a genetic match is.

## Summary

The core algorithm in plain terms:

1. For each mystery person, compare their DNA against everyone in the database
2. At each of 21 genetic markers, check if they share any values
3. Filter out people who share too few values (strangers) or share values too perfectly (same person)
4. Rank remaining candidates by how rare their shared values are
5. Return the top match as the predicted parent