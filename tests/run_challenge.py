# tests/run_challenge.py (final version for batch run)
"""
Batch-evaluates all submissions and outputs all_results.csv with review columns (initially empty).
"""

import os
import time
import pandas as pd
import importlib.util
import sys
from pathlib import Path
from datetime import datetime
import re

DEADLINE = datetime(2025, 12, 18)
SUBMISSIONS_DIR = Path("src/codechallenge2025/submissions")
DB_PATH = "data/str_database.csv"
QUERIES_PATH = "data/str_queries.csv"
GT_PATH = "data/ground_truth.csv"


def ensure_dataset():
    if not (os.path.exists(DB_PATH) and os.path.exists(QUERIES_PATH)):
        print("Dataset not found. Generating...")
        os.system("uv run src/codechallenge2025/dataset_generator.py")


def load_module_from_path(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[path.stem] = module
    spec.loader.exec_module(module)
    return module


def parse_submission_date(filename: str):
    match = re.search(r"_(\d{8})\.py$", filename)
    if not match:
        return None
    try:
        return datetime.strptime(match.group(1), "%Y%m%d")
    except ValueError:
        return None


def evaluate_results(participant_results: list) -> float:
    if not os.path.exists(GT_PATH):
        return 0.0
    gt = pd.read_csv(GT_PATH)
    gt_dict = dict(zip(gt["QueryID"], gt["TrueCounterpartID"]))
    correct = 0
    total = len(gt_dict)
    for result in participant_results:
        qid = result["query_id"]
        if qid not in gt_dict:
            continue
        true_id = gt_dict[qid]
        top = result.get("top_candidates", [])
        if top and top[0]["person_id"] == true_id:
            correct += 1
    return correct / total if total > 0 else 0.0


def score_auto(accuracy: float, runtime: float, is_late: bool) -> float:
    if accuracy <= 0.0:
        return 0.0

    base = accuracy * 100
    speed_bonus = 0
    late_penalty = 10 if is_late else 0

    if runtime < 300:
        speed_bonus = 20
    elif runtime < 600:
        speed_bonus = 10

    final = base + speed_bonus - late_penalty
    return max(0.0, min(100.0, final))


def process_submission(sub_file: Path):
    print(f"→ {sub_file.name}")
    sub_date = parse_submission_date(sub_file.name)
    is_late = sub_date.date() > DEADLINE.date() if sub_date else True

    try:
        mod = load_module_from_path(sub_file)
        start = time.time()
        results = mod.find_matches(DB_PATH, QUERIES_PATH)
        runtime = time.time() - start
        accuracy = evaluate_results(results)
        auto_score = score_auto(accuracy, runtime, is_late)
    except Exception as e:
        print(f"  ❌ Error: {e}")
        return {
            "filename": sub_file.name,
            "user": sub_file.stem.split("_")[0].replace("-", " ").title(),
            "submission_date": sub_date.strftime("%Y-%m-%d") if sub_date else "unknown",
            "is_late": is_late,
            "accuracy": 0.0,
            "runtime_sec": 0,
            "auto_score": 0.0,
            "review_score": "",  # ← empty for manual fill
            "features": "",  # ← e.g., "✅ Mutations, ✅ Microvariants"
            "comments": "",  # ← optional notes
            "error": str(e),
        }

    return {
        "filename": sub_file.name,
        "user": sub_file.stem.split("_")[0].replace("-", " ").title(),
        "submission_date": sub_date.strftime("%Y-%m-%d") if sub_date else "unknown",
        "is_late": is_late,
        "accuracy": round(accuracy, 6),
        "runtime_sec": round(runtime, 2),
        "auto_score": round(auto_score, 1),
        "review_score": "",  # ← TO BE FILLED BY REVIEWER
        "features": "",  # ← TO BE FILLED
        "comments": "",  # ← optional
        "error": "",
    }


def main():
    ensure_dataset()
    subs = sorted([f for f in SUBMISSIONS_DIR.glob("*.py") if f.name != "__init__.py"])
    if not subs:
        print("No submissions found.")
        return

    print(f"Processing {len(subs)} submissions...")
    results = [process_submission(f) for f in subs]

    df = pd.DataFrame(results)
    # Reorder columns for clarity
    cols = [
        "filename",
        "user",
        "submission_date",
        "is_late",
        "accuracy",
        "runtime_sec",
        "auto_score",
        "review_score",
        "features",
        "comments",
        "error",
    ]
    df = df[cols]
    df.to_csv("src/codechallenge2025/results/all_results.csv", index=False)
    print("\n✅ Saved to all_results.csv")
    print(
        "👉 Please fill in 'review_score' (0–20), 'features', and 'comments' manually."
    )


if __name__ == "__main__":
    main()
