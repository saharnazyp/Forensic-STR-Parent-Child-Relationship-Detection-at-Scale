# tests/update_leaderboard.py
"""
Reads manually reviewed all_results.csv and generates Leaderboard.md + leaderboard.json
Final score = auto_score (0–100) + review_score (0–20), capped at 120.
The 'features' column is no longer displayed in the leaderboard.
"""

import pandas as pd
import json
from datetime import datetime

CSV_INPUT = "src/codechallenge2025/results/all_results.csv"
LEADERBOARD_JSON = "leaderboard.json"
LEADERBOARD_MD = "Leaderboard.md"


def main():
    if not pd.io.common.file_exists(CSV_INPUT):
        print(f"❌ {CSV_INPUT} not found. Run `python tests/run_challenge.py` first.")
        return

    df = pd.read_csv(CSV_INPUT)

    # Ensure required columns exist
    required = [
        "user",
        "auto_score",
        "review_score",
        "accuracy",
        "runtime_sec",
        "submission_date",
    ]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing column in CSV: {col}")

    # Clean and compute scores
    df["auto_score"] = pd.to_numeric(df["auto_score"], errors="coerce").fillna(0)
    df["review_score"] = pd.to_numeric(df["review_score"], errors="coerce").fillna(0)
    df["total_score"] = (df["auto_score"] + df["review_score"]).clip(upper=120)

    # Sort by total score (desc), then runtime (asc for tie-break)
    df = df.sort_values(by=["total_score", "runtime_sec"], ascending=[False, True])
    df = df.reset_index(drop=True)
    df["rank"] = range(1, len(df) + 1)

    # Save full data to JSON
    json_records = df.to_dict(orient="records")
    with open(LEADERBOARD_JSON, "w") as f:
        json.dump(json_records, f, indent=2)

    # Generate clean Markdown leaderboard (no 'Features' column)
    md = "# 🏆 #codechallenge2025 Leaderboard\n\n"
    md += f"_Last updated: {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')}_\n\n"
    md += "| Rank | User | Final Score | Auto | Review | Accuracy | Time (s) | Date |\n"
    md += "|------|------|-------------|------|--------|----------|----------|------|\n"

    for _, row in df.iterrows():
        md += (
            f"| {row['rank']} | {row['user']} | **{row['total_score']:.1f}** | "
            f"{row['auto_score']:.1f} | {row['review_score']:.1f} | "
            f"{row['accuracy']:.1%} | {row['runtime_sec']:.2f} | "
            f"{row['submission_date']} |\n"
        )

    md += "\n---\n"
    md += "**Scoring**: Final = Auto (0–100) + Review (0–20). Maximum = 120.\n"
    md += "Auto score: based on accuracy, speed, and deadline compliance.\n"
    md += "Review score: assigned manually based on code quality and feature completeness.\n"

    with open(LEADERBOARD_MD, "w") as f:
        f.write(md)

    print(f"✅ Leaderboard updated from {CSV_INPUT}")
    print(f"   → {LEADERBOARD_MD} and {LEADERBOARD_JSON}")


if __name__ == "__main__":
    main()
