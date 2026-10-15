from pathlib import Path
import sys

import pandas as pd


RESULTS = Path("tests/results/differential_expression")

EXPECTED_FILES = [
    "DE_results.csv",
    "comparison.csv",
    "PCA_conditions.pdf",
    "volcano.pdf",
    "heatmap.pdf",
]


def main():
    failed = False

    # Check expected output files.
    for filename in EXPECTED_FILES:
        path = RESULTS / filename

        if not path.exists():
            print(f"FAIL: {path} does not exist")
            failed = True
        elif path.stat().st_size == 0:
            print(f"FAIL: {path} is empty")
            failed = True
        else:
            print(f"PASS: {path}")

    if failed:
        sys.exit(1)

    # Summarize differential-expression results.
    path = RESULTS / "DE_results.csv"
    df = pd.read_csv(path)

    required_columns = {"gene_id", "log2FoldChange", "padj"}

    if not required_columns.issubset(df.columns):
        print("FAIL: DE_results.csv is missing required columns")
        sys.exit(1)

    sig = df[df["padj"] < 0.05].copy()

    higher_in_b = (sig["log2FoldChange"] > 0).sum()
    higher_in_a = (sig["log2FoldChange"] < 0).sum()

    print(f"\nGenes tested: {len(df)}")
    print(f"Significant genes (padj < 0.05): {len(sig)}")
    print(f"Higher in Condition B: {higher_in_b}")
    print(f"Higher in Condition A: {higher_in_a}")

    print("\nSignificant genes:")
    if sig.empty:
        print("No significant genes detected.")
    else:
        print(
            sig[["gene_id", "log2FoldChange", "padj"]]
            .sort_values("padj")
            .to_string(index=False)
        )

    # Basic biological sanity check.
    if sig.empty:
        print("\nFAIL: No significant genes detected.")
        sys.exit(1)

    print("\nAll pipeline output checks passed.")
    print("Differential-expression summary completed.")


if __name__ == "__main__":
    main()
