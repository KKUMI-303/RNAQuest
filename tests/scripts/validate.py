from pathlib import Path
import sys


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

    print("\nAll DESeq2 outputs generated successfully.")


if __name__ == "__main__":
    main()
