#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd


def parse_file_mapping(items: List[str]) -> Dict[int, Path]:

    mapping: Dict[int, Path] = {}

    for item in items:
        if ":" not in item:
            raise ValueError(
                f"Invalid --files entry '{item}'. Expected format MODEL:filename"
            )
        model_str, file_str = item.split(":", 1)

        try:
            model_nr = int(model_str)
        except ValueError as exc:
            raise ValueError(
                f"Invalid model number '{model_str}' in --files entry '{item}'"
            ) from exc

        if model_nr in mapping:
            raise ValueError(f"Duplicate model number in --files: {model_nr}")

        mapping[model_nr] = Path(file_str)

    if not mapping:
        raise ValueError("No input files were provided.")

    return mapping


def resolve_columns(df: pd.DataFrame) -> Dict[str, str]:
    return {col.lower(): col for col in df.columns}


def standardize_subset(
    df: pd.DataFrame,
    snp_cols: List[str],
    position_col: str,
    pvalue_col: str,
    model_nr: int,
    file_label: str,
) -> pd.DataFrame:

    cols = resolve_columns(df)

    required = [*snp_cols, position_col, pvalue_col]
    missing = [c for c in required if c.lower() not in cols]
    if missing:
        raise ValueError(f"{file_label} is missing columns: {missing}")

    selected_cols = [cols[c.lower()] for c in snp_cols]
    selected_cols += [cols[position_col.lower()], cols[pvalue_col.lower()]]

    sub = df[selected_cols].copy()

    standardized_names = list(snp_cols) + ["position", "pvalue"]
    sub.columns = standardized_names

    for col in snp_cols:
        sub[col] = sub[col].astype(str).str.strip()

    sub["position"] = pd.to_numeric(sub["position"], errors="coerce")
    sub["pvalue"] = pd.to_numeric(sub["pvalue"], errors="coerce")

    sub = sub.dropna(subset=["position", "pvalue"])

    for col in snp_cols:
        sub = sub[sub[col].notna()]
        sub = sub[sub[col].astype(str).str.strip() != ""]

    sub["model"] = model_nr
    return sub


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Find overlapping SNPs across association result models."
    )
    parser.add_argument(
        "--analysis-name",
        required=True,
        help="Short label used in output filenames, e.g. logistic, snptest, imputed_snptest",
    )
    parser.add_argument(
        "--base-dir",
        default="data",
        help="Directory containing input files (default: data)",
    )
    parser.add_argument(
        "--out-dir",
        default="results",
        help="Directory where outputs will be written (default: results)",
    )
    parser.add_argument(
        "--files",
        nargs="+",
        required=True,
        help="Model:file mappings, e.g. 1:file1.txt 2:file2.txt",
    )
    parser.add_argument(
        "--snp-cols",
        nargs="+",
        required=True,
        help="One or more SNP identity columns used for grouping, e.g. snp OR alternate_ids rsid",
    )
    parser.add_argument(
        "--position-col",
        required=True,
        help="Column containing genomic position, e.g. bp or position",
    )
    parser.add_argument(
        "--pvalue-col",
        required=True,
        help="Column containing p-values, e.g. p or frequentist_add_pvalue",
    )

    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    file_map = parse_file_mapping(args.files)

    all_rows = []

    for model_nr, rel_path in sorted(file_map.items()):
        file_path = rel_path if rel_path.is_absolute() else base_dir / rel_path

        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        df = pd.read_csv(file_path, sep=r"\s+", engine="python")

        if df.empty:
            print(f"Skipping empty file: {file_path}")
            continue

        sub = standardize_subset(
            df=df,
            snp_cols=args.snp_cols,
            position_col=args.position_col,
            pvalue_col=args.pvalue_col,
            model_nr=model_nr,
            file_label=str(file_path),
        )

        if sub.empty:
            print(f"Skipping file with no valid rows after cleaning: {file_path}")
            continue

        all_rows.append(sub)

    if not all_rows:
        raise ValueError("No non-empty valid result files were found.")

    combined = pd.concat(all_rows, ignore_index=True)

    sort_cols = list(args.snp_cols) + ["position", "model", "pvalue"]
    combined = combined.sort_values(sort_cols).copy()

    combined_csv = out_dir / f"all_{args.analysis_name}_snps.csv"
    combined.to_csv(combined_csv, index=False)

    group_cols = list(args.snp_cols) + ["position"]

    summary = (
        combined.groupby(group_cols, as_index=False)
        .agg(
            models_identified=("model", lambda x: ", ".join(map(str, sorted(set(x))))),
            model_count=("model", lambda x: len(set(x))),
            min_pvalue=("pvalue", "min"),
        )
    )

    summary = summary.sort_values(
        by=["model_count", "min_pvalue"] + list(args.snp_cols),
        ascending=[False, True] + [True] * len(args.snp_cols),
    ).copy()

    summary_csv = out_dir / f"{args.analysis_name}_snp_intersections_summary.csv"
    summary.to_csv(summary_csv, index=False)

    print("Done.")
    print(f"Saved merged SNP list to: {combined_csv}")
    print(f"Saved grouped summary to: {summary_csv}")
    print()
    print(summary.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
