import argparse
import pandas as pd
from pathlib import Path
import re
import fnmatch

DATA_DIR = Path(__file__).resolve().parent / "GDC_download"


def normalize_tcga_id(identifier: str) -> str:
    """Normalize various TCGA barcode formats to TCGA-XX-YYYY."""
    m = re.match(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", identifier)
    if m:
        return m.group(0)
    if identifier.startswith("TCGA-"):
        parts = identifier.split("-")
        if len(parts) >= 3:
            return f"TCGA-{parts[1]}-{parts[2]}"
    return identifier


def extract_sample_id(path: Path) -> str:
    """Extract sample or aliquot identifier from a file path."""
    name = path.name
    m = re.search(r"(TCGA-[A-Z0-9-]+)", name)
    if m:
        return normalize_tcga_id(m.group(1))
    m = re.search(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", name)
    if m:
        return m.group(0)
    return normalize_tcga_id(path.stem)


def parse_mirna_quant(path: Path) -> pd.DataFrame:
    """Parse miRNA quantification file into wide format."""
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(path, sep=sep)
    df = df[['miRNA_ID', 'reads_per_million_miRNA_mapped']]
    wide = df.set_index('miRNA_ID').T
    wide.index = [extract_sample_id(path)]
    return wide


def parse_rppa(path: Path) -> pd.DataFrame:
    """Parse RPPA data file into wide format."""
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(path, sep=sep)
    df = df[['peptide_target', 'protein_expression']]
    wide = df.set_index('peptide_target').T
    wide.index = [extract_sample_id(path)]
    return wide


def parse_seg(path: Path) -> pd.DataFrame | None:
    """Parse segmentation file and summarize by chromosome."""
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(path, sep=sep)
    if 'Segment_Mean' in df.columns:
        agg = df.groupby('Chromosome')['Segment_Mean'].mean()
    elif 'Copy_Number' in df.columns:
        agg = df.groupby('Chromosome')['Copy_Number'].mean()
    else:
        return None
    result = agg.to_frame().T
    result.index = [extract_sample_id(path)]
    return result


def collect_data(data_dir: Path = DATA_DIR, pattern: str = "*") -> pd.DataFrame:
    """Iterate through the download directory and consolidate data."""
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    files = [f for f in data_dir.rglob("*") if fnmatch.fnmatch(f.name, pattern) and f.suffix in {".txt", ".tsv", ".csv"}]
    print("Matched files:")
    for f in files:
        print(f" - {f}")
    if not files:
        raise FileNotFoundError(f"No files matched pattern '{pattern}' in {data_dir}")

    mirna_frames: list[pd.DataFrame] = []
    rppa_frames: list[pd.DataFrame] = []
    seg_frames: list[pd.DataFrame] = []

    for f in files:
        name = f.name.lower()
        try:
            if "quantification" in name and "mirbase21" in name:
                mirna_frames.append(parse_mirna_quant(f))
            elif "rppa" in name:
                rppa_frames.append(parse_rppa(f))
            elif "seg" in name:
                parsed = parse_seg(f)
                if parsed is not None:
                    seg_frames.append(parsed)
        except Exception as e:
            print(f"Failed to parse {f}: {e}")

    combined: list[pd.DataFrame] = []
    if mirna_frames:
        all_columns = []
        for df in mirna_frames:
            all_columns.extend(df.columns.tolist())
        duplicates = {x for x in all_columns if all_columns.count(x) > 1}
        if duplicates:
            print("Duplicate miRNA_IDs in columns:", duplicates)
        mirna_concat = pd.concat(mirna_frames, axis=0, sort=False)
        mirna_concat = mirna_concat[~mirna_concat.index.duplicated(keep="first")]
        if mirna_concat.empty:
            raise ValueError("miRNA merge produced 0 rows")
        combined.append(mirna_concat)
    if rppa_frames:
        all_columns = []
        for df in rppa_frames:
            all_columns.extend(df.columns.tolist())
        duplicates = {x for x in all_columns if all_columns.count(x) > 1}
        if duplicates:
            print("Duplicate rppa_IDs in columns:", duplicates)
        rppa_concat = pd.concat(rppa_frames, axis=0, sort=False)
        rppa_concat = rppa_concat[~rppa_concat.index.duplicated(keep="first")]
        if rppa_concat.empty:
            raise ValueError("RPPA merge produced 0 rows")
        combined.append(rppa_concat)
    if seg_frames:
        all_columns = []
        for df in seg_frames:
            all_columns.extend(df.columns.tolist())
        duplicates = {x for x in all_columns if all_columns.count(x) > 1}
        if duplicates:
            print("Duplicate seg_IDs in columns:", duplicates)
        seg_concat = pd.concat(seg_frames, axis=0, sort=False)
        seg_concat = seg_concat[~seg_concat.index.duplicated(keep="first")]
        if seg_concat.empty:
            raise ValueError("Segmentation merge produced 0 rows")
        combined.append(seg_concat)

    if not combined:
        raise ValueError("No data parsed from matched files")

    final = pd.concat(combined, axis=1, sort=False).fillna(0)
    final.index = [normalize_tcga_id(idx) for idx in final.index]
    final.index.name = "sample_id"
    return final


def main():
    parser = argparse.ArgumentParser(description="Consolidate CRC data")
    parser.add_argument("--input", type=Path, default=DATA_DIR, help="Input directory")
    parser.add_argument("--pattern", default="*", help="Glob pattern for files")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "crc_consolidated.csv",
        help="Output CSV path",
    )
    args = parser.parse_args()

    data = collect_data(args.input, args.pattern)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_path = args.out.resolve()
    data.to_csv(out_path)
    print(f"Final shape: {data.shape}")
    print(f"Consolidated data written to {out_path}")


if __name__ == '__main__':
    main()
