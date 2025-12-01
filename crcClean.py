import argparse
import re
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).resolve().parent / "GDC_download"


def extract_sample_id(path: Path) -> str:
    """Extract sample or aliquot identifier from a file path."""

    name = path.name
    m = re.search(r"(TCGA-[A-Z0-9-]+)", name)
    if m:
        return m.group(1)
    m = re.search(
        r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", name
    )
    if m:
        return m.group(0)
    return path.stem


def parse_mirna_quant(path: Path) -> pd.DataFrame:
    """Parse miRNA quantification file into wide format.

    Multiple rows can exist per miRNA (e.g., isoforms). Values are summed so
    that each miRNA appears once per sample.
    """

    df = pd.read_csv(path, sep="\t")
    rpm_column_options = [
        "reads_per_million_miRNA_mapped",
        "reads_per_million_mirna_mapped",
        "RPM",
    ]
    rpm_column = next((c for c in rpm_column_options if c in df.columns), None)
    if rpm_column is None:
        raise ValueError(f"RPM column not found in {path.name}")

    df = df[["miRNA_ID", rpm_column]].rename(columns={rpm_column: "RPM"})
    collapsed = df.groupby("miRNA_ID", as_index=False)["RPM"].sum()
    wide = collapsed.set_index("miRNA_ID").T
    wide.index = [extract_sample_id(path)]
    wide.columns = [f"mirna__{c}" for c in wide.columns]
    return wide


def parse_rppa(path: Path) -> pd.DataFrame:
    """Parse RPPA data file into wide format."""

    df = pd.read_csv(path, sep="\t")
    df = df[["peptide_target", "protein_expression"]]
    wide = df.set_index("peptide_target").T
    wide.index = [extract_sample_id(path)]
    wide.columns = [f"rppa__{c}" for c in wide.columns]
    return wide


def parse_seg(path: Path) -> pd.DataFrame | None:
    """Parse segmentation file and summarize by chromosome."""

    df = pd.read_csv(path, sep="\t")
    if "Chromosome" not in df.columns:
        return None

    chromosome = df["Chromosome"].astype(str).str.replace("chr", "", regex=False)
    df = df.assign(Chromosome=chromosome)

    if "Segment_Mean" in df.columns:
        agg = df.groupby("Chromosome")["Segment_Mean"].mean()
    elif "Copy_Number" in df.columns:
        agg = df.groupby("Chromosome")["Copy_Number"].mean()
    else:
        return None

    agg.index = [f"seg__chr{c}" for c in agg.index]
    result = agg.to_frame().T
    result.index = [extract_sample_id(path)]
    return result


def collect_data(data_dir: Path = DATA_DIR) -> pd.DataFrame:
    """Iterate through the download directory and consolidate data."""

    data_dir = data_dir.expanduser().resolve()
    if not data_dir.exists():
        print(f"Data directory not found: {data_dir}")
        return pd.DataFrame()

    mirna_frames: list[pd.DataFrame] = []
    rppa_frames: list[pd.DataFrame] = []
    seg_frames: list[pd.DataFrame] = []

    for f in data_dir.rglob("*"):
        if f.suffix.lower() not in {".txt", ".tsv"}:
            continue
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

    combined = []
    if mirna_frames:
        mirna_concat = pd.concat(mirna_frames, axis=0, sort=False)
        mirna_concat = mirna_concat[~mirna_concat.index.duplicated(keep="first")]
        combined.append(mirna_concat)
        print(f"Loaded miRNA samples: {len(mirna_concat)}")

    if rppa_frames:
        rppa_concat = pd.concat(rppa_frames, axis=0, sort=False)
        rppa_concat = rppa_concat[~rppa_concat.index.duplicated(keep="first")]
        combined.append(rppa_concat)
        print(f"Loaded RPPA samples: {len(rppa_concat)}")

    if seg_frames:
        seg_concat = pd.concat(seg_frames, axis=0, sort=False)
        seg_concat = seg_concat[~seg_concat.index.duplicated(keep="first")]
        combined.append(seg_concat)
        print(f"Loaded segmentation samples: {len(seg_concat)}")

    if not combined:
        return pd.DataFrame()

    final = pd.concat(combined, axis=1, sort=False).fillna(0)
    final.index.name = "sample_id"
    return final


def main():
    parser = argparse.ArgumentParser(description="Consolidate CRC data")
    parser.add_argument(
        "--data-dir", type=Path, default=DATA_DIR,
        help="Directory containing downloaded data",
    )
    parser.add_argument(
        "--output", type=Path,
        default=Path(__file__).resolve().parent / 'crc_consolidated.csv',
        help="Path for the consolidated CSV output",
    )
    args = parser.parse_args()

    data = collect_data(args.data_dir)
    if not data.empty:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        data.to_csv(args.output)
        print(f"Consolidated data written to {args.output}")
    else:
        print('No data parsed. Nothing was written.')


if __name__ == '__main__':
    main()
