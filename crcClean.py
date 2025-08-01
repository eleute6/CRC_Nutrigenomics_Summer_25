import pandas as pd
from pathlib import Path
import re

DATA_DIR = Path('GDC_download')


def extract_sample_id(path: Path) -> str:
    """Extract sample or aliquot identifier from a file path."""
    name = path.name
    m = re.search(r"(TCGA-[A-Z0-9-]+)", name)
    if m:
        return m.group(1)
    m = re.search(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", name)
    if m:
        return m.group(0)
    return path.stem


def parse_mirna_quant(path: Path) -> pd.DataFrame:
    """Parse miRNA quantification file into wide format."""
    df = pd.read_csv(path, sep='\t')
    df = df[['miRNA_ID', 'reads_per_million_miRNA_mapped']]
    wide = df.set_index('miRNA_ID').T
    wide.index = [extract_sample_id(path)]
    return wide


def parse_rppa(path: Path) -> pd.DataFrame:
    """Parse RPPA data file into wide format."""
    df = pd.read_csv(path, sep='\t')
    df = df[['peptide_target', 'protein_expression']]
    wide = df.set_index('peptide_target').T
    wide.index = [extract_sample_id(path)]
    return wide


def parse_seg(path: Path) -> pd.DataFrame | None:
    """Parse segmentation file and summarize by chromosome."""
    df = pd.read_csv(path, sep='\t')
    if 'Segment_Mean' in df.columns:
        agg = df.groupby('Chromosome')['Segment_Mean'].mean()
    elif 'Copy_Number' in df.columns:
        agg = df.groupby('Chromosome')['Copy_Number'].mean()
    else:
        return None
    result = agg.to_frame().T
    result.index = [extract_sample_id(path)]
    return result


def collect_data() -> pd.DataFrame:
    """Iterate through the download directory and consolidate data."""
    mirna_frames = []
    rppa_frames = []
    seg_frames = []

    for f in DATA_DIR.rglob('*'):
        if f.suffix not in {'.txt', '.tsv'}:
            continue
        name = f.name.lower()
        try:
            if 'quantification' in name and 'mirbase21' in name:
                mirna_frames.append(parse_mirna_quant(f))
            elif 'rppa' in name:
                rppa_frames.append(parse_rppa(f))
            elif 'seg' in name:
                parsed = parse_seg(f)
                if parsed is not None:
                    seg_frames.append(parsed)
        except Exception as e:
            print(f"Failed to parse {f}: {e}")

    combined = []
    if mirna_frames:
        combined.append(pd.concat([df.reset_index(drop=True) for df in mirna_frames], sort=False))
    if rppa_frames:
        combined.append(pd.concat([df.reset_index(drop=True) for df in rppa_frames], sort=False))
    if seg_frames:
        combined.append(pd.concat([df.reset_index(drop=True) for df in seg_frames], sort=False))

    if not combined:
        return pd.DataFrame()
    final = pd.concat(combined, axis=1, sort=False).fillna(0)
    final.index.name = 'sample_id'
    return final


def main():
    data = collect_data()
    if not data.empty:
        data.to_csv('crc_consolidated.csv')
        print('Consolidated data written to crc_consolidated.csv')
    else:
        print('No data parsed. Nothing was written.')


if __name__ == '__main__':
    main()
