# CRC_Nutrigenomics_Summer_25

This repository contains scripts for handling multi-omic data from colon cancer studies.

- `crcClean.py` consolidates miRNA, RPPA, and copy number data into a single CSV file.
- `qvae.py` provides a Quantum Variational Autoencoder implementation using PyTorch and PennyLane.

To run the QVAE demo with synthetic data:

```bash
python qvae.py --epochs 5
```

If you have generated `crc_consolidated.csv` with `crcClean.py`, place it in the
repository root and run:

```bash
python qvae.py --csv crc_consolidated.csv --epochs 20
```
