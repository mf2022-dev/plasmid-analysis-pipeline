#!/usr/bin/env python3
"""
Validate and parse the input samplesheet for PlasmidScope pipeline.

Samplesheet format (CSV):
    sample,short_reads_1,short_reads_2,long_reads,assembly
    sample1,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,/path/to/long.fastq.gz,
    sample2,,,,/path/to/assembly.fasta
"""

import sys
import csv
import os
from pathlib import Path


def validate_samplesheet(samplesheet_path):
    """Validate and parse the input samplesheet."""

    required_columns = ['sample']
    optional_columns = ['short_reads_1', 'short_reads_2', 'long_reads', 'assembly']
    valid_extensions_fastq = ['.fastq', '.fastq.gz', '.fq', '.fq.gz']
    valid_extensions_fasta = ['.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz']

    samples = []
    errors = []

    with open(samplesheet_path, 'r') as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames

        # Check required columns
        for col in required_columns:
            if col not in header:
                errors.append(f"ERROR: Missing required column '{col}' in samplesheet")

        for row_idx, row in enumerate(reader, start=2):
            sample_id = row.get('sample', '').strip()

            if not sample_id:
                errors.append(f"ERROR: Row {row_idx}: Empty sample name")
                continue

            # Validate sample name (alphanumeric, underscore, hyphen only)
            if not all(c.isalnum() or c in '_-.' for c in sample_id):
                errors.append(
                    f"ERROR: Row {row_idx}: Sample name '{sample_id}' contains "
                    f"invalid characters. Use only alphanumeric, underscore, hyphen, or dot."
                )

            # Check for duplicate sample names
            existing_ids = [s['sample'] for s in samples]
            if sample_id in existing_ids:
                errors.append(f"ERROR: Row {row_idx}: Duplicate sample name '{sample_id}'")

            # Validate file paths
            sr1 = row.get('short_reads_1', '').strip()
            sr2 = row.get('short_reads_2', '').strip()
            lr = row.get('long_reads', '').strip()
            assembly = row.get('assembly', '').strip()

            has_short = bool(sr1 and sr2)
            has_long = bool(lr)
            has_assembly = bool(assembly)

            if not has_short and not has_long and not has_assembly:
                errors.append(
                    f"ERROR: Row {row_idx}: Sample '{sample_id}' must have at least "
                    f"short reads, long reads, or a pre-assembled genome."
                )

            # Validate short reads come in pairs
            if bool(sr1) != bool(sr2):
                errors.append(
                    f"ERROR: Row {row_idx}: Sample '{sample_id}' has only one short-read file. "
                    f"Both R1 and R2 are required."
                )

            # Validate file extensions
            if sr1 and not any(sr1.endswith(ext) for ext in valid_extensions_fastq):
                errors.append(
                    f"WARNING: Row {row_idx}: Short reads file '{sr1}' has unexpected extension."
                )
            if lr and not any(lr.endswith(ext) for ext in valid_extensions_fastq):
                errors.append(
                    f"WARNING: Row {row_idx}: Long reads file '{lr}' has unexpected extension."
                )
            if assembly and not any(assembly.endswith(ext) for ext in valid_extensions_fasta):
                errors.append(
                    f"WARNING: Row {row_idx}: Assembly file '{assembly}' has unexpected extension."
                )

            # Determine assembly mode capability
            if has_short and has_long:
                mode = 'hybrid'
            elif has_long:
                mode = 'long'
            elif has_short:
                mode = 'short'
            elif has_assembly:
                mode = 'preassembled'
            else:
                mode = 'unknown'

            samples.append({
                'sample': sample_id,
                'short_reads_1': sr1 if sr1 else None,
                'short_reads_2': sr2 if sr2 else None,
                'long_reads': lr if lr else None,
                'assembly': assembly if assembly else None,
                'mode': mode
            })

    # Report results
    if errors:
        print("\n".join(errors), file=sys.stderr)
        error_count = sum(1 for e in errors if e.startswith("ERROR"))
        if error_count > 0:
            sys.exit(1)

    # Summary
    print(f"Samplesheet validation complete:")
    print(f"  Total samples: {len(samples)}")
    print(f"  Hybrid mode:   {sum(1 for s in samples if s['mode'] == 'hybrid')}")
    print(f"  Long-read:     {sum(1 for s in samples if s['mode'] == 'long')}")
    print(f"  Short-read:    {sum(1 for s in samples if s['mode'] == 'short')}")
    print(f"  Pre-assembled: {sum(1 for s in samples if s['mode'] == 'preassembled')}")

    return samples


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: check_samplesheet.py <samplesheet.csv>", file=sys.stderr)
        sys.exit(1)

    validate_samplesheet(sys.argv[1])
