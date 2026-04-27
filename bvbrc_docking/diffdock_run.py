#!/usr/bin/env python
"""
Simplified DiffDock wrapper for BV-BRC docking service.

This script directly invokes DiffDock inference without intermediate
configuration layers. It handles:
  1. Creating the input CSV from PDB + SMILES file
  2. Running DiffDock inference (supports both 1.0 and 1.1/L)
  3. Post-processing results (combined PDBs, GNINA scoring, result.csv)

Usage:
    diffdock_run.py --pdb protein.pdb --ligands ligands.smi --outdir ./results
    diffdock_run.py --pdb protein.pdb --ligands ligands.smi --outdir ./results --version 1.0

Environment:
    BVDOCK_DIFFDOCK_DIR: Path to DiffDock installation directory
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Use existing utility functions from bvbrc_docking.utils
from bvbrc_docking.utils import (
    clean_pdb,
    comb_pdb,
    sdf2pdb,
    validate_smiles,
    cal_cnn_aff,
    run_and_save,
)


def create_input_csv(pdb_file, ligands_file, output_dir, version="1.1"):
    """
    Create DiffDock input CSV from PDB and SMILES file.

    Args:
        pdb_file: Path to cleaned protein PDB
        ligands_file: Tab-separated file with ID and SMILES columns
        output_dir: Output directory
        version: DiffDock version ("1.0" or "1.1")

    Returns:
        Tuple of (csv_path, list of (id, smiles) tuples)
    """
    csv_path = os.path.join(output_dir, 'input.csv')
    ligands = []
    failed = 0

    with open(ligands_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            ident, smiles = parts[0], parts[1]

            # Validate SMILES, try swapping if needed
            if not validate_smiles(smiles):
                if validate_smiles(ident):
                    ident, smiles = smiles, ident
                else:
                    failed += 1
                    print(f"Warning: Invalid SMILES for {ident}", file=sys.stderr)
                    continue

            ligands.append((ident, smiles))

    if failed:
        print(f"Warning: {failed} ligands failed SMILES validation", file=sys.stderr)

    if not ligands:
        print("Error: No valid ligands found", file=sys.stderr)
        sys.exit(1)

    # Write CSV - format differs between versions
    with open(csv_path, 'w') as f:
        if version == "1.0":
            # DiffDock 1.0 format: protein_path,ligand
            f.write('protein_path,ligand\n')
            for ident, smiles in ligands:
                f.write(f'{pdb_file},{smiles}\n')
        else:
            # DiffDock 1.1/L format: protein_path,ligand_description,complex_name,protein_sequence
            f.write('protein_path,ligand_description,complex_name,protein_sequence\n')
            for ident, smiles in ligands:
                f.write(f'{pdb_file},{smiles},{ident},\n')

    print(f"Created input CSV with {len(ligands)} ligands (format: v{version})")
    return csv_path, ligands


def run_diffdock(diffdock_dir, csv_path, output_dir, version="1.1",
                 batch_size=10, samples_per_complex=10, inference_steps=20,
                 log_handle=None):
    """
    Run DiffDock inference.

    Args:
        diffdock_dir: Path to DiffDock installation
        csv_path: Input CSV file
        output_dir: Output directory
        version: DiffDock version ("1.0" or "1.1")
        batch_size: Inference batch size
        samples_per_complex: Number of poses per ligand
        inference_steps: Number of diffusion steps
        log_handle: File handle for logging
    """
    log_path = os.path.join(output_dir, 'diffdock_log')

    if version == "1.0":
        # DiffDock 1.0 requires separate ESM embedding step
        print("Running ESM embedding preparation (DiffDock 1.0)...")

        # Step 1: Prepare ESM input
        cmd_prep = (
            f"python {diffdock_dir}/datasets/esm_embedding_preparation.py "
            f"--protein_ligand_csv {csv_path} "
            f"--out_file {output_dir}/prepared_for_esm.fasta"
        )
        with open(log_path, 'w') as log:
            run_and_save(cmd_prep, cwd=output_dir, output_file=log)

        # Step 2: Run ESM embeddings
        print("Generating ESM embeddings...")
        model_def = os.getenv("BVDOCK_ESM_MODEL", "esm2_t33_650M_UR50D")
        cmd_esm = (
            f"python {diffdock_dir}/esm/scripts/extract.py "
            f"{model_def} {output_dir}/prepared_for_esm.fasta {output_dir}/esm2_output "
            f"--repr_layers 33 --include per_tok --truncation_seq_length 30000"
        )
        with open(log_path, 'a') as log:
            run_and_save(cmd_esm, cwd=output_dir, output_file=log)

        # Step 3: Run inference
        print("Running DiffDock 1.0 inference...")
        cmd_dock = (
            f"python {diffdock_dir}/inference.py "
            f"--protein_ligand_csv {csv_path} "
            f"--out_dir {output_dir} "
            f"--esm_embeddings_path {output_dir}/esm2_output "
            f"--cache_path {output_dir}/cache "
            f"--model_dir {diffdock_dir}/workdir/paper_score_model "
            f"--confidence_model_dir {diffdock_dir}/workdir/paper_confidence_model "
            f"--inference_steps {inference_steps} "
            f"--samples_per_complex {samples_per_complex} "
            f"--batch_size {batch_size}"
        )
        with open(log_path, 'a') as log:
            run_and_save(cmd_dock, cwd=diffdock_dir, output_file=log)

    else:
        # DiffDock 1.1/L - ESM embeddings are generated internally
        cmd = [
            'python', '-u', '-m', 'inference',
            '--protein_ligand_csv', csv_path,
            '--out_dir', output_dir,
        ]

        if batch_size > 0:
            cmd.extend(['--batch_size', str(batch_size)])

        if samples_per_complex > 0:
            cmd.extend(['--samples_per_complex', str(samples_per_complex)])

        if inference_steps > 0:
            cmd.extend(['--inference_steps', str(inference_steps)])

        print(f"Running DiffDock 1.1/L: {cmd}")

        with open(log_path, 'w') as log:
            result = subprocess.run(
                cmd,
                cwd=diffdock_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
            )

        if result.returncode != 0:
            print(f"Error: DiffDock failed. See {log_path}", file=sys.stderr)
            sys.exit(1)

    print("DiffDock completed successfully")


def post_process(output_dir, protein_pdb, ligands, gnina_path, top_n=3, version="1.1"):
    """
    Post-process DiffDock results.

    Args:
        output_dir: DiffDock output directory
        protein_pdb: Path to protein PDB
        ligands: List of (id, smiles) tuples
        gnina_path: Path to GNINA binary
        top_n: Number of top poses to keep
        version: DiffDock version ("1.0" or "1.1")
    """
    print(f"Post-processing {len(ligands)} ligands...")

    if version == "1.0":
        # DiffDock 1.0 output format: index{N}_{path}____{ligand}/
        _post_process_v1(output_dir, protein_pdb, ligands, gnina_path, top_n)
    else:
        # DiffDock 1.1/L output format: {ligand_id}/
        _post_process_v11(output_dir, protein_pdb, ligands, gnina_path, top_n)

    print("Post-processing complete")


def _post_process_v1(output_dir, protein_pdb, ligands, gnina_path, top_n):
    """
    Post-process DiffDock 1.0 results.
    Output directories are named: index{N}_{protein_path}____{smiles}
    """
    import pandas as pd

    csv_path = os.path.join(output_dir, 'input.csv')
    input_df = pd.read_csv(csv_path)
    output_rows = []

    for i, row in input_df.iterrows():
        prot_path = row['protein_path'].replace('/', '-')
        smiles = row['ligand']
        result_path = f"{output_dir}/index{i}_{prot_path}____{smiles}"

        if not os.path.isdir(result_path):
            print(f"Warning: No results for index {i}", file=sys.stderr)
            continue

        for j in range(top_n):
            sdf_files = glob.glob(f"{glob.escape(result_path)}/rank{j+1}_confidence-*.sdf")
            if not sdf_files:
                continue

            sdf = sdf_files[0]
            score = float(os.path.basename(sdf).split("-")[1][:-4])

            try:
                lig_pdb = sdf2pdb(sdf)
                combined_pdb = comb_pdb(protein_pdb, lig_pdb)

                # Run GNINA if available
                scores = {'CNNscore': 'NA', 'CNNaffinity': 'NA', 'Vinardo': 'NA'}
                if gnina_path and os.path.exists(gnina_path):
                    mol = cal_cnn_aff(protein_pdb, sdf, gnina_path)
                    if mol:
                        scores = {
                            'CNNscore': mol.data.get('CNNscore', 'NA'),
                            'CNNaffinity': mol.data.get('CNNaffinity', 'NA'),
                            'Vinardo': mol.data.get('minimizedAffinity', 'NA'),
                        }

                output_rows.append({
                    'ident': f'ligand_{i}',
                    'rank': j + 1,
                    'score': score,
                    'lig_sdf': os.path.basename(sdf),
                    'comb_pdb': os.path.basename(combined_pdb) if combined_pdb else 'NA',
                    'CNNscore': scores['CNNscore'],
                    'CNNaffinity': scores['CNNaffinity'],
                    'Vinardo': scores['Vinardo'],
                })
            except Exception as e:
                print(f"Warning: Failed to process rank {j+1} for index {i}: {e}", file=sys.stderr)

    # Write combined results
    result_csv = os.path.join(output_dir, 'result.csv')
    with open(result_csv, 'w') as f:
        f.write('\t'.join([
            'ident', 'rank', 'score', 'lig_sdf', 'comb_pdb',
            'CNNscore', 'CNNaffinity', 'Vinardo'
        ]) + '\n')
        for row in output_rows:
            f.write('\t'.join([str(row[k]) for k in [
                'ident', 'rank', 'score', 'lig_sdf', 'comb_pdb',
                'CNNscore', 'CNNaffinity', 'Vinardo'
            ]]) + '\n')


def _post_process_v11(output_dir, protein_pdb, ligands, gnina_path, top_n):
    """
    Post-process DiffDock 1.1/L results.
    Output directories are named by ligand identifier.
    """

    for ident, smiles in ligands:
        result_path = os.path.join(output_dir, ident)

        if not os.path.isdir(result_path):
            print(f"Warning: No results for {ident}", file=sys.stderr)
            continue

        # Find and parse rank files
        by_rank = []
        for filename in os.listdir(result_path):
            m = re.match(r'rank(\d+)_confidence(-?[\d.]+)\.sdf', filename)
            if not m:
                continue

            rank = int(m.group(1))
            confidence = float(m.group(2))

            # Filter by rank and confidence
            if top_n > 0 and rank > top_n:
                continue
            if confidence > 100:
                continue

            sdf_path = os.path.join(result_path, filename)

            # Convert to PDB and create complex using utils.py functions
            try:
                lig_pdb = sdf2pdb(sdf_path)
                combined_pdb = comb_pdb(protein_pdb, lig_pdb)
                if combined_pdb:
                    by_rank.append({
                        'ident': ident,
                        'rank': rank,
                        'confidence': confidence,
                        'sdf_path': sdf_path,
                        'combined_pdb': combined_pdb,
                    })
            except Exception as e:
                print(f"Warning: Failed to process {filename}: {e}", file=sys.stderr)

        # Sort by rank
        by_rank.sort(key=lambda x: x['rank'])

        # Run GNINA scoring and write results
        result_csv = os.path.join(result_path, 'result.csv')
        with open(result_csv, 'w') as f:
            f.write('\t'.join([
                'ident', 'rank', 'score', 'lig_sdf', 'comb_pdb',
                'CNNscore', 'CNNaffinity', 'Vinardo'
            ]) + '\n')

            for entry in by_rank:
                # Run GNINA if available
                scores = {'CNNscore': 'NA', 'CNNaffinity': 'NA', 'Vinardo': 'NA'}
                if gnina_path and os.path.exists(gnina_path):
                    mol = cal_cnn_aff(protein_pdb, entry['sdf_path'], gnina_path)
                    if mol:
                        scores = {
                            'CNNscore': mol.data.get('CNNscore', 'NA'),
                            'CNNaffinity': mol.data.get('CNNaffinity', 'NA'),
                            'Vinardo': mol.data.get('minimizedAffinity', 'NA'),
                        }

                f.write('\t'.join([
                    entry['ident'],
                    str(entry['rank']),
                    str(entry['confidence']),
                    os.path.basename(entry['sdf_path']),
                    os.path.basename(entry['combined_pdb']),
                    str(scores['CNNscore']),
                    str(scores['CNNaffinity']),
                    str(scores['Vinardo']),
                ]) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Run DiffDock molecular docking',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    parser.add_argument('--pdb', required=True,
                        help='Input protein PDB file')
    parser.add_argument('--ligands', required=True,
                        help='Input ligands file (tab-separated: ID SMILES)')
    parser.add_argument('--outdir', required=True,
                        help='Output directory')

    # Optional arguments
    parser.add_argument('--diffdock-dir',
                        default=os.getenv('BVDOCK_DIFFDOCK_DIR'),
                        help='DiffDock installation directory')
    parser.add_argument('--batch-size', type=int, default=10,
                        help='Inference batch size (reduce for large proteins)')
    parser.add_argument('--samples-per-complex', type=int, default=10,
                        help='Number of pose samples per ligand')
    parser.add_argument('--inference-steps', type=int, default=20,
                        help='Number of diffusion steps')
    parser.add_argument('--top-n', type=int, default=3,
                        help='Keep top N poses per ligand')
    parser.add_argument('--gnina', default=None,
                        help='Path to GNINA binary (optional, for CNN scoring)')
    parser.add_argument('--skip-postprocess', action='store_true',
                        help='Skip post-processing step')
    parser.add_argument('--version', choices=['1.0', '1.1'], default='1.1',
                        help='DiffDock version (1.0 requires separate ESM step, 1.1/L has internal ESM)')

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.pdb):
        print(f"Error: PDB file not found: {args.pdb}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.ligands):
        print(f"Error: Ligands file not found: {args.ligands}", file=sys.stderr)
        sys.exit(1)

    if not args.diffdock_dir:
        print("Error: DiffDock directory not specified. Set BVDOCK_DIFFDOCK_DIR or use --diffdock-dir", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(args.diffdock_dir):
        print(f"Error: DiffDock directory not found: {args.diffdock_dir}", file=sys.stderr)
        sys.exit(1)

    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)

    # Determine GNINA path - check multiple locations
    gnina_path = args.gnina
    if not gnina_path:
        # Check common GNINA locations in order of preference
        gnina_candidates = [
            os.environ.get('BVDOCK_GNINA_PATH'),  # Explicit environment variable
            shutil.which('gnina'),  # In PATH (e.g., from conda environment)
            os.path.join(args.diffdock_dir, 'gnina'),  # Legacy location
        ]
        for candidate in gnina_candidates:
            if candidate and os.path.exists(candidate):
                gnina_path = candidate
                break

    if gnina_path:
        print(f"GNINA: {gnina_path}")
    else:
        print("GNINA: not found (CNN scoring will be skipped)")

    print("=== DiffDock Docking ===")
    print(f"Version: {args.version}")
    print(f"Protein: {args.pdb}")
    print(f"Ligands: {args.ligands}")
    print(f"Output:  {args.outdir}")
    print(f"DiffDock: {args.diffdock_dir}")
    print()

    # Step 1: Clean PDB
    protein_name = Path(args.pdb).stem
    cleaned_pdb = os.path.join(args.outdir, f'{protein_name}_clean.pdb')
    print("Cleaning PDB...")
    cleaned_pdb = clean_pdb(args.pdb, cleaned_pdb)

    # Step 2: Create input CSV
    print("Preparing ligand inputs...")
    csv_path, ligands = create_input_csv(cleaned_pdb, args.ligands, args.outdir, version=args.version)

    # Step 3: Run DiffDock
    print("Running DiffDock inference...")
    run_diffdock(
        diffdock_dir=args.diffdock_dir,
        csv_path=csv_path,
        output_dir=args.outdir,
        version=args.version,
        batch_size=args.batch_size,
        samples_per_complex=args.samples_per_complex,
        inference_steps=args.inference_steps,
    )

    # Step 4: Post-process
    if not args.skip_postprocess:
        print("Post-processing results...")
        post_process(
            output_dir=args.outdir,
            protein_pdb=cleaned_pdb,
            ligands=ligands,
            gnina_path=gnina_path,
            top_n=args.top_n,
            version=args.version,
        )

    print()
    print("=== Done ===")
    print(f"Results in: {args.outdir}")


if __name__ == '__main__':
    main()
