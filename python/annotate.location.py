#!/usr/bin/env python
import os, sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser
import logging
from Bio.SeqIO import parse

Log_Format = "%(levelname)s %(asctime)s - %(message)s"

logging.basicConfig(
    stream=sys.stdout,
    filemode="w",
    format=Log_Format,
    level=logging.INFO
)

logger = logging.getLogger("IDRs_Features")

proteins = {}


def get_windows(prot_vector: np.ndarray, window_size: int, step_size: int) -> np.ndarray:
   
    if window_size <= 0:
        raise ValueError("window_size must be > 0")
    if step_size != 1:
        # Safer fallback if step_size is ever changed
        return np.array([prot_vector[i:i + window_size]
                         for i in range(0, len(prot_vector) - window_size + 1, step_size)], dtype=int)
    num_windows = (len(prot_vector) - window_size) // step_size + 1
    strides = prot_vector.strides + (prot_vector.strides[-1],)
    windows = np.lib.stride_tricks.as_strided(
        prot_vector, shape=(num_windows, window_size), strides=strides
    )
    return windows


def process(args):
    logger.info("Program started.")

    # Load protein sequences (not used for labeling here, but kept as in your original code)
    records = parse(args.fasta, "fasta")
    for record in records:
        proteins[record.id] = str(record.seq)
    logger.info(f"Done reading Protein Fasta sequences, Total Proteins: {len(proteins.keys())}")

    # Read list of protein IDs
    with open(args.proteins, "r") as reader:
        total_proteins = [x.strip() for x in reader.readlines() if x.strip()]

    if args.position < 1 or args.position > len(total_proteins):
        logger.info(f"--position {args.position} is out of range (1..{len(total_proteins)}). Aborting.")
        sys.exit(0)

    protein_id = total_proteins[args.position - 1]
    if protein_id not in proteins:
        logger.info(f"Protein ID: {protein_id} doesn't exist in the loaded multiple FASTA. Aborting.")
        sys.exit(0)

    logger.info("Loading AlphaFold residues out")
    af_residues = pd.read_csv(args.alphafold)

    af_residues = af_residues[["Uniprot", "residue", "amino acid", "SASA", "RSA", "score"]]
    af_residues.columns = ["uniprot_id", "pos", "aa", "sasa", "rsa", "plddt"]

    out_file = os.path.join(args.output, f"{protein_id}.csv")
    curr_prot_mask = af_residues["uniprot_id"] == protein_id
    current_protein = af_residues.loc[curr_prot_mask, ["aa", "pos", "plddt", "rsa"]].copy()
    current_protein["uniprot_id"] = protein_id

    # Basic sanity checks
    current_protein["pos"] = current_protein["pos"].astype(int)
    current_protein["plddt"] = current_protein["plddt"].astype(float)

    if current_protein.shape[0] < args.window:
        logger.info(f"Protein length < window size ({args.window}); ignoring {protein_id}.")
        return

    current_protein = current_protein.sort_values(by="pos", ascending=True).reset_index(drop=True)

    prot_vector = current_protein["pos"].to_numpy(dtype=int)
    windows = get_windows(prot_vector=prot_vector, window_size=args.window, step_size=args.step)

    positions = current_protein["pos"].to_numpy(dtype=int)
    votes = {int(p): np.zeros(3, dtype=np.int64) for p in positions}
    plddt_by_pos = dict(zip(current_protein["pos"].astype(int), current_protein["plddt"].astype(float)))

    for pos_window in windows:
        if len(pos_window) < args.window:
            continue
        w_plddt = np.array([plddt_by_pos[int(p)] for p in pos_window], dtype=float)
        window_is_ordered = (w_plddt.mean() >= 70.0)
        for p in pos_window:
            p = int(p)
            plddt = plddt_by_pos[p]

            if window_is_ordered:
                if plddt >= 70.0:
                    votes[p][2] += 1  # ordered
                else:
                    votes[p][1] += 1  # intermediate
            else:
                if plddt < 50.0:
                    votes[p][0] += 1  # disordered
                else:
                    votes[p][1] += 1  # intermediate

    label_map = {0: "disordered", 1: "intermediate", 2: "ordered"}
    final_labels = []
    for p in positions:
        p = int(p)
        v = votes[p]
        if v.sum() == 0:
            plddt = plddt_by_pos[p]
            if plddt < 70.0:
                final_labels.append("intermediate")
            else:
                final_labels.append("ordered")
        else:
            final_labels.append(label_map[int(np.argmax(v))])

    current_protein["location"] = final_labels
    current_protein.loc[
        (current_protein['location'].isna()) & (current_protein['plddt'] < 70), "location"] = "intermediate"
    current_protein.loc[(current_protein['location'].isna()) & (current_protein['plddt'] >= 70), "location"] = "ordered"
    logger.info(f"Saving protein's output file: {out_file}")
    current_protein.to_csv(out_file, index=False)
    logger.info("All Done.")


def main():
    p = ArgumentParser(
        description=(
            "Annotate each protein residue as disordered/intermediate/ordered using AlphaFold pLDDT "
            "with a 30-aa sliding window context. Outputs CSV with uniprot_id, position, residue, location."
        )
    )
    p.add_argument("-f", "--fasta", help="Protein FASTA file (multi-fasta).", required=True)
    p.add_argument("-a", "--alphafold", help="AlphaFold residues output CSV.", required=True)
    p.add_argument("-w", "--window", default=30, type=int,
                   help="Sliding window size. Default: 30 residues.")
    p.add_argument("-s", "--step", default=1, type=int,
                   help="Step size for sliding window. Default: 1 residue.")
    p.add_argument("-i", "--position", required=True, type=int,
                   help="1-based index of uniprot_id in the proteins list file.")
    p.add_argument("-P", "--proteins", required=True,
                   help="Text file containing all uniprot_ids, one per line.")
    p.add_argument("-o", "--output", required=True,
                   help="Output directory to save results to.")

    if len(sys.argv) <= 1:
        p.print_help()
        return

    args = p.parse_args()
    os.makedirs(args.output, exist_ok=True)
    process(args)


if __name__ == "__main__":
    main()

