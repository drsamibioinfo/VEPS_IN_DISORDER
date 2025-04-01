#!/usr/bin/env python
import os, sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser, Namespace
import logging
from Bio.SeqIO import parse

Log_Format = "%(levelname)s %(asctime)s - %(message)s"

logging.basicConfig(
    stream=sys.stdout,
    filemode="w",
    format=Log_Format,
    level=logging.INFO)

logger = logging.getLogger("IDRs_Features")

proteins = {}


def get_windows(prot_vector: np.ndarray, window_size: int, step_size: int) -> np.ndarray:
    # Calculate the number of windows
    num_windows = (len(prot_vector) - window_size) // step_size + 1
    # Create a sliding window view of the vector
    strides = prot_vector.strides + (prot_vector.strides[-1],)
    windows = np.lib.stride_tricks.as_strided(prot_vector, shape=(num_windows, window_size), strides=strides)
    return windows


def get_blocks(windows, pos):
    blocks = []
    for idx, window in enumerate(windows):
        if pos in windows:
            blocks.append(idx)
    return blocks


def get_func_score(func, sequence, windows):
    from utils.sequence import Sequence
    score = 0.0
    try:
        for window in windows:
            current_seq = "".join([sequence[x] for x in window.tolist()])
            try:
                seq = Sequence(current_seq)
                val = getattr(seq, func)()
                if (func == "kappa") and val == -1:
                    continue
                else:
                    score += val
            except Exception:
                continue
    except Exception as e:
        return None
    finally:
        return score / len(windows)


def process(args):
    logger.info(f"Program started.")
    records = parse(args.fasta, "fasta")
    for record in records:
        proteins[record.id] = str(record.seq)
    logger.info(f"Done reading Protein Fasta sequences , Total Proteins : {len(proteins.keys())}")
    total_proteins = []
    with open(args.proteins,"r") as reader:
        total_proteins = reader.readlines()
        total_proteins = [x.replace("\n","") for x in total_proteins if len(x) > 0]

    protein_id = total_proteins[args.position - 1]
    if protein_id not in proteins.keys():
        logger.info(f"Protein ID: {protein_id} doesn't exist in the loaded multiple Fasta. Aborting....")
        sys.exit(0)
    logger.info(f"Loading Alphafold residues out")
    af_residues = pd.read_csv(args.alphafold)
    af_residues = af_residues[["Uniprot", "residue", "amino acid", "SASA", "RSA", "score"]]
    af_residues.columns = ["uniprot_id", "pos", "aa", "sasa", "rsa", "plddt"]
    out_file = os.path.join(args.output, f"{protein_id}.csv")
    curr_prot_mask = af_residues['uniprot_id'] == protein_id
    current_protein: pd.DataFrame = af_residues.loc[curr_prot_mask, ['aa', 'pos', 'plddt', 'rsa']].copy()
    current_protein["uniprot_id"] = protein_id
    if current_protein.shape[0] < 30:
        logger.info(f"protein is lower than 30 residues, ignoring it...")
        return
    prot_vector = current_protein.sort_values(by='pos', ascending=True)['pos'].to_numpy()
    windows = get_windows(prot_vector=prot_vector, window_size=args.window, step_size=args.step)
    for pos_window in windows:
        if len(pos_window) < 30:
            continue
        plddt_window = current_protein.loc[current_protein['pos'].isin(pos_window), "plddt"].to_numpy()
        pos_mask = current_protein['pos'].isin(pos_window)
        intermediate_mask = current_protein["plddt"] >= 50
        idr_mask = current_protein["plddt"] < 50
        ordered_mask = current_protein["plddt"] >= 70

        if np.sum(plddt_window) < (70 * len(pos_window)):
            current_protein.loc[pos_mask & idr_mask, "location"] = "disordered"
            current_protein.loc[pos_mask & intermediate_mask, "location"] = "intermediate"
        else:
            current_protein.loc[pos_mask & idr_mask, "location"] = "intermediate"
            current_protein.loc[pos_mask & ordered_mask, "location"] = "ordered"

    logger.info(f"Saving protein's output file: {out_file}")
    current_protein.loc[(current_protein['location'].isna()) & (current_protein['plddt'] < 50), "location"] = "disordered"
    current_protein.loc[(current_protein['location'].isna()) & (current_protein['plddt'] < 70), "location"] = "intermediate"
    current_protein.loc[(current_protein['location'].isna()) & (current_protein['plddt'] > 70), "location"] = "ordered"
    current_protein.to_csv(out_file, index=False)
    logger.info("All Done.")


def main():
    p = ArgumentParser(
        description="This program will annotate each protein residue into disordered,intermediate or ordered structural region."
                    "The program takes a multiple fasta file and alphafold models file and create a CSV file with uniprot_id, position, residue, region")
    p.add_argument("-f", "--fasta", help="The protein fasta file to use", required=True)
    p.add_argument("-a", "--alphafold", help="The alphafold residues out.", required=True)
    p.add_argument("-w", "--window", default=30, type=int, help="The window size used to slide over the entire protein "
                                                                "sequence. Defaults : 30 residues")
    p.add_argument("-s", "--step", default=1, type=int, help="The step size used to slide over the entire protein "
                                                             "sequence with the window size. Defaults: 1 residue")
    p.add_argument("-i", "--position", help="The position of uniprot_id in the file to use", required=True,type=int)
    p.add_argument("-P","--proteins",help="The text file containing all uniprot_ids one per line",required=True)
    p.add_argument("-o", "--output", help="The output directory to save results to")
    if len(sys.argv) <= 1:
        p.print_help()
        return
    args = p.parse_args()
    process(args)


if __name__ == '__main__':
    main()
