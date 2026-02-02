import os, sys
import numpy as np
import pandas as pd
import logging
from argparse import ArgumentParser, Namespace
import traceback as tr
from collections import defaultdict
from sklearn import metrics as m
from itertools import combinations

random_state = 420
# Logging
Log_Format = "%(levelname)s %(asctime)s - %(message)s"
logging.basicConfig(
    stream=sys.stdout,
    filemode="w",
    format=Log_Format,
    level=logging.INFO)

logger = logging.getLogger("Compare AUCs")


VEPS = ['GenoCanyon', 'LRT', 'PhD_SNP', 'PANTHER', 'EVE', 'BayesDel',
        'MutationTaster', 'MVP', 'Polyphen2_HumDiv', 'fitCons', 'AlphaMissense',
        'phastCons', 'PROVEAN', 'CAPICE', 'EVmutation_epistatic',
        'ClinPred', 'MOIpred_recesive', 'PrimateAI', 'COSMIS', 'M-CAP',
        'mutationTCN', 'MOIpred_dominant', 'MetaRNN', 'Rhapsody', 'ESCOTT',
        'phyloP', 'PonP2', 'GERP++', 'MetaLR', 'fathmm-XF', 'Eigen-pred', 'SIFT4G',
        'SPRI', 'VESPAl', 'Envision', 'ESM-1b', 'PAPI', 'fathmm-MKL', 'VARITY_R', 'GEMME',
        'SIFT', 'SuSPect', 'SNPred', 'MutFormer', 'NetDiseaseSNP', 'PonPS', 'DeepSAV', 'VEST4',
        'SNPs&GO', 'LASSIE', 'REVEL', 'MOIpred_benign', 'LINSIGHT', 'S3D-PROF', 'InMeRF', 'popEVE',
        'Eigen', 'SNPs&GO3D', 'ProGen2', 'Grantham', 'MutScore', 'MutPred', 'Tranception', 'FATHMM', 'CONDEL',
        'SIGMA', 'SNAP2', 'SiPhy', 'MetaSVM', 'MISTIC', 'UNEECON', 'CADD', 'gMVP', 'MutationAssessor', 'iGEMME',
        'DeMaSk', 'DANN', 'EVmutation_independent', 'ESM-1v', 'DEOGEN2', 'CPT', 'MPC', 'VARITY_ER',
        'Polyphen2_HumVar', 'sequence_unet', 'LIST-S2', 'DeepSequence', 'AlphScore', 'BLOSUM62']



enriched_veps = []

def change_to_numeric(num):
    try:
        return float(num)
    except Exception as e:
        logger.error(str(e))
        return 0.0


def load_variants(mutations):
    variants = None
    chunks = pd.read_csv(mutations, compression="gzip", low_memory=False,
                         chunksize=100_000)
    for chunk in chunks:
        if variants is None:
            variants = chunk
        else:
            variants = pd.concat([variants, chunk])
    return variants


def calculate_directionality(variants):
    logger.info(f"Started Calculating Directionality")
    directions = {}
    for vep in VEPS:
        if vep == 'ESM-1v':
            continue
        curr_data = variants[['ESM-1v', vep]].dropna(subset=['ESM-1v', vep])
        corr_mat = curr_data.corr(method="spearman")
        val = corr_mat.loc[corr_mat.index == 'ESM-1v', vep].item()
        inverted = val > 0
        directions[vep] = inverted
    return directions


def calculate_thresholds(variants, directions):
    optimal_results = {}
    for vep in enriched_veps:
        try:
            vep_data = variants[['significance', vep]].dropna(subset=[vep])
            pos_lbl = 1 - directions[vep]
            fpr, tpr, thresholds = m.roc_curve(vep_data['significance'], vep_data[vep],
                                               pos_label=pos_lbl)
            optimal_idx = np.argmax(tpr - fpr)

            optimal_results[vep] = thresholds[optimal_idx]
        except Exception as e:
            logger.error(str(e))
            continue

    return optimal_results


def calculate_agreement(variants, thresholds, directions):
    agreement = []
    for vep in enriched_veps:
        vep_data = variants[['significance', vep]].dropna(subset=[vep])
        if not vep in thresholds.keys():
            continue
        vep_threshold = thresholds[vep]
        is_inverted = directions[vep]
        if is_inverted:
            vep_data['agreement'] = np.where(vep_data[vep] < vep_threshold, 1, 0)
        else:
            vep_data['agreement'] = np.where(vep_data[vep] > vep_threshold, 1, 0)
        try:
            #vep_rating_data , _ = rating.aggregate_raters(vep_data[['significance', 'agreement']], n_cat=2)
            kappa = m.cohen_kappa_score(vep_data['significance'],vep_data['agreement'])
            agreement.append((vep, kappa))
        except Exception as e:
            logger.error(f"{vep} : {str(e)}")
            continue
    return agreement


def calculate_agreement_regions(variants, thresholds, directions):
    agreement = []
    for region in ['disordered', 'intermediate', 'ordered']:
        for vep in enriched_veps:
            vep_data = variants.loc[variants['location'] == region,
            ['significance', vep]].dropna(subset=[vep])
            if not vep in thresholds.keys():
                continue
            vep_threshold = thresholds[vep]
            is_inverted = directions[vep]
            if is_inverted:
                vep_data['agreement'] = np.where(vep_data[vep] < vep_threshold, 1, 0)
            else:
                vep_data['agreement'] = np.where(vep_data[vep] > vep_threshold, 1, 0)
            try:
                #vep_rating_data, _ = rating.aggregate_raters(vep_data[['significance', 'agreement']], n_cat=2)
                kappa = m.cohen_kappa_score(vep_data['significance'],vep_data['agreement'])
                agreement.append((vep, region, kappa))
            except Exception as e:
                logger.error(f"{vep} : {str(e)}")
                continue
    return agreement


def calculate_pairwise_agreement(variants, thresholds, directions):
    agreement = []
    groups = combinations(enriched_veps, 2)
    groups = [x for x in groups]
    logger.info(f"Pairwise Groups count: {len(groups)}")
    for counter, group in enumerate(groups):
        logger.info(f"Working on: {counter+1} / {len(groups)}")
        first_vep, second_vep = group
        group_data = variants[[first_vep, second_vep]].copy().dropna(subset=group)
        # The optimal threshold
        if first_vep not in thresholds.keys() or second_vep not in thresholds.keys():
            continue
        first_threshold = thresholds[first_vep]
        second_threshold = thresholds[second_vep]
        # the direction
        first_is_inverted = directions[first_vep]
        second_is_inverted = directions[second_vep]
        if first_is_inverted:
            group_data['first_agreement'] = np.where(group_data[first_vep] < first_threshold, 1, 0)
        else:
            group_data['first_agreement'] = np.where(group_data[first_vep] > first_threshold, 1, 0)

        if second_is_inverted:
            group_data['second_agreement'] = np.where(group_data[second_vep] < second_threshold, 1, 0)
        else:
            group_data['second_agreement'] = np.where(group_data[second_vep] > second_threshold, 1, 0)
        try:
            #group_data_rating , _ = rating.aggregate_raters(group_data[['first_agreement', 'second_agreement']], n_cat=2)
            kappa = m.cohen_kappa_score(group_data['first_agreement'],group_data['second_agreement'])
            agreement.append((first_vep, second_vep, kappa))
        except Exception as e:
            logger.error(f"({first_vep} , {second_vep}) : {str(e)}")
            continue
    return agreement


def regions_pairwise_agreement(variants, thresholds, directions):
    agreement = []
    groups = combinations(enriched_veps, 2)
    groups = [x for x in groups]
    logger.info(f"Pairwise Groups count: {len(groups)}")
    for counter,group in enumerate(groups):
        logger.info(f"Working on: {counter+1} / {len(groups)}")
        first_vep, second_vep = group
        for region in ['disordered', 'intermediate', 'ordered']:
            group_data = variants.loc[variants['location'] == region, [first_vep, second_vep]].copy().dropna(subset=group)
            # The optimal threshold
            if first_vep not in thresholds.keys() or second_vep not in thresholds.keys():
                continue
            first_threshold = thresholds[first_vep]
            second_threshold = thresholds[second_vep]
            # the direction
            first_is_inverted = directions[first_vep]
            second_is_inverted = directions[second_vep]
            if first_is_inverted:
                group_data['first_agreement'] = np.where(group_data[first_vep] < first_threshold, 1, 0)
            else:
                group_data['first_agreement'] = np.where(group_data[first_vep] > first_threshold, 1, 0)

            if second_is_inverted:
                group_data['second_agreement'] = np.where(group_data[second_vep] < second_threshold, 1, 0)
            else:
                group_data['second_agreement'] = np.where(group_data[second_vep] > second_threshold, 1, 0)
            try:
                #group_data_rating , _ = rating.aggregate_raters(group_data[['first_agreement', 'second_agreement']],                                           n_cat=2)
                kappa = m.cohen_kappa_score(group_data['first_agreement'],group_data['second_agreement'])
                agreement.append((first_vep, second_vep, region, kappa))
            except Exception as e:
                logger.error(f"({first_vep} , {second_vep}) : {str(e)}")
                continue
    return agreement


def process(args):
    global enriched_veps
    logger.info(f"Loading missense dataset: {args.mutations}")
    variants = load_variants(args.mutations)
    logger.info(f"Variants were loaded: {variants.shape[0]}")
    logger.info("Coercing non-numeric to numeric values..")
    for vep in VEPS:
        variants[vep] = variants[vep].replace("-", None)
        variants[vep] = variants[vep].replace(r'[^0-9.]', None, regex=True)
        variants[vep] = variants[vep].astype(float)
    logger.info(f"Checking enriched VEPs")
    enriched_data = variants[VEPS].dropna(thresh=0.75*variants.shape[0],axis=1)
    enriched_veps = [x for x in enriched_data.columns.tolist() if x in VEPS]
    logger.info(f"Calculating Predictors Directionality")
    directions = calculate_directionality(variants)
    logger.info(f"VEPs directionality was calculated. Saving it..")
    out_file = os.path.join(args.output, f"veps.directions.csv")
    directions_df = pd.DataFrame(directions.items(), columns=['vep', 'direction'])
    directions_df.to_csv(out_file, index=False)
    logger.info(f"Calculating optimal thresholds")
    thresholds = calculate_thresholds(variants, directions)
    logger.info(f"Optimal Thresholds were calculated. Saving the results....")
    out_file = os.path.join(args.output, "veps.optimal.thresholds.csv")
    thresholds_df = pd.DataFrame(thresholds.items(), columns=['vep', 'threshold'])
    thresholds_df.to_csv(out_file, index=False)
    logger.info("Now Calculating Inter-rating agreement among predictors across structural regions")
    logger.info("First calculating their agreement with the ground truth")
    ground_agreement = calculate_agreement(variants, thresholds, directions)
    ground_agreement_df = pd.DataFrame(ground_agreement, columns=['vep', 'kappa'])
    out_file = os.path.join(args.output, f"veps.agreement.significance.csv")
    ground_agreement_df.to_csv(out_file, index=False)
    logger.info(f"Calculating the agreement with the ground truth across structural regions")
    regions_agreement = calculate_agreement_regions(variants, thresholds, directions)
    regions_agreement_df = pd.DataFrame(regions_agreement, columns=['vep', 'region', 'kappa'])
    out_file = os.path.join(args.output, "veps.agreement.regions.significance.csv")
    regions_agreement_df.to_csv(out_file, index=False)
    logger.info(f"Calculating pairwise agreement of predictors")
    ground_agreement = calculate_pairwise_agreement(variants, thresholds, directions)
    logger.info(f"Pairwise agreement among predictors has finished. saving results....")
    ground_agreement_df = pd.DataFrame(ground_agreement, columns=['first_vep', 'second_vep', 'kappa'])
    out_file = os.path.join(args.output, f"veps.pairwise.agreement.csv")
    ground_agreement_df.to_csv(out_file, index=False)
    logger.info(f"Calculating pairwise agreement across structural regions")
    rpa = regions_pairwise_agreement(variants, thresholds, directions)
    logger.info(f"Regions pairwise agreement calculation has finished. saving the results")
    regions_pairwise_agreement_df = pd.DataFrame(rpa, columns=['first_vep', 'second_vep',
                                                                                      'region', 'kappa'])
    out_file = os.path.join(args.output, f"veps.pairwise.agreement.regions.csv")
    regions_pairwise_agreement_df.to_csv(out_file, index=False)
    logger.info("All Done.")


def main():
    p = ArgumentParser(description="This script will calculate the inter-predictors agreement using Fleiss Kappa")
    p.add_argument("-m", "--mutations", help="The missense mutations dataset",required=True)
    p.add_argument("-o", "--output", help="The output directory where to save the results",required=True)
    if len(sys.argv) <= 1:
        p.print_help()
        return
    args = p.parse_args()
    process(args)


if __name__ == '__main__':
    main()
