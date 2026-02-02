import sys, os
import numpy as np
import pandas as pd
from sklearn import metrics as m
from argparse import ArgumentParser, Namespace
import logging
import traceback as tr

Log_Format = "%(levelname)s %(asctime)s - %(message)s"

logging.basicConfig(
    stream=sys.stdout,
    filemode="w",
    format=Log_Format,
    level=logging.INFO)

logger = logging.getLogger("Metrics")

VEPS = ['GenoCanyon', 'LRT', 'PhD_SNP', 'PANTHER', 'EVE', 'BayesDel', 'MutationTaster',
        'MVP', 'Polyphen2_HumDiv', 'fitCons', 'AlphaMissense', 'phastCons', 'PROVEAN', 'CAPICE',
        'EVmutation_epistatic', 'ClinPred', 'MOIpred_recesive', 'PrimateAI', 'COSMIS', 'M-CAP',
        'mutationTCN', 'MOIpred_dominant', 'MetaRNN', 'Rhapsody', 'ESCOTT', 'phyloP', 'PonP2', 'GERP++', 'MetaLR',
        'fathmm-XF',
        'Eigen-pred', 'SIFT4G', 'SPRI', 'VESPAl',
        'Envision', 'ESM-1b', 'PAPI', 'fathmm-MKL', 'VARITY_R', 'GEMME', 'SIFT', 'SuSPect', 'SNPred',
        'MutFormer', 'NetDiseaseSNP', 'PonPS', 'DeepSAV', 'VEST4', 'SNPs&GO', 'LASSIE', 'REVEL',
        'MOIpred_benign', 'LINSIGHT', 'S3D-PROF', 'InMeRF', 'popEVE', 'Eigen', 'SNPs&GO3D', 'ProGen2',
        'Grantham', 'MutScore', 'MutPred', 'Tranception', 'FATHMM', 'CONDEL', 'SIGMA', 'SNAP2', 'SiPhy',
        'MetaSVM', 'MISTIC', 'UNEECON', 'CADD', 'gMVP', 'MutationAssessor', 'iGEMME', 'DeMaSk', 'DANN',
        'EVmutation_independent', 'ESM-1v', 'DEOGEN2', 'CPT', 'MPC', 'VARITY_ER', 'Polyphen2_HumVar',
        'sequence_unet', 'LIST-S2', 'DeepSequence', 'AlphScore', 'BLOSUM62']

inverted_veps = []


def calculate_directionality(variants: pd.DataFrame, enriched_veps: list):
    global inverted_veps
    logger.info(f"Started Calculating Directionality")
    directions = {}
    for vep in enriched_veps:
        if vep == 'ESM-1v':
            directions['ESM-1v'] = 1
            inverted_veps.append(vep)
            continue
        curr_data = variants[['ESM-1v', vep]].dropna(subset=['ESM-1v', vep])
        corr_mat = curr_data.corr(method="spearman")
        val = corr_mat.loc[corr_mat.index == 'ESM-1v', vep].item()
        inverted = val > 0
        directions[vep] = inverted
        if inverted:
            inverted_veps.append(vep)
    return directions


def calculate_optimal_thresholds_per_region(data: pd.DataFrame,directions) -> pd.DataFrame:
    results = []
    cols = ['vep', 'threshold', 'region']
    existing_veps = [x for x in VEPS if x in data.columns.tolist()]
    for region in ['disordered', 'intermediate', 'ordered']:
        region_data = data[data['location'] == region]
        for vep in existing_veps:
            try:
                vep_data = region_data[[vep, "significance"]].dropna(subset=[vep])
                y, x = vep_data.pop('significance'), vep_data
                pos_label = 1 - directions[vep]
                fpr, tpr, thresholds = m.roc_curve(y, x, pos_label=pos_label)
                optimal_idx = np.argmax(tpr - fpr)
                results.append((vep, thresholds[optimal_idx], region))
            except Exception as e:
                logger.error(f"VEP: {vep} , Error: {str(e)}")
                tr.print_exc()
    return pd.DataFrame(results, columns=cols)


def calculate_optimal_thresholds(data: pd.DataFrame,directions) -> pd.DataFrame:
    results = []
    cols = ['vep', 'threshold']
    existing_veps = [x for x in VEPS if x in data.columns.tolist()]
    for vep in existing_veps:
        try:
            vep_data = data[[vep, "significance"]].dropna(subset=[vep])
            y, x = vep_data.pop('significance'), vep_data
            pos_label = 1 - directions[vep]
            fpr, tpr, thresholds = m.roc_curve(y, x, pos_label=pos_label)
            optimal_idx = np.argmax(tpr - fpr)
            results.append((vep, thresholds[optimal_idx]))
        except Exception as e:
            logger.error(f"VEP: {vep} , Error: {str(e)}")
            tr.print_exc()
    return pd.DataFrame(results, columns=cols)


def calcluate_metrics_per_region(data, regions_optimal,directions):
    metrics = []
    existing_veps = [x for x in VEPS if x in data.columns.tolist()]
    if 'EVE' in existing_veps:
        existing_veps.remove('EVE')
    logger.info(f"Before Culling: {data.shape[0]}")
    data.dropna(subset=existing_veps, how='any', inplace=True)
    logger.info(f"After Culling: {data.shape[0]}")
    for region in ['ordered', 'intermediate', 'disordered']:
        logger.info(f"Working on Structural Region: {region}...")
        for vep in existing_veps:
            try:
                pos_label = 1 - directions[vep]
                is_inverted = directions[vep]
                vep_data = data.loc[data['location'] == region, [vep, 'significance']]
                vep_data = vep_data.dropna(subset=[vep])
                y, x = vep_data.pop('significance'), vep_data
                optimal = regions_optimal.loc[(regions_optimal['vep'] == vep) & (regions_optimal['region'] == region),
                "threshold"].item()
                logger.info(f"VEP: {vep} , optimal: {optimal} , region: {region}")
                if is_inverted:
                    y_pred = x <= optimal
                else:
                    y_pred = x >= optimal
                tn, fp, fn, tp = m.confusion_matrix(y, y_pred).ravel()
                recall = m.recall_score(y, y_pred)
                precision = m.precision_score(y, y_pred)
                specificity = tn / (tn + fp)
                fpr = fp / (fp + tn)
                fnr = fn / (fn + tp)
                harmonic_bcc = 2 * (specificity * recall) / (specificity + recall)
                f1_score = (2 * recall * precision) / (recall + precision)
                precision_scores, recall_scores, _ = m.precision_recall_curve(y, x,
                                                                              pos_label=pos_label)

                roc_auc_score = m.roc_auc_score(y,y_pred)
                mcc = m.matthews_corrcoef(y, y_pred)
                pr_auc = m.auc(recall_scores, precision_scores)
                bacc = m.balanced_accuracy_score(y, y_pred)
                adj_bacc = m.balanced_accuracy_score(y, y_pred, adjusted=True)
                metrics.append(
                    (region, vep, tn, fp, fn, tp, recall, precision, specificity, fpr, fnr, harmonic_bcc, mcc,
                     pr_auc, bacc, adj_bacc, f1_score,roc_auc_score))

            except Exception as e:
                logger.error(f"VEP: {vep}, Error: {str(e)}")
                tr.print_exc()
    return pd.DataFrame(metrics,
                        columns=['region', 'vep', 'tn', 'fp', 'fn', 'tp', 'recall', 'precision', 'specificity', 'fpr',
                                 'fnr', 'harmonic_bacc', 'mcc',
                                 'pr_auc', 'bacc', 'adj_bacc', 'f1_score','auc'])


def process(args: Namespace):
    logger.info("Process Started")
    logger.info("Loading the missense variants dataset")
    data = pd.read_csv(args.missense, compression="gzip", low_memory=False)
    logger.info(f"Performing Threshing to remove predictors with 30% missing")
    # Dropping predictors with 30% missing
    thresh = 0.70 * data.shape[0]
    data.dropna(axis=1, thresh=thresh, inplace=True)
    existing_veps = [x for x in VEPS if x in data.columns.tolist()]
    logger.info("Missense variants dataset was loaded successfully")
    for vep in existing_veps:
        data[vep] = data[vep].replace("-", None)
        data[vep] = data[vep].replace(r'[^0-9.]', None, regex=True)
        data[vep] = data[vep].astype(float)
    data["significance"] = data["significance"].astype(int)
    logger.info(f"Calculating VEPs directionality")
    directions = calculate_directionality(data, existing_veps)
    logger.info(
        f"Calculating the optimal threshold of VEPs for different structural locations and across all structural regions")
    optimal_thresholds_without_regions = calculate_optimal_thresholds(data,directions)
    logger.info(f"Saving the optimal thresholds data")
    optimum_out = os.path.join(args.output, f"veps.optimal.thresholds.csv")
    optimal_thresholds_without_regions.to_csv(optimum_out, index=False)
    # calculating optimal thresholds per region as well
    logger.info("Keeping only rows where all VEPs do have scores for ")
    data.dropna(inplace=True)
    logger.info(f"Calculating different veps metrics per protein structural region")
    metrics_without_regions = calculate_metrics(data, optimal_thresholds_without_regions,directions)
    metrics_without_regions.to_csv(os.path.join(args.output, f"metrics.all.csv"), index=False)

    ## Calculate the same metrics but on region basis
    regions_optimal = calculate_optimal_thresholds_per_region(data,directions)
    regions_metrics = calcluate_metrics_per_region(data, regions_optimal,directions)
    regions_metrics.to_csv(
        os.path.join(args.output,f"regions.metrics.csv"),
        index=False
    )
    logger.info("All Done.")


def calculate_metrics(data, optimal_thresholds,directions) -> pd.DataFrame:
    metrics = []
    existing_veps = [x for x in VEPS if x in data.columns.tolist()]
    if 'EVE' in existing_veps:
        existing_veps.remove('EVE')
    logger.info(f"Before Culling: {data.shape[0]}")
    data.dropna(subset=existing_veps, how='any', inplace=True)
    logger.info(f"After Culling: {data.shape[0]}")
    for region in ['ordered', 'intermediate', 'disordered']:
        logger.info(f"Working on Structural Region: {region}...")
        for vep in existing_veps:
            try:
                pos_label = 1 - directions[vep]
                is_inverted = directions[vep]
                vep_data = data.loc[data['location'] == region, [vep, 'significance']]
                vep_data = vep_data.dropna(subset=[vep])
                y, x = vep_data.pop('significance'), vep_data
                optimal = optimal_thresholds.loc[(optimal_thresholds['vep'] == vep), "threshold"].item()
                logger.info(f"VEP: {vep} , optimal: {optimal}")
                if is_inverted:
                    y_pred = x <= optimal
                else:
                    y_pred = x >= optimal
                tn, fp, fn, tp = m.confusion_matrix(y, y_pred).ravel()
                recall = m.recall_score(y, y_pred)
                precision = m.precision_score(y, y_pred)
                specificity = tn / (tn + fp)
                fpr = fp / (fp + tn)
                fnr = fn / (fn + tp)
                harmonic_bcc = 2 * (specificity * recall) / (specificity + recall)
                f1_score = (2 * recall * precision) / (recall + precision)
                precision_scores, recall_scores, _ = m.precision_recall_curve(y, x,
                                                                              pos_label=pos_label)
                avg_precisions = np.mean(precision_scores)
                avg_recalls = np.mean(recall_scores)
                avg_precision = m.average_precision_score(y,x,pos_label=pos_label)
                mcc = m.matthews_corrcoef(y, y_pred)
                pr_auc = m.auc(recall_scores, precision_scores)
                bacc = m.balanced_accuracy_score(y, y_pred)
                adj_bacc = m.balanced_accuracy_score(y, y_pred, adjusted=True)
                roc_auc_score = m.roc_auc_score(y, x)
                metrics.append(
                    (region, vep, tn, fp, fn, tp, recall, precision, specificity, fpr, fnr, harmonic_bcc, mcc,
                     pr_auc, bacc, adj_bacc, f1_score, avg_precision,
                     avg_precisions, avg_recalls,roc_auc_score))

            except Exception as e:
                logger.error(f"VEP: {vep}, Error: {str(e)}")
                tr.print_exc()
    return pd.DataFrame(metrics,
                        columns=['region', 'vep', 'tn', 'fp', 'fn', 'tp', 'recall', 'precision', 'specificity', 'fpr',
                                 'fnr', 'harmonic_bacc', 'mcc',
                                 'pr_auc', 'bacc', 'adj_bacc', 'f1_score','avg_precision',
                                 'm_precision','m_recall','auc'])


def main():
    p = ArgumentParser(
        description="This program will calculate different confusion matrix metrics for all VEPS out there segregated by different "
                    "protein structural regions")
    p.add_argument("-m", "--missense", help="The missense variants dataset", required=True)
    p.add_argument("-o", "--output", help="The output directory to save the results", required=True)
    if len(sys.argv) <= 1:
        p.print_help()
        return
    args = p.parse_args()
    process(args)


if __name__ == '__main__':
    main()
