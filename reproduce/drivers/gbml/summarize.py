"""Summarize source-feature observability and score denominators without cell IDs."""
import argparse
import json
from pathlib import Path
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); root = a.root
    if a.output.exists():
        raise FileExistsError(a.output)
    q = json.loads((root/'qualification/source_qualification.json').read_text())['result']
    feature = q['source_feature_ids']
    if not isinstance(feature, str) or q['source_features'] != 1:
        raise ValueError('This summary is defined for the original single-feature GBML181 input')
    pos = int(feature.split('_')[1])
    variants = pd.read_csv(root/'variants/MT_variants.tsv', sep='\t', dtype={'Cell': str, 'ID': str})
    variants = variants[variants.Cell_reads > 10]; present = set(variants.Cell)
    donors = [x for x in (root/'variants/donors.txt').read_text().splitlines() if x in present]
    targets = [x for x in (root/'variants/targets.txt').read_text().splitlines() if x in present]
    coverage = pd.read_csv(root/'variants/Coverage.tsv.gz', sep='\t', index_col=0).loc[pos]
    hits = variants[variants.ID == feature].set_index('Cell'); groups = {}
    for name, cells in (('source', donors), ('target', targets)):
        observed = [x for x in cells if coverage[x] >= 2]
        positive = [x for x in observed if x in hits.index and hits.loc[x, 's_reads'] > 0]
        groups[name] = {'observed_cells': len(observed), 'positive_cells': len(positive),
                        'positive_consensus_UMIs': int(hits.loc[positive, 's_reads'].sum())}
    stats = pd.read_csv(root/'qualification/source_feature_statistics.csv')
    stat = stats[stats.ID == feature].iloc[0]
    scores = pd.read_csv(root/'scores/score.tsv', sep='\t', dtype={'RNA_score': str})
    counts = scores.RNA_score.value_counts(); ties = counts[counts > 1]
    result = {'status': 'COMPLETE', 'feature': feature, 'observed_and_positive': groups,
              'odds_ratio': float(stat['odds ratio']), 'p_value': float(stat['p_value']), 'q_value': float(stat['q_values']),
              'RNA_tie_blocks': len(ties), 'RNA_cells_in_tie_blocks': int(ties.sum()), 'N': len(scores),
              'Neff_zero': int((scores.Neff == 0).sum()),
              'Neff_zero_observable': int(((scores.Neff == 0) & scores.observable).sum()),
              'Neff_zero_unobservable': int(((scores.Neff == 0) & ~scores.observable).sum())}
    a.output.parent.mkdir(parents=True, exist_ok=True); a.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
