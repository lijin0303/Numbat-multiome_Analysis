#!/usr/bin/env python3

import re
try:
    import numpy as np
    import pandas as pd
except ImportError as e:
    print("Missing required package. Please install numpy and pandas: pip install numpy pandas")
    raise e

def parse_coordinates(bin_name):
    """Parse bin name like 'chr7:100000001-100200000' into seqnames, start, end"""
    match = re.match(r'(.+):(\d+)-(\d+)', bin_name)
    if match:
        seqnames = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        return seqnames, start, end
    else:
        raise ValueError(f"Cannot parse bin name: {bin_name}")

def main():
    # Load CONGAS results
    results_path = './congas_results.npy'
    params = np.load(results_path, allow_pickle=True).item()
    CNA = params['CNA']  # (K, segments)
    K, segments = CNA.shape
    
    # Use cluster assignments (RNA preferred)
    cluster_assignments = params.get('cluster_assignments_rna', None)
    if cluster_assignments is None:
        cluster_assignments = params.get('cluster_assignments_atac', None)
    if cluster_assignments is None:
        raise ValueError('No cluster assignments found in CONGAS results!')
    
    # Load bin segments
    bin_segments = pd.read_csv('./bin_segments.csv')
    bin_names = bin_segments['bin'].tolist()
    
    # Create CNV segments output (one per segment, assigned cluster only)
    cnv_segments = []
    for seg_idx, bin_name in enumerate(bin_names):
        assigned_cluster = cluster_assignments[seg_idx]
        cn = CNA[assigned_cluster, seg_idx]
        if cn == 2:
            continue  # Only output CNV segments
        seqnames, start, end = parse_coordinates(bin_name)
        cnv_type = 'amp' if cn > 2 else 'del'
        cnv_segments.append({
            'seqnames': seqnames,
            'start': start,
            'end': end,
            'cnv': cnv_type,
            'cluster': assigned_cluster,
            'copy_number': cn
        })
    
    # Create DataFrame and save
    cnv_df = pd.DataFrame(cnv_segments)
    if len(cnv_df) > 0:
        cnv_df = cnv_df.sort_values(['seqnames', 'start', 'cluster'])
        cnv_df[['seqnames', 'start', 'end', 'cnv', 'cluster']].to_csv(
            './cnv_segments_formatted.tsv', 
            sep='\t', 
            index=False
        )
        cnv_df.to_csv('./cnv_segments_detailed.tsv', sep='\t', index=False)
        print(f"Generated CNV segments file with {len(cnv_df)} entries (collapsed by assigned cluster)")
        print(cnv_df[['seqnames', 'start', 'end', 'cnv', 'cluster']].head(10))
    else:
        print("No copy number variations found!")

if __name__ == "__main__":
    main() 