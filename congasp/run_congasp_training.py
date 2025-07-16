import sys
sys.path.append('/app')
import argparse
import os
import pandas as pd
import torch
import numpy as np
import gc
from congas.Interface import Interface
from congas.models.LatentCategorical import LatentCategorical
from pyro.optim import ClippedAdam
from pyro.infer import TraceGraph_ELBO


def main():
    parser = argparse.ArgumentParser(description='Run CONGASp training with custom bin count data.')
    parser.add_argument('--separator', type=str, required=True, help='String to separate ATAC vs RNA cells (e.g., "146p")')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the output .npy file')
    parser.add_argument('--input_dir', type=str, default='.', help='Directory containing input CSV files (default: current directory)')
    args = parser.parse_args()

    try:
        sep = args.separator
        out_dir = args.output_dir
        in_dir = args.input_dir
        os.makedirs(out_dir, exist_ok=True)

        # Load bin count matrix (bins x cells)
        print(f"Loading data from {in_dir}...")
        bin_counts = pd.read_csv(os.path.join(in_dir, 'bin_counts_for_python.csv'), index_col=0)
        # Load ploidy file
        ploidy_df = pd.read_csv(os.path.join(in_dir, 'bin_ploidy.csv'))
        pld = torch.tensor(ploidy_df['ploidy'].values, dtype=torch.float32)

        # Split columns based on separator string
        atac_cols = [col for col in bin_counts.columns if sep in col]
        rna_cols = [col for col in bin_counts.columns if sep not in col]
        bin_counts_atac = bin_counts[atac_cols]
        bin_counts_rna = bin_counts[rna_cols]

        # Convert to torch tensors
        print("Converting to tensors...")
        data_rna = torch.tensor(bin_counts_rna.values, dtype=torch.float32)
        data_atac = torch.tensor(bin_counts_atac.values, dtype=torch.float32)
        norm_factor_rna = data_rna.sum(dim=0)
        norm_factor_atac = data_atac.sum(dim=0)
        segments = data_rna.shape[0]
        data_dict = {
            'data_rna': data_rna,
            'norm_factor_rna': norm_factor_rna,
            'data_atac': data_atac,
            'norm_factor_atac': norm_factor_atac,
            'pld': pld,
            'segments': segments
        }

        # Clear pandas dataframes to free memory
        del bin_counts, bin_counts_atac, bin_counts_rna, ploidy_df
        gc.collect()

        # Model setup
        print("Setting up model...")
        model = LatentCategorical
        optimizer = ClippedAdam
        loss = TraceGraph_ELBO
        interface = Interface(model, optimizer, loss)
        interface.initialize_model(data_dict)
        param_dict = {
            'K': 3,
            'theta_shape_rna': torch.ones(segments) * 20.0,
            'theta_rate_rna': torch.ones(segments) * 0.5,
            'theta_shape_atac': torch.ones(segments) * 10.0,
            'theta_rate_atac': torch.ones(segments) * 0.5,
            'lambda': 0.5,
            'multiome': False,
            'equal_mixture_weights': True,
            'likelihood_rna': 'NB',
            'likelihood_atac': 'NB',
            'nb_size_init_rna': torch.ones(segments) * 10.0,
            'nb_size_init_atac': torch.ones(segments) * 5.0,
            'hidden_dim': 5
        }
        interface.set_model_params(param_dict)
        
        print("Running training...")
        interface.run(steps=200)
        params = interface.learned_parameters()
        np.save(os.path.join(out_dir, 'congas_results.npy'), params)
        print(f"Saved results to {os.path.join(out_dir, 'congas_results.npy')}")
        
        # Clean up
        del interface, data_dict, params
        gc.collect()
        
    except Exception as e:
        print(f"Error processing {in_dir}: {str(e)}")
        raise

if __name__ == '__main__':
    main() 