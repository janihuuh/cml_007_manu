## Init
import os
import numpy as np
import numpy.random as random
import pandas as pd

from scvi.dataset.dataset import GeneExpressionDataset
from scvi.dataset.csv import CsvDataset
from scvi.inference import UnsupervisedTrainer
from scvi.models import SCANVI, VAE
from scvi.inference.autotune import auto_tune_scvi_model

from umap import UMAP

import torch
import scanpy as sc
import louvain

import logging
import pickle
from hyperopt import hp


# %matplotlib inline

use_cuda = True
n_epochs_all = None
save_path = ''
show_plot = True
os.chdir("/scratch/cs/csb/projects/NordCML/")


## CML (12 samples)
cml_706_dg   = CsvDataset(filename='results/scvi/input_files/706_dg.csv', save_path='', sep=',', new_n_genes=False)
cml_706_3mo  = CsvDataset(filename='results/scvi/input_files/706_3mo.csv', save_path='', sep=',', new_n_genes=False)
cml_706_12mo = CsvDataset(filename='results/scvi/input_files/706_12mo.csv', save_path='', sep=',', new_n_genes=False)

cml_716_dg   = CsvDataset(filename='results/scvi/input_files/716_dg.csv', save_path='', sep=',', new_n_genes=False)
cml_716_3mo  = CsvDataset(filename='results/scvi/input_files/716_3mo.csv', save_path='', sep=',', new_n_genes=False)
cml_716_12mo = CsvDataset(filename='results/scvi/input_files/716_12mo.csv', save_path='', sep=',', new_n_genes=False)

cml_720_dg   = CsvDataset(filename='results/scvi/input_files/720_dg.csv', save_path='', sep=',', new_n_genes=False)
cml_720_3mo  = CsvDataset(filename='results/scvi/input_files/720_3mo.csv', save_path='', sep=',', new_n_genes=False)
cml_720_12mo = CsvDataset(filename='results/scvi/input_files/720_12mo.csv', save_path='', sep=',', new_n_genes=False)

cml_730_dg   = CsvDataset(filename='results/scvi/input_files/730_dg.csv', save_path='', sep=',', new_n_genes=False)
cml_730_3mo  = CsvDataset(filename='results/scvi/input_files/730_3mo.csv', save_path='', sep=',', new_n_genes=False)
cml_730_12mo = CsvDataset(filename='results/scvi/input_files/730_12mo.csv', save_path='', sep=',', new_n_genes=False)



all_dataset = GeneExpressionDataset()
all_dataset.populate_from_per_batch_list(Xs = [cml_706_dg.X, cml_706_3mo.X, cml_706_12mo.X,
                                               cml_716_dg.X, cml_716_3mo.X, cml_716_12mo.X,
                                               cml_720_dg.X, cml_720_3mo.X, cml_720_12mo.X,
                                               cml_730_dg.X, cml_730_3mo.X, cml_730_12mo.X])

## Train, save and fin
vae      = VAE(all_dataset.nb_genes, n_batch=all_dataset.n_batches, n_labels=all_dataset.n_labels, n_hidden=128, n_latent=30, n_layers=2, dispersion='gene')
trainer  = UnsupervisedTrainer(vae, all_dataset, train_size=1.0)
trainer.train(n_epochs=100)
torch.save(trainer.model.state_dict(), 'results/scvi/results/cml_oneshot.pkl')
# trainer.model.load_state_dict(torch.load('results/scvi/output/fhrb1680_rm_oneshot.pkl'))

## Sample posterior to get latent representation and save those embeddings
full = trainer.create_posterior(trainer.model, all_dataset, indices=np.arange(len(all_dataset)))

latent, batch_indices, labels = full.sequential().get_latent()
batch_indices = batch_indices.ravel()

np.savetxt("results/scvi/results/cml_oneshot_latent.csv", latent, delimiter=",")
np.savetxt("results/scvi/results/ccml_oneshot_indices.csv", batch_indices, delimiter=",")
