"""
Script to write out model classification train/test splits and predictions. I didn't remember to write these out in the original modeling workflow, so I have to do it here. 
"""
import sys
import argparse
import os
from glob import glob
import numpy as np
import pandas as pd
from skopt import load as skload
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder


# function to import data for model classification
def import_data(folder: str, transformation: str) -> pd.DataFrame:

  # specify path to sumstats
  sumstats_paths = glob(folder + "/*sumstats.csv")

  # only concatenate if there is more than one file
  if len(sumstats_paths) > 1:
    sumstats = pd.concat((pd.read_csv(f) for f in sumstats_paths), ignore_index=True)
  else:
    sumstats = pd.read_csv(sumstats_paths)
  

  # remove the columns we're not using for analysis (per-deme pi, the third hill number, morans I, IDs)
  sumstats = sumstats.loc[:,~sumstats.columns.str.startswith('pi_pop')]
  sumstats = sumstats.loc[:,~sumstats.columns.str.startswith('sfs_h3')]
  sumstats = sumstats.drop(columns = "morans_i")
  sumstats = sumstats.loc[:,~sumstats.columns.str.contains("_id")]

  sumstats["transformation"] = transformation


  return sumstats

# path order must correspond with transformation order
def combine_sumstats(folder_list: list, transformation_list = list) -> pd.DataFrame:
  
  ss_list = []
  for p in range(len(folder_list)):
    ss = import_data(folder_list[p], transformation=transformation_list[p])
    ss_list.append(ss)

  sumstats = pd.concat(ss_list, ignore_index = True)

  return sumstats

def main():

  parser = argparse.ArgumentParser(description='Script to write out model classification train/test split and predictions')
  parser.add_argument('--species', type=str, help='Species three letter code')
  parser.add_argument('--cv_path', type=str, help='Path to cross validation results file')
  parser.add_argument('--sims_dir', type=str, help='Directory containing simulation results')
  parser.add_argument('--out_dir', type=str, help='Output directory path')
  parser.add_argument('--seed', type=int, help='Seed for random number generator', default = 1234)

  args = parser.parse_args()

  species = args.species
  cv_path = args.cv_path
  sims_dir = args.sims_dir
  out_dir = args.out_dir
  seed = int(args.seed)

  # load data from file paths
  searchcv = skload(cv_path)

  np.random.seed(seed)

  # read in the sumstats
  flist = glob(sims_dir + f"/*{species}*")
  tlist = ["linear", "threshold"]
  sumstats = combine_sumstats(folder_list = flist, transformation_list=tlist)

  # set predictor and response matrices
  X = sumstats.drop(columns = "transformation")
  y1 = sumstats["transformation"]

  # Encode the labels into unique integers
  encoder = LabelEncoder()
  y = encoder.fit_transform(np.ravel(y1))

  # get a seed based on the global seed for input into function
  split_seed = np.random.randint(1, 1e6)

  X_train, X_test, y_train, y_test = train_test_split(
    X, y, 
    test_size=0.2, 
    random_state=split_seed)
  
  # write out train/test split
  X_train.to_csv(os.path.join(out_dir, f"X_train_{species}.csv"), index = False)
  X_test.to_csv(os.path.join(out_dir, f"X_test_{species}.csv"), index = False)
  np.savetxt(os.path.join(out_dir, f"y_train_{species}.csv"), y_train, delimiter = ",")
  np.savetxt(os.path.join(out_dir, f"y_test_{species}.csv"), y_test, delimiter = ",")

  # fit the model
  best_pipeline = searchcv.best_estimator_
  best_pipeline.fit(X_train, y_train)
  test_pred_bin = best_pipeline.predict(X_test)
  test_pred_prob = best_pipeline.predict_proba(X_test)

  # write out predictions
  np.savetxt(os.path.join(out_dir, f"test_pred_bin_{species}.csv"), test_pred_bin, delimiter = ",")
  np.savetxt(os.path.join(out_dir, f"test_pred_prob_{species}.csv"), test_pred_prob, delimiter = ",")


if __name__ == "__main__":
  main()