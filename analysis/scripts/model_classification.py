# script to perform model classification using XGBoost


## setup
from sys import argv
from glob import glob
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.pipeline import make_pipeline
from sklearn.metrics import accuracy_score, roc_auc_score
from skopt import BayesSearchCV
from skopt import dump as skdump, load as skload
from skopt.space import Real, Integer
# for CPU machines:
## mamba install -c conda-forge py-xgboost-cpu
# for GPU machines (when I run remote):
## mamba install -c conda-forge py-xgboost-gpu
from xgboost import XGBClassifier

# name arguments 
linear_folder = argv[1]
threshold_folder = argv[2]
## where to write tuning results
cv_out = argv[3]
## has to be "cpu" or "cuda"
device = argv[4]
## seed
s = int(argv[5])


## set random seed
np.random.seed(s)

## Import and clean data
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



################################
########## MAIN ################
################################


def main():
  ## read and clean data
  flist = [linear_folder, threshold_folder]

  tlist = ["linear", "threshold"]

  sumstats = combine_sumstats(folder_list = flist, transformation_list=tlist)
  
  ## test-train split

  # set predictor and response matrices
  X = sumstats.drop(columns = "transformation")
  y = sumstats["transformation"]

  # Encode the labels into unique integers
  encoder = LabelEncoder()
  y = encoder.fit_transform(np.ravel(y))

  # get a seed based on the global seed for input into function
  split_seed = np.random.randint(1, 1e6)

  X_train, X_test, y_train, y_test = train_test_split(
    X, y, 
    test_size=0.2, 
    random_state=split_seed
    )
  
  ## Create modeling pipeline
  pipe = make_pipeline(
    StandardScaler(),
    
    XGBClassifier(tree_method = "hist", device = device, eval_metric = "auc")
  )

  ## Tune

  ## set hyperparameters
  hparams = {
    "xgbclassifier__n_estimators": Integer(100, 2000),
    "xgbclassifier__eta": Real(1e-3, 1, "log-uniform"),
    "xgbclassifier__gamma": Real(1e-10, 1.5, "log-uniform"),
    "xgbclassifier__max_depth": Integer(1, 15),
    "xgbclassifier__subsample": Real(0.1, 1, "uniform"),
    "xgbclassifier__colsample_bytree": Real(0.1, 1.0, "uniform"),
    "xgbclassifier__alpha": Real(1e-3, 0.5, "uniform")
  }

  # get a seed based on the global seed for input into function
  cv_seed = np.random.randint(1, 1e6)


  searchcv = BayesSearchCV(
    pipe,
    search_spaces=hparams,
    n_iter = 50,
    cv = 5,
    scoring = "roc_auc",
    n_jobs = 5,
    random_state=cv_seed
  )

  np.int = int # have to do this because skopt hasn't totally updated their use of np.int vs int
  searchcv.fit(X_train, y_train)

  cv_out_path = cv_out

  skdump(searchcv, cv_out_path)



if __name__=="__main__":
  main()