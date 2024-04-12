# script to perform model classification using XGBoost


## setup
from sys import argv
from glob import glob
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import mean_absolute_error, mean_squared_error
from skopt import BayesSearchCV
from skopt import dump as skdump, load as skload
from skopt.space import Real, Integer
# for CPU machines:
## mamba install -c conda-forge py-xgboost-cpu
# for GPU machines (when I run remote):
## mamba install -c conda-forge py-xgboost-gpu
from xgboost import XGBRegressor

# name arguments 
folder = argv[1]
## transformation model to import
transformation = argv[2]
## species
species = argv[3]
## has to be "cpu" or "cuda"
device = argv[4]
## seed
s = int(argv[5])
## number of admixture populations
num_admix = int(argv[6])


## set random seed
np.random.seed(s)

## Import and clean data
# function to import data for model classification
def import_data(folder: str, transformation: str, num_admix: int) -> pd.DataFrame:

  # specify path to sumstats
  sumstats_paths = glob(f"{folder}/*{transformation}*sumstats.csv")

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


  # read in params to get max_k and Nm
  params_paths = glob(f"{folder}/*{transformation}*params.csv")

  if len(params_paths) > 1:
    params = pd.concat((pd.read_csv(f) for f in params_paths), ignore_index=True)
  else:
    params = pd.read_csv(params_paths)

  params = params.loc[:, ["max_k", "Nm"]]
  

  # read in demography summary to get admix pop sizes
  
  demo_paths = glob(f"{folder}/*{transformation}*demo_summary.csv")

  if len(demo_paths) > 1:
    demo = pd.concat((pd.read_csv(f) for f in demo_paths), ignore_index=True)
  else:
    demo = pd.read_csv(demo_paths)

  # filter for the columns that contain admixture population sizes
  admix_cols = []
  for i in range(1, num_admix + 1):
    admix_cols.append(f"admix_n_a{i}")

  demo = demo.loc[:, admix_cols]

  y = pd.concat([params, demo], axis = 1)

  X = sumstats

  return X, y

# get evaluation statistics
def get_eval_stats(y_emp, y_pred, eval_stat):
  
  stat_list = []
  for i in range(y_pred.shape[1]):
    emp = y_emp.iloc[:,i]

    pred = y_pred[:,i]

    if eval_stat == "mae":
      stat = mean_absolute_error(emp, pred)

      stat_list.append(stat)
    elif eval_stat == "rmse":
      stat = mean_squared_error(emp, pred)
      stat_list.append(stat)
    else:
      print("eval_stat should be mae or rmse")

  stat_df = pd.Series(stat_list, index=y_emp.columns)

  return(stat_df)



################################
########## MAIN ################
################################


def main():

  X, y = import_data(folder, transformation = transformation, num_admix = num_admix)

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
    XGBRegressor(tree_method = "hist", device = device, eval_metric = "mae")
  )

  ## Tune

  ## set hyperparameters
  hparams = {
    "xgbregressor__n_estimators": Integer(100, 2000),
    "xgbregressor__eta": Real(1e-2, 1, "log-uniform"),
    "xgbregressor__gamma": Real(1e-3, 1.5, "log-uniform"),
    "xgbregressor__max_depth": Integer(1, 15),
    "xgbregressor__subsample": Real(0.1, 1, "uniform"),
    "xgbregressor__colsample_bytree": Real(0.1, 1.0, "uniform"),
    "xgbregressor__alpha": Real(1e-3, 0.5, "uniform")
  }

  # get a seed based on the global seed for input into function
  cv_seed = np.random.randint(1, 1e6)

  searchcv = BayesSearchCV(
    pipe,
    search_spaces=hparams,
    n_iter = 100,
    cv = 5,
    scoring = "neg_mean_absolute_error",
    n_jobs = 5,
    random_state=cv_seed
  )

  np.int = int # have to do this because skopt hasn't totally updated their use of np.int vs int
  searchcv.fit(X_train, y_train)

  cv_out_path = f"reg_cv_tuning_{species}.pkl"

  # write tuning out
  skdump(searchcv, cv_out_path)

  # write train data
  X_train.to_csv(f"x_train_param_est_{species}.csv", index=False)
  y_train.to_csv(f"y_train_param_est_{species}.csv", index=False)

  # write test data
  X_test.to_csv(f"x_test_param_est_{species}.csv")
  y_test.to_csv(f"y_test_param_est_{species}.csv")


if __name__=="__main__":
  main()