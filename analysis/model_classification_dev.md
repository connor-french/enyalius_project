# Model classification in python dev


## Setup

``` python
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
```

## Set random seed

``` python
seed np.random.seed(1111)
```

## Import and clean data

``` python
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


flist = ["output/simulations/sims_ihe_linear", "output/simulations/sims_ihe_hinge"]

tlist = ["linear", "threshold"]

sumstats = combine_sumstats(folder_list = flist, transformation_list=tlist)
```

## Split into train/test split

``` python
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
  random_state=split_seed)
```

## Baseline model

Creating a baseline model to compare the tuned model to.

> From this first go with some somewhat arbitrary values, the accuracy
> and the AUC are 0.82 on the training data (they differ very slightly).
> The test accuracy and AUC are 0.76 (they also differ slightly).

``` python
pipe_baseline = make_pipeline(
  StandardScaler(),
  XGBClassifier(n_estimators=100, max_depth=4, learning_rate=0.5, objective='binary:logistic')
)

pipe_baseline.fit(X_train, y_train)

# predict to training
y_pred = pipe_baseline.predict(X_train)

train_acc = accuracy_score(y_train, y_pred)
train_auc = roc_auc_score(y_train, y_pred)

# predict to testing
y_pred_test = pipe_baseline.predict(X_test)

test_acc = accuracy_score(y_test, y_pred_test)
test_auc = roc_auc_score(y_test, y_pred_test)
```

## Create modeling pipeline

Using `make_pipeline` so I ensure that centering and scaling is applied
consistently across CV folds, training, and testing data.

``` python
pipe = make_pipeline(
  StandardScaler(),
  XGBClassifier(tree_method = "hist", device = "cuda", updater = "grow_gpu_hist", eval_metric = "auc")
)

pipe
```

## Tune

I’m using Bayesian optimization to tune model hyperparameters.

``` python
hparams = {
  "xgbclassifier__n_estimators": Integer(100, 2000),
  "xgbclassifier__eta": Real(1e-3, 1, "log-uniform"),
  "xgbclassifier__gamma": Real(1e-10, 1.5, "log-uniform"),
  "xgbclassifier__max_depth": Integer(1, 15),
  "xgbclassifier__subsample": Real(0.1, 1, "log-uniform"),
  "xgbclassifier__colsample_bytree": Real(0.01, 0.2, "log-uniform"),
  "xgbclassifier__alpha": Real(0, 0.5, "log-uniform")
}
```

``` python
# get a seed based on the global seed for input into function
cv_seed = np.random.randint(1, 1e6)


searchcv = BayesSearchCV(
  pipe_baseline,
  search_spaces=hparams,
  n_iter = 50,
  cv = 5,
  scoring = "roc_auc",
  random_state=cv_seed
)

np.int = int # have to do this because skopt hasn't totally updated their use of np.int vs int
searchcv.fit(X_train, y_train)
```

## Write tuning results to file

``` python
cv_out_path = "cv_results.pkl"

skdump(searchcv, out_path)
```
