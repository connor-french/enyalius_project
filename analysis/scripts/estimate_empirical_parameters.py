import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotnine import ggplot, aes, geom_point, geom_segment, geom_hline, labs, theme_bw, theme, element_blank, lims
import pickle
from mapie.regression import MapieRegressor
from skopt import load as skload
from os import path

# the scope has also widened to include writing out test predictions, since I didn't include them with my original simulation scripts

def read_empirical(path):

  empirical_df = pd.read_csv(path)

  # only keep relevant summary statistics
  empirical_df = empirical_df.loc[:, empirical_df.columns.str.startswith(('pi_a', 'pi_sd_a', 'taj', 'sfs', 'dxy', 'ibd'))]

  # remove columns that contain "fst" or "h3"
  empirical_df = empirical_df.loc[:, ~empirical_df.columns.str.contains('fst|h3')]

  # remove dxy from column names
  empirical_df = empirical_df.rename(columns=lambda x: x.replace('_dxy', '') if x.startswith('ibd') else x)
  
  return empirical_df



def plot_pred_int(pred_int, param, prior_range):
  
  y1 = pred_int["lower"]
  y2 = pred_int["upper"]


  ylims = prior_range.copy()
  # in case interval is outside of prior range
  if y1[0] < prior_range[0]:
    ylims[0] = y1[0]

  if y2[0] > prior_range[1]:
    ylims[1] = y2[0]


  emp_predplot_df = pd.DataFrame(
    {
      "x1": 1,
      "x2": 1,
      "y1": y1,
      "y2": y2,
      "center": pred_int["center"]
    }
    
  )

  p = (
      ggplot(emp_predplot_df) +
        geom_segment(aes(x = "x1", xend = "x2", y = "y1", yend = "y2"), size = 3) +
        geom_point(aes(x = "x1", y = "center"), size = 8) +
        geom_hline(yintercept = prior_range[0], linetype = "dashed") +
        geom_hline(yintercept = prior_range[1], linetype = "dashed") +
        lims(y = (ylims[0], ylims[1])) +
        labs(y = f"{param}", title = f"95% prediction interval of {param}") +
        theme_bw() +
        theme(axis_text_x = element_blank(), 
        axis_title_x=element_blank())
  )
  
  return p


# function to perform the whole predictive interval workflow
def calc_pred_interval(x_train, y_train, emp_data, estimator, response, prior_range):
  
  # fit CV folds
  mapie_cv = MapieRegressor(estimator, cv = 5, n_jobs = 5, agg_function = "median")
  mapie_cv.fit(x_train, y_train[response])

  # predict to empirical data
  emp_pred = mapie_cv.predict(emp_data, alpha = [0.05])

  # convert to a DataFrame
  emp_pred = pd.DataFrame(
    {
      "center": emp_pred[0],
      "lower": emp_pred[1][0][0],
      "upper": emp_pred[1][0][1],
      "response": response
    }
  )

  # plot interval
  pred_plot = plot_pred_int(pred_int = emp_pred, param = response, prior_range = prior_range)

  return mapie_cv, emp_pred, pred_plot


def main(emp_path, X_train, X_test, y_train, searchcv_path, prior_ranges, responses, output_path):
  
  emp = read_empirical(emp_path)
  X_train = pd.read_csv(X_train)
  X_test = pd.read_csv(X_test)
  y_train = pd.read_csv(y_train)
  searchcv = skload(searchcv_path)
  best_pipeline = searchcv.best_estimator_
  best_pipeline.fit(X_train, y_train)

  # write out predictions to the test set as a csv
  test_pred = best_pipeline.predict(X_test)

  #pd.DataFrame(test_pred).to_csv(path.join(output_path, "test_preds.csv"), index=False)
  # write out using numpy instead
  np.savetxt(path.join(output_path, "test_preds.csv"), test_pred, delimiter = ",")


  for i, response in enumerate(responses):
    prior_range = prior_ranges[i]
    mapie, emp_pred, pred_plot = calc_pred_interval(x_train=X_train, y_train=y_train, emp_data=emp, estimator=best_pipeline, response=response, prior_range=prior_range)
    pred_plot.save(path.join(output_path, f"pred_plot_{response}.png"))
    with open(path.join(output_path, f"mapie_{response}.pkl"), "wb") as f:
      pickle.dump(mapie, f)
    with open(path.join(output_path, f"emp_pred_{response}.csv"), "wb") as f:
      emp_pred.to_csv(f, index=False)

  cv_results = pd.DataFrame(searchcv.cv_results_)

  # write out cross-validation results
  cv_results.to_csv(path.join(output_path, "cv_results.csv"))

  # plot rank test score vs hyperparameters
  param_vars = [var for var in cv_results.columns if var.startswith("param_xgbregressor")]

  for var in param_vars:
      plt.scatter(cv_results[var], cv_results["rank_test_score"])
      plt.xlabel(var)
      plt.ylabel("rank_test_score")
      plt.savefig(path.join(output_path, f"ranktest_{var}.png"))
      plt.close()

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("--emp_path", help="Path to the empirical data file")
  parser.add_argument("--X_train", help="Path to the X_train data file")
  parser.add_argument("--X_test", help="Path to the X_test data file")
  parser.add_argument("--y_train", help="Path to the y_train data file")
  parser.add_argument("--searchcv_path", help="Path to the searchcv data file")
  parser.add_argument("--responses", nargs='+', help="List of responses", type=str)
  parser.add_argument("--prior_ranges", nargs='+', help="List of prior ranges")
  parser.add_argument("--output_path", help="Path to the output directory")
  args = parser.parse_args()
  prior_ranges = [list(map(float, range_str.split(','))) for range_str in args.prior_ranges]

  main(args.emp_path, args.X_train, args.X_test, args.y_train, args.searchcv_path, prior_ranges, args.responses, args.output_path)


