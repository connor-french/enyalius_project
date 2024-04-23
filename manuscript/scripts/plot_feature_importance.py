import argparse
import os
import pandas as pd
from xgboost import XGBClassifier, XGBRegressor, plot_importance
from skopt import BayesSearchCV
from skopt import load as skload
from matplotlib import pyplot

# script to plot feature importance for model classification and regression
def main():
  # Parse command line arguments
  parser = argparse.ArgumentParser()
  parser.add_argument('--species', type=str, help='Which species was inferred')
  parser.add_argument('--model', type=str, help='Which model was used. "classification" or "regression"')
  #parser.add_argument('--param_list', nargs='+', help='List of responses in multioutput regression', type=str)
  parser.add_argument('--X_train_path', type=str, help='Path to training set predictors (csv file)')
  parser.add_argument('--y_train_path', type=str, help='Path to training set for responses (csv file)')
  parser.add_argument('--searchcv_path', type=str, help='Path to model cross-validation object')
  parser.add_argument('--output_folder', type=str, help='Path to output folder')
  args = parser.parse_args()

  X_train = pd.read_csv(args.X_train_path)
  y_train = pd.read_csv(args.y_train_path)

  searchcv = skload(args.searchcv_path)
  best_pipeline = searchcv.best_estimator_
  best_pipeline.fit(X_train, y_train)

  if args.species == "cat":
    species_name = "E. catenatus"
    bar_col = "#CCA200"
  elif args.species == "ihe":
    species_name = "E. iheringii"
    bar_col = "#7CCE00"
  elif args.species == "per":
    species_name = "E. perditus"
    bar_col = "darkgreen"
  
  feature_names = list(X_train.columns)
  feature_importances = best_pipeline._final_estimator.feature_importances_
  
  # Sort feature importances in descending order
  sorted_indices = feature_importances.argsort()[::-1]
  sorted_importances = feature_importances[sorted_indices]
  sorted_names = [feature_names[i] for i in sorted_indices]

  # Plot the top 10 feature importances
  pyplot.barh(sorted_names[:10][::-1], sorted_importances[:10][::-1], color=bar_col)
  pyplot.xlabel('Importance (total gain)')
  pyplot.ylabel('Variable')
  pyplot.title(rf'{args.model} variable importance for ${species_name}$')
  pyplot.tight_layout()


  pyplot.savefig(os.path.join(args.output_folder, f'feature_importance_{args.model}_{args.species}.png'), dpi=600)


if __name__ == "__main__":
  main()
