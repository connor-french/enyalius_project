
# script for catenatus classification
plot_feature_importance.py \
--species "cat" \
--model "Classification" \
--X_train_path "../../analysis/output/model_classification/results_cat/X_train_cat.csv" \
--y_train_path "../../analysis/output/model_classification/results_cat/y_train_cat.csv" \
--searchcv_path "../../analysis/output/model_classification/results_cat/cv_class_results_cat.pkl" \
--output_folder "../figures/feature_importance"

# script for catenatus regression
plot_feature_importance.py \
--species "cat" \
--model "Regression" \
--X_train_path "../../analysis/output/parameter_estimation/reg_tuning_cat/x_train_cat.csv" \
--y_train_path "../../analysis/output/parameter_estimation/reg_tuning_cat/y_train_cat.csv" \
--searchcv_path "../../analysis/output/parameter_estimation/reg_tuning_cat/reg_cv_tuning_cat.pkl" \
--output_folder "../figures/feature_importance"


# script for iheringii classification
plot_feature_importance.py \
--species "ihe" \
--model "Classification" \
--X_train_path "../../analysis/output/model_classification/results_ihe/X_train_ihe.csv" \
--y_train_path "../../analysis/output/model_classification/results_ihe/y_train_ihe.csv" \
--searchcv_path "../../analysis/output/model_classification/results_ihe/ihe_class_cv_results.pkl" \
--output_folder "../figures/feature_importance"


# script for iheringii regression
plot_feature_importance.py \
--species "ihe" \
--model "Regression" \
--X_train_path "../../analysis/output/parameter_estimation/reg_tuning_ihe/X_train_ihe.csv" \
--y_train_path "../../analysis/output/parameter_estimation/reg_tuning_ihe/y_train_ihe.csv" \
--searchcv_path "../../analysis/output/parameter_estimation/reg_tuning_ihe/ihe_cv_reg_results.pkl" \
--output_folder "../figures/feature_importance"