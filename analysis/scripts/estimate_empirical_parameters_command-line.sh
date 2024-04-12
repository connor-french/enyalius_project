# command line arguments used to estimate empirical parameters for each species
# I copy-pasted into the command line, but I guess this could be run as a bash script

## E. catenatus
python3 estimate_empirical_parameters.py \
--emp_path "../output/empirical_sumstats/catenatus_sumstats.csv" \
--X_train "../output/parameter_estimation/reg_tuning_cat/x_train_cat.csv" \
--X_test "../output/parameter_estimation/reg_tuning_cat/x_test_cat.csv" \
--y_train "../output/parameter_estimation/reg_tuning_cat/y_train_cat.csv" \
--searchcv_path "../output/parameter_estimation/reg_tuning_cat/reg_cv_tuning_cat.pkl" \
--responses "max_k" "Nm" "admix_n_a1" "admix_n_a2" "admix_n_a3" \
--prior_ranges "100,5000" "0.1,5" "5e5,2e6" "5e4,1e6" "5e4,1e6" \
--output_path "../output/parameter_estimation/reg_tuning_cat" 

## E. iheringii
python3 estimate_empirical_parameters.py \
--emp_path "../output/empirical_sumstats/iheringii_sumstats.csv" \
--X_train "../output/parameter_estimation/reg_tuning_ihe/x_train_ihe.csv" \
--X_test "../output/parameter_estimation/reg_tuning_ihe/x_test_ihe.csv" \
--y_train "../output/parameter_estimation/reg_tuning_ihe/y_train_ihe.csv" \
--searchcv_path "../output/parameter_estimation/reg_tuning_ihe/ihe_cv_reg_results.pkl" \
--responses "max_k" "Nm" "admix_n_a1" "admix_n_a2" \
--prior_ranges "100,5000" "0.1,5" "1e6,3e6" "5e4,5e5"  \
--output_path "../output/parameter_estimation/reg_tuning_ihe" 
