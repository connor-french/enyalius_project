# command line arguments for write_model_classification_preds.py
# "Usage: python3 write_model_classification_preds.py --species --cv_path --sims_dir --out_dir --seed"

## iheringii
python3 write_model_classification_preds.py \
--species "ihe" \
--cv_path "../output/model_classification/results_ihe/ihe_class_cv_results.pkl" \
--sims_dir "../output/simulations" \
--out_dir "../output/model_classification/results_ihe" \
--seed 1112

# catenatus
python3 write_model_classification_preds.py \
--species "cat" \
--cv_path "../output/model_classification/results_cat/cv_class_results_cat.pkl" \
--sims_dir "../output/simulations" \
--out_dir "../output/model_classification/results_cat" \
--seed 21212121