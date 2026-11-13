import argparse, json, glob, os
import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.metrics import log_loss

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--integrated-tsvs", required=True, help="Путь к integrated_variants.tsv development кейсов (можно с wildcard *)")
    p.add_argument("--labels-tsv", required=True, help="development_labels.tsv")
    p.add_argument("--out-model", required=True, help="Куда сохранить модель (напр., model.json)")
    p.add_argument("--out-threshold", required=True, help="Куда сохранить порог abstention (напр., threshold.txt)")
    args = p.parse_args()

    paths = glob.glob(args.integrated_tsvs)
    dfs = [pd.read_csv(p, sep="\t", dtype=str) for p in paths if os.path.getsize(p) > 0]
    df = pd.concat(dfs, ignore_index=True)

    labels = pd.read_csv(args.labels_tsv, sep="\t", dtype=str)
    df['y'] = 0
    
    for _, row in labels.iterrows():
        chrom = str(row['chr']).replace('chr', '').replace('CHR', '')
        mask = (
            (df['case_id'] == str(row['case_id'])) &
            (df['chrom'].astype(str).str.replace('chr', '').str.replace('CHR', '') == chrom) &
            (df['pos'].astype(str) == str(row['pos'])) &
            (df['ref'].str.upper() == str(row['ref']).upper()) &
            (df['alt'].str.upper() == str(row['alt']).upper())
        )
        df.loc[mask, 'y'] = 1

    feature_cols = ["feat_functional", "feat_clinical", "feat_population_rarity", "feat_phen2gene", "feat_lirical_lr", "feat_zygosity_compat"]
    X = df[feature_cols].fillna(0.0)
    y = df['y']
    case_ids = df['case_id']

    unique_cases = case_ids.unique()
    max_probs_negatives = []
    
    print(f"Запуск LOCO-CV на {len(unique_cases)} кейсах...")
    for test_case in unique_cases:
        train_mask = case_ids != test_case
        test_mask = case_ids == test_case
        
        X_train, y_train = X[train_mask], y[train_mask]
        X_test, y_test = X[test_mask], y[test_mask]
        
        scale_pos_weight = (y_train == 0).sum() / max((y_train == 1).sum(), 1)
        
        model = xgb.XGBClassifier(
            objective='binary:logistic', eval_metric='logloss', max_depth=3,
            learning_rate=0.05, min_child_weight=5, gamma=0.2,
            subsample=0.8, colsample_bytree=0.8, n_estimators=200,
            scale_pos_weight=scale_pos_weight, random_state=42
        )
        model.fit(X_train, y_train, verbose=False)
        
        probs = model.predict_proba(X_test)[:, 1]
        
        max_probs_negatives.append(probs.max())

    threshold = float(np.percentile(max_probs_negatives, 75))
    print(f"Калиброванный порог abstention: {threshold:.4f}")

    final_scale_pos_weight = (y == 0).sum() / max((y == 1).sum(), 1)
    final_model = xgb.XGBClassifier(
        objective='binary:logistic', eval_metric='logloss', max_depth=3,
        learning_rate=0.05, min_child_weight=5, gamma=0.2,
        subsample=0.8, colsample_bytree=0.8, n_estimators=200,
        scale_pos_weight=final_scale_pos_weight, random_state=42
    )
    final_model.fit(X, y, verbose=False)

    final_model.save_model(args.out_model)
    with open(args.out_threshold, 'w') as f:
        f.write(f"{threshold:.4f}")
    
    print(f"model saved: {args.out_model}")
    print(f"threshold saved: {args.out_threshold}")

if __name__ == "__main__":
    main()