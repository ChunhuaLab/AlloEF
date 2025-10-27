from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import make_scorer, matthews_corrcoef


def build_models():
    xgb = XGBClassifier(
        tree_method='hist', random_state=42,
        colsample_bytree=1.0, learning_rate=0.083,
        max_depth=10, n_estimators=200, subsample=1.0
    )

    lgbm = LGBMClassifier(
        is_unbalance=True, random_state=42,
        colsample_bytree=0.72, learning_rate=0.1,
        max_depth=10, n_estimators=200, subsample=0.82
    )

    rf = RandomForestClassifier(
        random_state=42, max_depth=14, min_samples_leaf=1,
        min_samples_split=2, n_estimators=198
    )

    return {
        "xgb": xgb,
        "lgbm": lgbm,
        "rf": rf,
        "ensemble": VotingClassifier(
            estimators=[('xgb', xgb), ('lgbm', lgbm), ('rf', rf)],
            voting='soft'
        )
    }


def train_evaluate_model(model, X_train, y_train,
                         X_independent, y_independent, selector, report_file):

    mcc_scorer = make_scorer(matthews_corrcoef)
    cv_scores = cross_val_score(model, X_train, y_train, cv=5, scoring=mcc_scorer)

    report_file.write(f"Cross-validation MCC: {np.mean(cv_scores):.4f}\n")
    report_file.write(f"Scores: {cv_scores}\n")

    model.fit(X_train, y_train)

    X_independent_selected = X_independent[:, selector.support_]
    evaluate_model(model, X_independent_selected, y_independent, report_file, "Independent Test Set")