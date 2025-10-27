from sklearn.metrics import (accuracy_score, roc_auc_score, matthews_corrcoef,
                             f1_score, recall_score, precision_score,
                             confusion_matrix, classification_report,
                             average_precision_score)


def evaluate_model(model, X, y, file_obj, set_name="Dataset"):
    y_pred = model.predict(X)
    y_prob = model.predict_proba(X)[:, 1]

    metrics = {
        "Accuracy": accuracy_score(y, y_pred),
        "AUC": roc_auc_score(y, y_prob),
        "AUPRC": average_precision_score(y, y_prob),
        "MCC": matthews_corrcoef(y, y_pred),
        "F1": f1_score(y, y_pred),
        "Recall": recall_score(y, y_pred),
        "Precision": precision_score(y, y_pred)
    }

    tn, fp, fn, tp = confusion_matrix(y, y_pred).ravel()
    metrics["Specificity"] = tn / (tn + fp)

    file_obj.write(f"\n===== {set_name} Evaluation =====\n")
    for name, value in metrics.items():
        file_obj.write(f"{name}: {value:.4f}\n")

    file_obj.write("\nConfusion Matrix:\n")
    file_obj.write(str(confusion_matrix(y, y_pred)))

    file_obj.write("\n\nClassification Report:\n")
    file_obj.write(classification_report(y, y_pred))