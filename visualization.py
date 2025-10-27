import os
import shap
import matplotlib.pyplot as plt


def analyze_shap_values(model, X, feature_names, output_dir):
    os.makedirs(output_dir, exist_ok=True)


    explainer = shap.Explainer(model.named_estimators['xgb'])
    shap_values = explainer(X)


    plt.figure()
    shap.summary_plot(shap_values, X, feature_names=feature_names, plot_type="bar", show=False)
    plt.savefig(os.path.join(output_dir, "shap_summary_bar.png"), bbox_inches='tight')
    plt.close()


    plt.figure()
    shap.summary_plot(shap_values, X, feature_names=feature_names, show=False)
    plt.savefig(os.path.join(output_dir, "shap_summary.png"), bbox_inches='tight')
    plt.close()


    plt.figure()
    shap.plots.waterfall(shap_values[0], show=False)
    plt.savefig(os.path.join(output_dir, "shap_waterfall.png"), bbox_inches='tight')
    plt.close()
