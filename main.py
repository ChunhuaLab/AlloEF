import os
import sys
import time
import traceback
import argparse
import pandas as pd
import numpy as np
from preprocessing import load_data, preprocess_data, oversample_data, feature_selection
from modeling import build_models, train_evaluate_model
from evaluation import evaluate_model
from visualization import analyze_shap_values


def predict_single_file(model_path, feature_selector_path, input_file, output_file=None):
    """
    Predict using a single input file

    Args:
        model_path: Path to saved model
        feature_selector_path: Path to saved feature selector
        input_file: Path to input CSV/Excel file
        output_file: Path to output file (optional)
    """
    try:
        print(f"Loading input file: {input_file}")

        # Load input data
        if input_file.endswith('.csv'):
            df_input = pd.read_csv(input_file)
        elif input_file.endswith('.xlsx') or input_file.endswith('.xls'):
            df_input = pd.read_excel(input_file)
        else:
            raise ValueError("Input file must be CSV or Excel format")

        print(f"Loaded {len(df_input)} samples from {input_file}")

        # Load trained model and feature selector
        import joblib
        model = joblib.load(model_path)
        selector = joblib.load(feature_selector_path)

        print("Loaded trained model and feature selector")

        # Preprocess input data (without labels)
        X_input, _, feature_names = preprocess_data(df_input, has_labels=False)

        # Apply feature selection
        X_input_selected = selector.transform(X_input)

        # Make predictions
        y_pred = model.predict(X_input_selected)
        y_pred_proba = model.predict_proba(X_input_selected)

        # Create results dataframe
        results_df = df_input.copy()
        results_df['Predicted_Label'] = y_pred
        results_df['Prediction_Probability_Class_0'] = y_pred_proba[:, 0]
        results_df['Prediction_Probability_Class_1'] = y_pred_proba[:, 1]

        # Save results
        if output_file is None:
            base_name = os.path.splitext(input_file)[0]
            output_file = f"{base_name}_predictions.xlsx"

        if output_file.endswith('.csv'):
            results_df.to_csv(output_file, index=False)
        else:
            results_df.to_excel(output_file, index=False)

        print(f"Predictions saved to: {output_file}")

        # Print summary
        positive_count = np.sum(y_pred == 1)
        total_count = len(y_pred)
        print(f"\nPrediction Summary:")
        print(f"Total samples: {total_count}")
        print(f"Predicted positive: {positive_count}")
        print(f"Predicted negative: {total_count - positive_count}")
        print(f"Positive rate: {positive_count / total_count:.2%}")

        return results_df

    except Exception as e:
        print(f"Error in single file prediction: {str(e)}")
        print(f"Traceback: {traceback.format_exc()}")
        return None


def train_model():
    """
    Train the model using training and independent test data
    """
    start_time = time.time()
    report_path = "training_report.txt"

    try:
        # File path configuration
        train_file = "./train.xlsx"
        independent_test_file = "./independent_test.xlsx"

        print("Loading and preprocessing training data...")
        # Load and preprocess data
        df_train = load_data(train_file)
        X_train, y_train, feature_names = preprocess_data(df_train)

        print("Loading and preprocessing independent test data...")
        df_independent = load_data(independent_test_file)
        X_independent, y_independent, _ = preprocess_data(df_independent)

        print("Applying oversampling...")
        X_train_res, y_train_res = oversample_data(X_train, y_train)

        print("Performing feature selection...")
        X_train_selected, selected_features, selector = feature_selection(
            X_train_res, y_train_res, feature_names
        )

        print("Building models...")
        models = build_models()
        ensemble_model = models["ensemble"]

        print("Training and evaluating model...")
        with open(report_path, "w") as report_file:
            train_evaluate_model(
                ensemble_model,
                X_train_selected,
                y_train_res,
                X_independent,
                y_independent,
                selector,
                report_file
            )

        print("Analyzing SHAP values...")
        os.makedirs("results", exist_ok=True)
        analyze_shap_values(
            ensemble_model,
            X_train_selected,
            selected_features,
            "results/shap"
        )

        # Save trained model and feature selector
        import joblib
        os.makedirs("models", exist_ok=True)
        joblib.dump(ensemble_model, "models/trained_model.pkl")
        joblib.dump(selector, "models/feature_selector.pkl")
        print("Model and feature selector saved to models/ directory")

    except Exception as e:
        with open(report_path, "a") as f:
            f.write(f"Exception: {traceback.format_exc()}\n")
        print(f"Error: {traceback.format_exc()}")

    finally:
        print(f"Total training time: {time.time() - start_time:.2f}s")


def main():
    parser = argparse.ArgumentParser(
        description="Machine Learning Pipeline for Protein Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Train the model
  python main.py --mode train

  # Predict using a single file
  python main.py --mode predict --input ./Input_data/1WQW_A_feature.csv

  # Predict with custom output file
  python main.py --mode predict --input ./Input_data/1WQW_A_feature.csv --output ./results/predictions.xlsx
        """
    )

    parser.add_argument(
        "--mode",
        choices=["train", "predict"],
        required=True,
        help="Mode: 'train' to train the model, 'predict' to predict using existing model"
    )

    parser.add_argument(
        "--input",
        type=str,
        help="Input file path for prediction (CSV or Excel format)"
    )

    parser.add_argument(
        "--output",
        type=str,
        help="Output file path for predictions (optional)"
    )

    parser.add_argument(
        "--model",
        type=str,
        default="models/trained_model.pkl",
        help="Path to trained model file (default: models/trained_model.pkl)"
    )

    parser.add_argument(
        "--selector",
        type=str,
        default="models/feature_selector.pkl",
        help="Path to feature selector file (default: models/feature_selector.pkl)"
    )

    # Handle legacy command line format: python main.py input_file
    if len(sys.argv) == 2 and not sys.argv[1].startswith('--'):
        input_file = sys.argv[1]
        print(f"Legacy format detected. Predicting for file: {input_file}")

        # Check if model files exist
        model_path = "models/trained_model.pkl"
        selector_path = "models/feature_selector.pkl"

        if not os.path.exists(model_path) or not os.path.exists(selector_path):
            print("Error: Trained model not found. Please train the model first:")
            print("python main.py --mode train")
            return

        predict_single_file(model_path, selector_path, input_file)
        return

    args = parser.parse_args()

    if args.mode == "train":
        print("Starting model training...")
        train_model()

    elif args.mode == "predict":
        if not args.input:
            print("Error: --input argument is required for prediction mode")
            parser.print_help()
            return

        if not os.path.exists(args.input):
            print(f"Error: Input file does not exist: {args.input}")
            return

        # Check if model files exist
        if not os.path.exists(args.model) or not os.path.exists(args.selector):
            print("Error: Trained model not found. Please train the model first:")
            print("python main.py --mode train")
            return

        print(f"Starting prediction for file: {args.input}")
        predict_single_file(args.model, args.selector, args.input, args.output)


if __name__ == "__main__":
    main()