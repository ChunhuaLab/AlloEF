import pandas as pd
import numpy as np
import string
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import MinMaxScaler, LabelEncoder
from tensorflow.keras.utils import to_categorical
from boruta import BorutaPy
from sklearn.ensemble import RandomForestClassifier
from imblearn.over_sampling import SVMSMOTE

def load_data(file_path):
    return pd.read_excel(file_path)


def preprocess_data(df):
    # Ensure STRUCTURE column is of string type
    df['STRUCTURE'] = df['STRUCTURE'].astype(str)

    # Extract the first letter
    def extract_first_letter(structure):
        first_letter = structure[0]
        return first_letter if first_letter in string.ascii_letters else "Other"

    df['FirstLetter'] = df['STRUCTURE'].apply(extract_first_letter)

    # Label encode and one-hot encode the first letter
    label_encoder = LabelEncoder()
    encoded_letters = label_encoder.fit_transform(df['FirstLetter'])
    one_hot_encoded_letters = to_categorical(encoded_letters)
    one_hot_columns = label_encoder.classes_
    one_hot_df = pd.DataFrame(one_hot_encoded_letters, columns=one_hot_columns)
    one_hot_df.index = df.index

    # Select numerical columns
    non_numerical_columns = ['ProteinID', 'label', 'FirstLetter', 'STRUCTURE', 'x', 'y', 'z', 'ResidueInfo',
                             'Residue Number']
    numerical_features = df.columns.difference(non_numerical_columns)

    def process_cell(cell):
        if isinstance(cell, str) and "," in cell:
            numbers = cell.split(',')
            return [float(number.strip()) for number in numbers]
        return [float(cell)]

    # Process columns with commas
    comma_columns = [col for col in numerical_features if df[col].astype(str).str.contains(',').any()]
    for col in comma_columns:
        df[col] = df[col].astype(str).apply(process_cell)
        df[col] = df[col].apply(lambda x: sum(x) / len(x) if len(x) > 1 else x[0])

    # Fill missing values with mean in features
    imputer = SimpleImputer(strategy='mean')
    df[numerical_features] = imputer.fit_transform(df[numerical_features])

    # Normalize numerical features
    scaler = MinMaxScaler()
    normalized_features = scaler.fit_transform(df[numerical_features])
    normalized_df = pd.DataFrame(normalized_features, columns=numerical_features)
    normalized_df.index = df.index

    # Combine the DataFrame
    final_df = pd.concat([df[['ProteinID']], df['label'].fillna(0), one_hot_df, normalized_df], axis=1)

    X = final_df.drop(columns=["label", "ProteinID"])
    y = df["label"].fillna(0)

    return X, y, final_df.columns.drop(["label", "ProteinID"])

def oversample_data(X_train, y_train, random_state=42):
    svmsmote = SVMSMOTE(random_state=random_state, n_jobs=-1)
    X_train_res, y_train_res = svmsmote.fit_resample(X_train, y_train)
    return X_train_res, y_train_res

def feature_selection(X_train_res, y_train_res, feature_names):
    # 定义随机森林分类器用于 Boruta
    rf = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1)
    feat_selector = BorutaPy(rf, n_estimators='auto', random_state=42, max_iter=100)

    # 确保 X_train_res 是一个 NumPy 数组
    if isinstance(X_train_res, pd.DataFrame):
        X_train_res = X_train_res.values

    # 训练 Boruta 特征选择器
    feat_selector.fit(X_train_res, y_train_res)

    # 选出所有重要的特征
    selected_mask = feat_selector.support_

    # 使用掩码选择实际特征
    selected_features = [feature_names[i] for i in range(len(feature_names)) if selected_mask[i]]

    # 创建选择后的训练集
    X_train_selected = X_train_res[:, selected_mask]

    return X_train_selected, selected_features, feat_selector
