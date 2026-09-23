import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xgboost as xgb
import joblib
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error

# Function used to convert deltaG ---into---> Kd
def calculate_dissociation_constant(deltaG):
    T = 300  # temperature in K
    R = 0.001987  # Universal gas constant in kcal/(mol.K)
    BindingAffinity = np.exp(deltaG / (R * T))
    return BindingAffinity

# ===============================================
# PREDICTING Kd FROM TRAINED ML MODEL
# ===============================================
df = pd.read_csv('test-file.tsv', sep ='\t')

# Separate features (X) and target (y)
X = df.iloc[:, 1:]
y = df.iloc[:, -1].values 

# Load the model
loaded_model =joblib.load('xgboost_model.pkl')
predictions = loaded_model.predict(X)
print(f"\nPrediction probabilities: {predictions}")

# map the binding energy to dissociation constant 
predicted_Kd = calculate_dissociation_constant(predictions)

# map the dissociation constant back to bNAbs in column 0
bNAbs = df.iloc[:, 0].values
results_df = pd.DataFrame({'bNAb': bNAbs,'Actual_DeltaG': y, 'Predicted_DeltaG': predictions, 'Predicted_Dissociation_constant': predicted_Kd})
results_df.to_csv('XGBoost_predicted_Kd.tsv', sep='\t', index=False)
print(results_df)
print("\nPredictions saved to 'predicted_bNAb_results.csv'")

