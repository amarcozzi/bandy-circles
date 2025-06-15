import pandas as pd
import numpy as np
from pathlib import Path
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go

# --- 1. Load and Prepare Data ---
data_path = Path("/Volumes/T7 Shield/bandy-circles/grid-search-combined")
quantity = "heat_flux"
df = pd.read_csv(data_path / f"ideal_curvatures_{quantity}.csv")

# Create scaled curvature
df['line_length'] = df['circle_radius']
df['curvature'] = df['curvature'] * df['line_length']

df = df[df['wind_speed'] <= 3.1]  # Filter for wind speed <= 3 m/s

print(f"Filtered data shape: {df.shape}")
if df.empty:
    raise ValueError("Filtering resulted in an empty DataFrame. Check filter conditions.")

# --- 2. EDA (Correlation) ---
print("Correlation Matrix:")
# Choose 'curvature' or 'scaled_curvature'
target_var = 'curvature'
predictors = ['wind_speed', 'treatment_fuel_height', 'circle_radius']
print(df[[target_var] + predictors].corr())

# --- 3. Regression Modeling ---

# Model 1: Linear
formula1 = f"{target_var} ~ wind_speed + treatment_fuel_height + circle_radius"
model1 = smf.ols(formula1, data=df).fit()
print("\n--- Model 1: Linear ---")
print(model1.summary())

# Model 2: Linear + Interactions
# Or just 2-way:
formula2 = f"{target_var} ~ (wind_speed + treatment_fuel_height + circle_radius)**2"
model2 = smf.ols(formula2, data=df).fit()
print("\n--- Model 2: Linear + Interactions (up to 3-way) ---")
print(model2.summary())

# --- 4. Save the Model ---
import pickle

# Save the model with all its details
with open(data_path / 'curvature_model2.pkl', 'wb') as f:
    pickle.dump(model2, f)

# Save just the coefficients as a JSON file for easier access
import json

coef_dict = {
    'intercept': model2.params[0],
    'wind_speed': model2.params[1],
    'treatment_fuel_height': model2.params[2],
    'circle_radius': model2.params[3],
    'wind_speed:treatment_fuel_height': model2.params[4],
    'wind_speed:circle_radius': model2.params[5],
    'treatment_fuel_height:circle_radius': model2.params[6]
}

with open(data_path / 'model2_coefficients.json', 'w') as f:
    json.dump(coef_dict, f, indent=4)

print(f"Model saved to {data_path / 'curvature_model2.pkl'}")
print(f"Coefficients saved to {data_path / 'model2_coefficients.json'}")
print("Model coefficients:")
for key, value in coef_dict.items():
    print(f"{key}: {value:.4f}")
