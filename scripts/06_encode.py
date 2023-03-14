import os
import sys
import joblib
import numpy as np

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "fragmentembedding"))

from encoder import start_encoder

h5_data = os.path.join(root, "..", "data", "data_for_training.h5")

# parameters found in a 50 trials experiment
skip_optuna_with_params = {
    "activation": "swish",
    "layers_num": 4,
    "dropout_prob": 0.1,
    "learning_rate": 0.0004056362326296727,
    "epsilon": 4.743926632146026e-07,
}
# skip_optuna_with_params = None

feature_extractor, hyperparams = start_encoder(
    h5_data,
    trials=100,
    epochs=100,
    timeout=3600 * 24,
    skip_optuna_with_params=skip_optuna_with_params,
)

print("Saving feature extractor")

feature_extractor.save(os.path.join(root, "..", "results", "encoder_model"))
joblib.dump(
    hyperparams, os.path.join(root, "..", "results", "encoder_model_hyperparams.joblib")
)
