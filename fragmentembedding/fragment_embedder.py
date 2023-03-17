import os
import joblib
import numpy as np
import onnxruntime as rt
import sys


CHUNKSIZE = 1024

root = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root)


class FragmentEmbedder(object):
    def __init__(self, models_dir=None):
        if models_dir is None:
            models_dir = os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "..", "results"
            )
        self.models_dir = os.path.abspath(models_dir)
        self.morgan_desc = joblib.load(
            os.path.join(models_dir, "morgan_descriptor.joblib")
        )
        self.physchem_desc = joblib.load(
            os.path.join(models_dir, "physchem_descriptor.joblib")
        )

    def _chunker(self, l, n):
        for i in range(0, len(l), n):
            yield l[i : i + n]

    def encoder_inference(self, X):
        sess = rt.InferenceSession(os.path.join(self.models_dir, "encoder_model.onnx"))
        input_name = sess.get_inputs()[0].name
        output_name = sess.get_outputs()[0].name
        output_data = sess.run(
            [output_name], {input_name: np.array(X, dtype=np.float32)}
        )
        Y = np.array(output_data[0])
        return Y

    def transform(self, smiles):
        X = None
        for smiles_chunk in self._chunker(smiles, CHUNKSIZE):
            X_0 = self.morgan_desc.transform(smiles_chunk)
            X_1 = self.physchem_desc.transform(smiles_chunk)
            X_i = np.hstack([X_0, X_1])
            X_o = self.encoder_inference(X_i)
            if X is None:
                X = X_o
            else:
                X = np.vstack([X, X_o])
        return X
