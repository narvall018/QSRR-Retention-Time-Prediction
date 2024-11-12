# QSRR-Retention-Time-Prediction 🔬

> QSRR model for predicting retention times in liquid chromatography.

[![Python Version](https://img.shields.io/badge/python-3.6%2B-blue)]()
[![License](https://img.shields.io/badge/license-MIT-green.svg)]()

## 📖 Documentation 

Detailed documentation available at: [https://narvall018.github.io/QSRR-Retention-Time-Prediction](https://narvall018.github.io/QSRR-Retention-Time-Prediction)

## 🚀 Installation

```bash
pip install git+https://github.com/narvall018/QSRR-Retention-Time-Prediction.git
```

### Dependencies
- RDKit
- NumPy  
- Pandas
- scikit-learn
- joblib
- tqdm

## 💡 Usage

Initialize the predictor

```python
from qsrr_predictor import RTPredictor

# Example molecules
smiles_test = [
   "CC(=O)OC1=CC=CC=C1C(=O)O",         # Aspirin
   "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",     # Caffeine
   "CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F"  # Celecoxib
]

# Initialize and predict
predictor = RTPredictor()
results = predictor.predict(smiles_test)
print(results)
```

# Example output:
Initialization successful!
Number of selected features: 91
Distribution of descriptors:
- ECFP: 47
- MACCS: 6
- MQN: 1
- RDKit: 37

Computing molecular descriptors...
100%|████████████████████████████████████| 3/3 [00:00<00:00, 336.66it/s]


| Index | Original_SMILES | Canonical_SMILES | Predicted_RT |
|-------|----------------|------------------|--------------|
| 0 | CC(=O)OC1=CC=CC=C1C(=O)O | CC(=O)Oc1ccccc1C(=O)O | 9.26 |
| 1 | CN1C=NC2=C1C(=O)N(C(=O)N2C)C | Cn1c(=O)c2c(ncn2C)n(C) | 5.08 |
| 2 | CC1=CC=C(C=C1)C2=CC(=NN2C3=... | Cc1ccc(-c2cc(C(F)(F)... | 10.33 |

## 🧬 Molecular Descriptors

The model uses four types of descriptors:

| Type | Description | Count |
|------|-------------|-------|
| ECFP | Circular molecular fingerprints (radius 3) | 47 |
| MACCS | Predefined structural keys | 6 |
| MQN | Atom/bond based descriptors | 1 |
| RDKit | Physicochemical properties | 37 |

## 📊 Results Format

The output DataFrame contains:
- Original_SMILES: Input SMILES
- Canonical_SMILES: Standardized SMILES
- Predicted_RT: Predicted retention time (minutes)

## ⚠️ Important Notes

1. Ensure input SMILES are valid
2. Predictions are optimal for molecules similar to the training set
3. Times are based on specific chromatographic conditions

## 📫 Support & Contact

For any questions or issues:
- Open an issue on GitHub
- Contact me : julien.sade@u-pec.fr
- Check the [documentation](https://narvall018.github.io/QSRR-Retention-Time-Prediction)

## 📄 License

This project is licensed under MIT - see the [LICENSE](LICENSE) file for details.