# QSRR-Retention-Time-Prediction

QSRR model development for retention time prediction in liquid chromatography.
## Documentation

Read the results of the QSRR model created : [https://narvall018.github.io/QSRR-Retention-Time-Prediction](https://narvall018.github.io/QSRR-Retention-Time-Prediction)

## Installation

```bash
pip install git+https://github.com/narvall018/QSRR-Retention-Time-Prediction.git
```

## Example Usage and Results

```python
from qsrr_predictor import RTPredictor

# List of SMILES to predict
smiles_test = [
    "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirine
    "CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F",  # Celecoxib
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caféine
]

# Initialize and predict
predictor = RTPredictor()
results = predictor.predict(smiles_test)
```

### Output:

```
Nombre de features sélectionnées: 91
Répartition des descripteurs:
- ECFP: 47
- MACCS: 6
- MQN: 1
- RDKit: 37

Conversion des SMILES en molécules...
100%|███████████████████████████████████████████| 3/3 [00:00<00:00, 1959.65it/s]

Calcul des descripteurs...
100%|███████████████████████████████████████████| 3/3 [00:00<00:00, 2950.96it/s]

Résultats des prédictions:
                                     Original_SMILES                                    Canonical_SMILES  Predicted_RT
0                       CC(=O)OC1=CC=CC=C1C(=O)O                          CC(=O)Oc1ccccc1C(=O)O      9.322365
1  CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)... Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)...     10.266563
2                   CN1C=NC2=C1C(=O)N(C(=O)N2C)C                     Cn1c(=O)c2c(ncn2C)n(C)c1=O      4.982333
```

### Results Explanation:

The model predicts retention times in minutes for each input molecule:
- Aspirin: 9.32 minutes
- Celecoxib: 10.27 minutes
- Caffeine: 4.98 minutes

The output includes:
1. Original SMILES input
2. Canonical SMILES (standardized format)
3. Predicted retention time in minutes

## Model Features

The predictor uses multiple types of molecular descriptors:
- ECFP (Extended-Connectivity Fingerprints): 47 features
- MACCS keys: 6 features
- MQN descriptors: 1 feature
- RDKit descriptors: 37 features

## Try Your Own Molecules

You can predict retention times for your own molecules by providing SMILES strings:

```python
your_smiles = [
    "CC(=O)NC1=CC=C(O)C=C1"  # Paracetamol
]

predictor = RTPredictor()
results = predictor.predict(your_smiles)
print(results)
```