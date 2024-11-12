import os
import joblib
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, rdMolDescriptors, Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from tqdm import tqdm
import warnings

warnings.filterwarnings('ignore', category=FutureWarning)

class RTPredictor:
    """
    Prédicteur de temps de rétention pour les composés chimiques.
    
    Cette classe utilise un modèle de régression pré-entraîné pour prédire 
    les temps de rétention à partir de descripteurs moléculaires calculés via RDKit.
    
    Attributs:
        model : Le modèle de régression chargé
        selected_features (list): Liste des descripteurs sélectionnés
        needed_bits (list): Indices des bits ECFP nécessaires
        needed_maccs (list): Indices des clés MACCS nécessaires
        needed_mqn (list): Indices des descripteurs MQN nécessaires
        needed_rdkit (list): Noms des descripteurs RDKit nécessaires
    
    Example:
        >>> predictor = RTPredictor()
        >>> smiles_list = ['CCO', 'CCC']
        >>> results = predictor.predict(smiles_list)
        >>> print(results)
    """
    
    def __init__(self):
        """
        Initialise le prédicteur en chargeant le modèle et les descripteurs.
        
        Le constructeur charge le modèle LASSO et la liste des descripteurs depuis 
        les fichiers .joblib stockés dans le dossier 'models'.
        
        Raises:
            Exception: Si le chargement des fichiers modèles échoue
        """
        current_dir = os.path.dirname(os.path.abspath(__file__))
        model_path = os.path.join(current_dir, 'models', 'ridge_all_descriptors_LASSO.joblib')
        features_path = os.path.join(current_dir, 'models', 'features_all_descriptors_LASSO.joblib')
        
        try:
            self.model = joblib.load(model_path)
            features = joblib.load(features_path)
            self.selected_features = [feat.replace('_x', '') for feat in features]
            
            self.needed_bits = [int(f.split('_')[1]) for f in self.selected_features if f.startswith('Bit_')]
            self.needed_maccs = [int(f.split('_')[1]) for f in self.selected_features if f.startswith('MACCS_')]
            self.needed_mqn = [int(f.split('_')[1]) for f in self.selected_features if f.startswith('MQN_')]
            self.needed_rdkit = [f for f in self.selected_features if not any(f.startswith(p) for p in ['Bit_', 'MACCS_', 'MQN_'])]

            print("Initialisation réussie!")
            print("Nombre de features sélectionnées:", len(self.selected_features))
            print("Répartition des descripteurs:")
            print(f"- ECFP: {len(self.needed_bits)}")
            print(f"- MACCS: {len(self.needed_maccs)}")
            print(f"- MQN: {len(self.needed_mqn)}")
            print(f"- RDKit: {len(self.needed_rdkit)}")

        except Exception as e:
            raise Exception(f"Erreur lors du chargement: {str(e)}")

    def prepare_molecules(self, smiles_list):
        """
        Prépare les molécules à partir d'une liste de SMILES.
        
        Args:
            smiles_list (list): Liste de chaînes SMILES à traiter
            
        Returns:
            pandas.DataFrame: DataFrame contenant :
                - original_smiles: SMILES d'origine
                - canonical_smiles: SMILES canonique
                - mol: Objets moléculaires RDKit
        
        Notes:
            Les SMILES invalides sont ignorés silencieusement
        """
        mols = []
        valid_smiles = []
        invalid_smiles = []

        for i, smiles in enumerate(smiles_list):
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    canon_smiles = Chem.MolToSmiles(mol, canonical=True)
                    mols.append(mol)
                    valid_smiles.append(canon_smiles)
                else:
                    invalid_smiles.append((i, smiles))
            except:
                invalid_smiles.append((i, smiles))

        return pd.DataFrame({
            'original_smiles': smiles_list,
            'canonical_smiles': valid_smiles,
            'mol': mols
        })

    def compute_descriptors(self, smiles_list):
        """
        Calcule tous les descripteurs moléculaires nécessaires.
        
        Cette fonction calcule quatre types de descripteurs :
        - ECFP (Extended Connectivity Fingerprints)
        - MACCS (Molecular ACCess System) keys
        - MQN (Molecular Quantum Numbers)
        - Descripteurs RDKit standards
        
        Args:
            smiles_list (list): Liste de chaînes SMILES
            
        Returns:
            tuple: (descriptors_df, canonical_smiles) où
                - descriptors_df est un DataFrame des descripteurs calculés
                - canonical_smiles est la liste des SMILES canoniques
        """
        mol_df = self.prepare_molecules(smiles_list)
        all_descriptors = []
        
        print("\nCalcul des descripteurs moléculaires...")
        for mol in tqdm(mol_df['mol']):
            mol_descriptors = {}
            
            # ECFP
            if self.needed_bits:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=1024)
                fp_bits = list(fp.GetOnBits())
                for i in self.needed_bits:
                    mol_descriptors[f'Bit_{i}'] = 1 if i in fp_bits else 0
            
            # MACCS
            if self.needed_maccs:
                maccs = MACCSkeys.GenMACCSKeys(mol)
                for i in self.needed_maccs:
                    mol_descriptors[f'MACCS_{i}'] = maccs[i]
            
            # MQN
            if self.needed_mqn:
                mqn = rdMolDescriptors.MQNs_(mol)
                for i, idx in enumerate(self.needed_mqn, 1):
                    mol_descriptors[f'MQN_{idx}'] = mqn[idx-1]
            
            # RDKit
            if self.needed_rdkit:
                calc = MoleculeDescriptors.MolecularDescriptorCalculator(self.needed_rdkit)
                rdkit_desc = calc.CalcDescriptors(mol)
                for name, value in zip(self.needed_rdkit, rdkit_desc):
                    mol_descriptors[name] = value
            
            all_descriptors.append(mol_descriptors)
        
        descriptors_df = pd.DataFrame(all_descriptors)
        descriptors_df = descriptors_df[self.selected_features]
        
        return descriptors_df, mol_df['canonical_smiles']

    def predict(self, smiles_list):
        """
        Prédit les temps de rétention pour une liste de SMILES.
        
        Cette méthode effectue le pipeline complet de prédiction :
        1. Préparation des molécules
        2. Calcul des descripteurs
        3. Normalisation
        4. Prédiction
        
        Args:
            smiles_list (list): Liste de chaînes SMILES
            
        Returns:
            pandas.DataFrame: DataFrame contenant :
                - Original_SMILES: SMILES d'origine
                - Canonical_SMILES: SMILES canonique
                - Predicted_RT: Temps de rétention prédit
        
        Example:
            >>> predictor = RTPredictor()
            >>> results = predictor.predict(['CCO', 'CCC'])
            >>> print(results['Predicted_RT'])
        """
        descriptors_df, canon_smiles = self.compute_descriptors(smiles_list)

        # Normalisation
        X_normalized = descriptors_df.copy()
        for col in X_normalized.columns:
            if X_normalized[col].std() > 0:
                mean = X_normalized[col].mean()
                std = X_normalized[col].std()
                X_normalized[col] = (X_normalized[col] - mean) / std
            else:
                X_normalized[col] = 0

        # Prédiction
        predictions = self.model.predict(X_normalized)

        return pd.DataFrame({
            'Original_SMILES': smiles_list,
            'Canonical_SMILES': canon_smiles,
            'Predicted_RT': predictions
        })