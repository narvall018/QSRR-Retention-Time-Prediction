# QSRR_predictor/predictor.py

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, rdMolDescriptors, Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from tqdm import tqdm
import joblib
import os

class RTPredictor:
    def __init__(self):
        """Initialisation du prédicteur"""
        # Chargement du modèle et des features
        self.model = joblib.load('models/ridge_all_descriptors_LASSO.joblib')
        self.selected_features = [feat.replace('_x', '') for feat in joblib.load('models/features_all_descriptors_LASSO.joblib')]
        
        # Catégoriser les descripteurs nécessaires
        self.needed_bits = [int(f.split('_')[1]) for f in self.selected_features if f.startswith('Bit_')]
        self.needed_maccs = [int(f.split('_')[1]) for f in self.selected_features if f.startswith('MACCS_')]
        self.needed_mqn = [int(f.split('_')[1]) for f in self.selected_features if f.startswith('MQN_')]
        self.needed_rdkit = [f for f in self.selected_features if not any(f.startswith(p) for p in ['Bit_', 'MACCS_', 'MQN_'])]
        
        print("Nombre de features sélectionnées:", len(self.selected_features))
        print("Répartition des descripteurs:")
        print(f"- ECFP: {len(self.needed_bits)}")
        print(f"- MACCS: {len(self.needed_maccs)}")
        print(f"- MQN: {len(self.needed_mqn)}")
        print(f"- RDKit: {len(self.needed_rdkit)}")

    def canonicalize_smiles(self, smiles):
        """Convertit un SMILES en sa forme canonique"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                return Chem.MolToSmiles(mol, canonical=True)
            else:
                return None
        except:
            return None

    def prepare_molecules(self, smiles_list):
        """Préparation des molécules à partir des SMILES"""
        mols = []
        valid_smiles = []
        invalid_smiles = []
        
        print("\nConversion des SMILES en molécules...")
        for i, smiles in enumerate(tqdm(smiles_list)):
            canon_smiles = self.canonicalize_smiles(smiles)
            if canon_smiles is not None:
                mol = Chem.MolFromSmiles(canon_smiles)
                if mol is not None:
                    mols.append(mol)
                    valid_smiles.append(canon_smiles)
                else:
                    invalid_smiles.append((i, smiles))
            else:
                invalid_smiles.append((i, smiles))
        
        if invalid_smiles:
            print("\nSMILES invalides détectés:")
            for idx, smiles in invalid_smiles:
                print(f"Index {idx}: {smiles}")
        
        return pd.DataFrame({
            'original_smiles': smiles_list,
            'canonical_smiles': valid_smiles,
            'mol': mols
        })

    def compute_ecfp(self, df, column='mol'):
        """Calcule uniquement les bits ECFP nécessaires"""
        if not self.needed_bits:
            return pd.DataFrame()
            
        print("Calcul des ECFP6 sélectionnés...")
        fps = []
        for mol in tqdm(df[column]):
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=1024)
            fp_bits = list(fp.GetOnBits())
            fps.append([1 if i in fp_bits else 0 for i in self.needed_bits])
        
        return pd.DataFrame(fps, columns=[f'Bit_{i}' for i in self.needed_bits])

    def compute_maccs(self, df, column='mol'):
        """Calcule uniquement les clés MACCS nécessaires"""
        if not self.needed_maccs:
            return pd.DataFrame()
            
        print("Calcul des MACCS sélectionnés...")
        maccs_list = []
        for mol in tqdm(df[column]):
            maccs = MACCSkeys.GenMACCSKeys(mol)
            maccs_list.append([maccs[i] for i in self.needed_maccs])
        
        return pd.DataFrame(maccs_list, columns=[f'MACCS_{i}' for i in self.needed_maccs])

    def compute_mqn(self, df, column='mol'):
        """Calcule uniquement les MQN nécessaires"""
        if not self.needed_mqn:
            return pd.DataFrame()
            
        print("Calcul des MQN sélectionnés...")
        mqn_list = []
        for mol in tqdm(df[column]):
            mqn = rdMolDescriptors.MQNs_(mol)
            mqn_list.append([mqn[i-1] for i in self.needed_mqn])
        
        return pd.DataFrame(mqn_list, columns=[f'MQN_{i}' for i in self.needed_mqn])

    def compute_rdkit_descriptors(self, df, column='mol'):
        """Calcule uniquement les descripteurs RDKit nécessaires"""
        if not self.needed_rdkit:
            return pd.DataFrame()
            
        print("Calcul des descripteurs RDKit sélectionnés...")
        calc = MoleculeDescriptors.MolecularDescriptorCalculator(self.needed_rdkit)
        
        desc_list = []
        for mol in tqdm(df[column]):
            desc = calc.CalcDescriptors(mol)
            desc_list.append(desc)
        
        return pd.DataFrame(desc_list, columns=self.needed_rdkit)

    def compute_descriptors(self, smiles_list):
        """Calcule uniquement les descripteurs nécessaires"""
        # Préparation des molécules
        mol_df = self.prepare_molecules(smiles_list)
        
        # Calcul des descripteurs par type
        desc_parts = []
        
        if self.needed_bits:
            desc_parts.append(self.compute_ecfp(mol_df))
        if self.needed_maccs:
            desc_parts.append(self.compute_maccs(mol_df))
        if self.needed_mqn:
            desc_parts.append(self.compute_mqn(mol_df))
        if self.needed_rdkit:
            desc_parts.append(self.compute_rdkit_descriptors(mol_df))
        
        # Combinaison des descripteurs
        descriptors_df = pd.concat(desc_parts, axis=1)
        
        # S'assurer que les colonnes sont dans le bon ordre
        descriptors_df = descriptors_df[self.selected_features]
        
        print(f"\nNombre de descripteurs calculés: {descriptors_df.shape[1]}")
        return descriptors_df, mol_df['canonical_smiles']

    def normalize_features(self, X):
        """Normalisation des features"""
        X_normalized = X.copy()
        for col in X_normalized.columns:
            if X_normalized[col].std() > 0:
                mean = X_normalized[col].mean()
                std = X_normalized[col].std()
                X_normalized[col] = (X_normalized[col] - mean) / std
            else:
                X_normalized[col] = 0
        return X_normalized

    def predict(self, smiles_list):
        """Pipeline complète de prédiction"""
        # Calcul des descripteurs
        descriptors_df, canon_smiles = self.compute_descriptors(smiles_list)

        # Normalisation
        descriptors_normalized = self.normalize_features(descriptors_df)

        # Renommer les colonnes pour correspondre exactement au modèle
        original_features = joblib.load('models/features_all_descriptors_LASSO.joblib')
        descriptors_normalized.columns = original_features

        # Prédiction
        predictions = self.model.predict(descriptors_normalized)

        # Création du DataFrame final
        results_df = pd.DataFrame({
            'Original_SMILES': smiles_list,
            'Canonical_SMILES': canon_smiles,
            'Predicted_RT': predictions
        })

        return results_df