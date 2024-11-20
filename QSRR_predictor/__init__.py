import os

# Définir le chemin vers les modèles comme une variable de package
MODELS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models')

from .predictor import RTPredictor