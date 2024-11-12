from setuptools import setup, find_packages

setup(
    name="QSRR_predictor",
    version="0.1",
    packages=find_packages(),
    package_data={
        'QSRR_predictor.models': ['*.joblib'],
    },
    install_requires=[
        'numpy',
        'pandas',
        'scikit-learn',
        'rdkit',
        'joblib',
        'tqdm',
    ],
)