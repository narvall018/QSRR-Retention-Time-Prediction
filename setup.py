from setuptools import setup, find_packages

setup(
    name="qsrr_predictor",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "rdkit",
        "scikit-learn",
        "tqdm",
        "joblib"
    ],
    include_package_data=True,
    package_data={
        "qsrr_predictor": ["models/*.joblib"],
    },
    author="Narvall018",
    description="QSRR model for retention time prediction in liquid chromatography",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/narvall018/QSRR-Retention-Time-Prediction",
)