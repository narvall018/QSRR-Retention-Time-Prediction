# Conclusion

This study developed a QSRR model for retention time prediction in reversed-phase liquid chromatography. The analysis of different molecular descriptor types revealed that combining multiple descriptors enhances prediction accuracy. The LASSO feature selection method, coupled with Ridge regression, achieved an R² of 0.925 on the test set with an RMSE of 0.83 minutes. The model identified physicochemical properties that govern chromatographic retention. Lipophilicity (MolLogP) emerged as the primary factor controlling analyte-stationary phase interactions. Surface area descriptors (PEOE_VSA) and hydrogen bonding capacity (NHOHCount) quantify additional retention mechanisms through molecular polarity. The integration of multiple descriptor types captured complementary aspects of molecular structure affecting retention behavior. While individual descriptor sets (RDKit 2D, ECFP, MACCS) provided partial characterization, their combination improved prediction performance. This indicates that retention time depends on multiple molecular features beyond standard physicochemical properties {cite}`kaliszan2007qsrr`. The analysis demonstrated consistent model performance across the retention time range, with prediction errors typically within ±2 minutes. This level of accuracy supports compound identification by providing retention time windows for suspect screening {cite}`schymanski2014identifying`. The model constitutes a quantitative tool for retention time prediction when reference standards are not available.


# Data availability

The Jupyter notebooks used to generate the results in this study are available at the following link: [https://github.com/narvall018/QSRR-Retention-Time-Prediction/tree/main/notebook](https://github.com/narvall018/QSRR-Retention-Time-Prediction/tree/main/notebook).

The data utilized for model training and testing can be found in the `data` folder, and the final model is saved in the `models` folder within the same repository.


