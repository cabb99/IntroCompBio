Optimization of Acute Myeloid Leukemia Predictions with a Five-Fold Cross-Validated Genetic Algorithm
======================================================================================================
##### Carlos Bueno<sup>1</sup>, Luiza Ferreira<sup>2</sup>, John Gawedzinski<sup>3</sup>, Sangheon Han<sup>3</sup>, Sohyun Park<sup>3</sup>, Trenton Piepergerdes<sup>3</sup>, Amina A. Qutub<sup>1,3</sup>
##### Department of <sup>1</sup>Systems, Synthetic and Physical Biology, <sup>2</sup>Chemistry and <sup>3</sup>Bioengineering, Rice University

Introduction: Acute Myeloid Leukemia (AML) is a rapidly progressing cancer of myeloid cells within the bone marrow and blood. The crowd-sourced DREAM 9 Challenge seeks to improve prognosis of AML based on a reverse phase protein array (RPPA) dataset as well as clinical variables (e.g., age at diagnosis, gender) provided by the M.D. Anderson Cancer Center. RPPA data contains expression levels of 271 different proteins and phosphoproteins for 191 AML patients, and a 74-patient test set. Using this proteomic data with a combination of clinical variables, we developed a unique genetic algorithm approach to predict four patient outcomes: complete response (pCR), relapse (pRelapse), remission duration (Remission), and overall survival (OS). We then determined proteins pertinent to the performance of each outcome prediction. Overall, this work provides a new method to optimize models of clinical predictions and offers insight into key factors that contribute to leukemia progression. 

Materials and Methods: Categorical variables (e.g., gender, prior infection) in the training set were converted to binary variables. The training set was divided into 5 groups twice. A Five-Fold Cross-Validated Genetic Algorithm was developed to search the variable space and to find the optimal combination of methods for an accurate prediction. Balanced accuracy was used as a prediction function for the genetic algorithm. Accuracy of the genetic algorithm was assessed by using three different prediction methods: linear regression, random forest classifier and random forest regressor. We submitted an average of 25 repetitions and the best repetition of each method. The final scores are measured as the Euclidean distance to the values from a withheld 74-patient test set.

Results and Discussion: Random forest classifier provided the best predictor for patients’ Remission and OS duration (Table 1). The proteins PRKCD, ZNF346, STK11 and ARC were repeatedly found as predictive of Remission. Proteins NOTCH1, EIF2S1 and PRKCB were highly implicated for OS. pCR and pRelapse were scored the best by the random forest regressor. pCR included NPM1 and CDK1 as pertinent genes. Genes RPS6KB1, ERBB2, SMAD5, BAK1, DLX1, KDR, GSKA_B, MDM4, MAPT and TP53 were the important parameters to predict pRelapse. 

|Prediction Method|pCR Score|pRelapse Score|Remission Score|OS Score|Overall Score|
|-----------------|---------|--------------|---------------|--------|-------------|
|Best Linear Regression|	5.71|	4.42|	3.87|	5.15|	9.68|
|Average Linear Regression|	4.62|	2.99|	4.47|	5.34|	8.87|
|Best Random Forest Classifier|	5.83|	3.16|	3.54|	3.71|	8.38|
|Average Random Forest Classifier|	4.58|	3.03|	3.87|	3.91|	7.77|
|Best Random Forest Regressor|	3.80|	3.02|	4.50|	5.00|	8.29|
|Average Random Forest Regressor|	3.99|	2.88|	4.24|	5.74|	8.67|

_Table 1. Scores obtained with the genetic algorithm using different prediction methods. pCR: complete response probability; pRelapse: relapse probability; Remission: remission duration. OS: overall survival duration. Best predictors are highlighted._

Conclusions: The Five-Fold Cross-Validated Genetic Algorithm we developed allow us to screen pertinent proteins associated with AML, and to predict leukemia patient outcomes. Some proteins that were not associated with AML played a role in the patient outcome predictions. This may shed a light on identifying new biomarkers and candidate target genes for AML prognosis and treatment. More broadly, the algorithm we developed may be utilized as a way to handle categorical and continuous high-dimensional clinical data for predictive models with small training sample sizes; and applied to screen patients by identifying relevant biomarkers for disease prognosis.


Installation Requirements
-------------------------
python2

pandas

numpy

xlrw

scipy

scikit-learn
