# AlloEF: An Ensemble Model for Protein Allosteric Site Identification Based on Transfer Entropy and Energetic Frustration 

Here, we develop AlloEF, an effective method for protein allosteric site prediction, which adopts a soft-voting classifier with LightGBM, Random Forest and XGBoost combined.

Authors: Jilong Zhang, Xiaohan Sun, Zhixiang Wu, Jingjie Su, Xinyu Zhang, Chunhua Li

The process includes three steps: feature extraction, prediction integration and prediction.

The following is an example of an allosteric protein (PDB ID: 1WQW) to introduce the implementation process of the AlloEF method to predict protein allosteric sites.  The forecasting process is done on Linux. MATLAB software is required in the process of extracting feature features.

## Step 1 feature extraction

* Python version: 3.10
numpy>=1.21.0
pandas>=1.3.0
scikit-learn>=1.0.0
xgboost>=1.5.0
lightgbm>=3.2.0
shap>=0.40.0
imbalanced-learn>=0.8.0
boruta>=0.3.0
tensorflow>=2.6.0
matplotlib>=3.5.0
openpyxl>=3.0.0
pyyaml>=6.0
Matlab R2023a

1、Here, we write the corresponding PDB file of protein 1WQW as 1WQW.pdb, and then extract the sequence of the corresponding chain (A chain) and save it as 1WQW_A. pdb. After that, run the "protein_information.py" program to extract protein information and obtain 1WQW_A_protein_info.csv file for easy calculation.

2、Run "physicochemical_features.py" (command line: python ./physico_feature.py) Note that the python program file is in the same folder as the downloaded txt file of the Aaindex database (which can be downloaded directly using the download_aaindex.py program). "aaindex1.txt" and "1WQW_A. pdb" and "1WQW_A_protein_info.csv" as input files, then get "1WQW_Aaa.csv".

3、Access the server or build an environment in the Linux system, directly run "pssm.py" to calculate the pssm matrix and extract the data to generate it into the "1WQW_A.csv", the database used at this time is uniref50, E-value = 0.001, -num_iterations = 3 Other options use the default parameters.

4、Secondary structure
a. Download and install the DSSP tool in "https://swift.cmbi.umcn.nl/gv/dssp"; Or download and compile it directly in the Linux system.
b. Run the "dssp.py" script file (command line: python dssp.py) and finally get the output file as "1WQW_A.csv".

5、 CX/DPX 
Download and install the psaia.exe program on the https://sourceforge.net/projects/psaia/ website.
a. Run the "psaia.exe";
b. Step by step through the "Structure Analyser" tab control, then find "Analysis Types" and check both "Analyse as Bound" and "Analyse by Chain". All parameters are set to default;
c. Enter the pdb file "1WQW_A.pdb" into the program and click "Run" to get the result, i.e. the files "1WQW_A _unbound.tbl" and "1WQW_A _bound.tbl".

6. CavityPlus 2022
a. Log in to the CavityPlus website（http://pkumdl.cn:8000/cavityplus/） and upload 1WQW_A.pdb file, which will be calculated and output
b. Run the "get_ispocket.py" script command for the output folder to get the 1WQW_A_ispocket.csv at the end

7.Network-based topological features and protein dynamics features
The topological features of the network, including centrality and mediation, are calculated by complex networks, and the dynamic features include residue transfer entropy, and the residue fluctuation and residue sensitivity are calculated by dGNM, GNM, and ANM, respectively. The above are calculated by MATLAB software, you need to put "MAIN_gnmte_zjl.m", "pdbread.m", "GNM_sel_zjl.m", "comnetbet.m", "anmselect_zjl.m", "ANM_calPRS_zjl.m", etc. in the same directory, and run "MAIN_gnmte_zjl.m" directly. All of the above features can be calculated directly and a 1WQW_A_dynamics.csv can be obtained

8.Protein energetic frustration feature extraction
Log in to the Frustratometer Server website(http://frustratometer.qb.fcen.uba.ar/), upload the 1WQW_A.pdb file, and output the result after calculation. For the output folder, find the file at the end of the pdb_singleresidue, extract the information of the last column of "FrstIndex" and integrate it into the 1WQW_A_ frustration.csv.

## Step 2 prediction integration
After all the feature extraction is completed, integrate all the feature files mentioned above and put them together in the initial 1WQW_A_protein_info.csv to generate 1WQW_Aall.csv file (note that the feature naming is consistent with the content of the main text). Run SNFA(Spatial Neighborhood Feature Aggregation).py script to get the final input file 1WQW_A_allfeature.csv.

## Step 3 prediction
Place 1WQW_A_allfeature.csv in the Input_data folder and place it in the same path as the main.py, preprocessing.py, modeling.py, and evaluation.py files.
Command line: "python ./ main.py ./Input_data/1WQW_A _feature.csv", the output is ". results/predictions.xlsx", where the predicted site is stored,


Help
For any questions, please contact us at chunhuali@bjut.edu.cn.
