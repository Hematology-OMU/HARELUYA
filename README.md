# Personalized Strategy for Allo-HSCT Guided by Machine Learning
<p>
  To clarify the clinical value of the allogeneic hematopoietic stem cell transplanation (allo-HSCT) strategy guided by machine learning models, we retrospectively investigated the prognostic effects of a personalized transplant procedure recommended by machine learning models using THE JAPANESE TRANSPLANT REGISTRY UNIFIED MANAGEMENT PROGRAM (TRUMP) data. Furthermore, we developed a interactive web application tool "HARELUYA", which plots personalized prediction overall survival (OS) curves after allo-HSCT for ten transplant procedures.</p>

<p>
  The directory of "Analysis" includes the source codes to analyze this study. The main analysis was performed using R. On the other hand, the random search for hyper-parameters  and calculating the predictive probabilities for each patient in the DeepSurv model was done using Python.
The directory of "data" under "Analysis" includes "analysed_df.csv". As an example, only three patients' data is inputted into "analysed_df.csv" here.
</p>

<p>
  "HARELUYA" is a calculator of personalized overall survival (OS) curves for ten allogeneic hematopoietic stem cell transplantation (allo-HSCT) procedures using random survival forest (RSF).
The "hareluya.R" under the directory of "HARELUYA" is the source code of the web application "HARELUYA".
</p>

#####About "analysed_df.csv:
<p>The column of "analysed_df.csv"<br/>
-study_ID: Unique patient identifier<br/>
-.SCT.Year: year at allo-HSCT<br/>
integer<br/>

<li>Prognostic predictors;</li>
-.Age: Recipent's age at allo-HSCT<br/>
integer<br/>
-.Sex.Mismatch3: Relationship between donor's and recipent's sex<br/>
factor;<br/>
1: Female donor to male recipent<br/>
0: Other than 1<br/>
-.PS24: Performance status at allo-HCT<br/>
integer; range: 0-4<br/>
-dx_to_sct_day_trump: Number of days from diagnosis to allo-HSCT<br/>
integer<br/>
-.Analysed_Disease: Underlying hematological disease<br/>
factor; "ALL", "AML", "MDS", "ML", "ATL" or "MPN"<br/>
-.Disease.satatus.all: Disease status at allo-HSCT<br/>
factor; "CR" or "nonCR"<br/>
-R_CMVAntiB: Recipent's cytomegalovirus serostatus<br/>  
factor; "yes" or "no"<br/>
- .HCT.CI: Hematopoietic cell transplantation-specific comorbidity index (HCT-CI) at allo-HSCT<br/>
factor;<br/>
0: 0<br/>
1: 1<br/>
2: 2<br/>
3: Equal to or higher than 3<br/>
- Tx_pattern: Transplatation procedure<br/>
factor;<br/>
1: RIC/MRD/non-PTCY [R-MRD]<br/>
2: RIC/MUD/non-PTCY [R-MUD]<br/>
3: RIC/UCB/non-PTCY [R-UCB]<br/>
4: RIC/Haplo/PTCY [R-Haplo-CY]<br/>
5: RIC/Haplo/non-PTCY [R-Haplo-nCY]<br/>
6: MAC/MRD/non-PTCY [M-MRD]<br/>
7: MAC/MUD/non-PTCY [M-MUD]<br/>
8: MAC/UCB/non-PTCY [M-UCB]<br/>
9: MAC/Haplo/PTCY [M-Haplo-CY]<br/>
10: MAC/Haplo/non-PTCY [M-Haplo-nCY]<br/>
</p>

<li>Outcome;</li>
-.OS: Observed status<br/>
integer;<br/>
1: death<br/>
0: censoring<br/>
-.MonthOS: Observed follow-up months<br/>
numeric<br/>

# Requirements
This application requires the following to run:
<p>
<li>R version: 4.1.0</li>
<li>Python version: 3.9.6</li>
</p>

<p>
package version<br/>
<li>R package:</li>
shiny: ver. 1.7.1<br/>
randomForestSRC: ver. 2.12.0<br/>
survivalmodels: ver. 0.1.8<br/>
mlr3: ver. 0.11.0<br/>
timeROC: ver. 0.4<br/>

<li>Pyhton library:</li>
pycox: ver. 0.2.2<br/>
torch: ver. 1.9.0<br/>
torchtuples: ver. 0.2.0<br/>
</p>

# License
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.<br/> 
© 2025, Hiroshi Okamura
