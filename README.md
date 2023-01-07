# Post-Radical-Cystectomy
Data analysis of large patient/encounter records into usable format for statistical investigation

The baseline patient/encounter data was to be processed in the Stata tool, but the Excel file was far too large. Instead, the data in ‘csv’ format is analyzed using the Pandas library in Python. Such csv data is converted into a DataFrame object, which is two-dimensional, size-mutable and in a tabular format that is easy to modify. All csv data is provided initially as a set of de-identified patient records.

The script reads in the input procedure and encounter csv files. The following steps detail the exact procedures by which the data was manipulated and formulated:

Step 1: Within the encounter file, all encounters without any end date were dropped. Any rows with empty values for the end_date field are excluded. The final list of these excluded encounters is also reflected onto the procedure file, whereas the same encounters are dropped.

Step 2: Within the procedure file, per each patient (patient_id), the date for the radical cystectomy is found. All patients (and all their corresponding encounters) are paired, and the script searches the code field for any matches to a given set of ICD-10, ICD-9 and CPT codes. 

Key (& = and; OR / , = or)

(ICD-10)

{0TRB07Z, 0TRB47Z}

OR

(ICD-10)

{ (0TTB0ZZ, 0TRB07Z, 0TTB3ZZ, 0TTB4ZZ, 0TTB8ZZ, 0TRB47Z, 8E0W4CZ) &
(0T16078, 0T16079, 0T160J8, 0T160J9,  0T160K8, 0T160K9, 0T160Z8, 0T160Z9, 0T17078, 0T17079, 0T170J8, 0T170J9, 0T170K8, 0T170K9, 0T170Z8, 0T170Z9, 0T18078, 0T18079, 0T180J8, 0T180J9, 0T180K8, 0T180K9, 0T180Z8, 0T180Z9, 0T1607A, 0T1607C, 0T160JA, 0T160JC, 0T160KA, 0T160KC, 0T160ZA, 0T160ZC, 0T1707A, 0T1707C, 0T170JA, 0T170JC, 0T170KA, 0T170KC, 0T170ZA, 0T170ZC, 0T1807A, 0T1807C, 0T180JA, 0T180JC, 0T180KA, 0T180KC, 0T180ZA, 0T180ZC, 0T1B07C, 0T1607D, 0T160JD, 0T160KD, 0T160ZD, 0T1707D, 0T170JD, 0T170KD, 0T170ZD, 0T1807D, 0T180JD, 0T180KD, 0T180ZD, 0TRB07Z, 0T16478, 0T16479, 0T164J8, 0T164J9, 0T164K8, 0T164K9, 0T164Z8, 0T164Z9, 0T17478, 0T17479, 0T174J8, 0T174J9, 0T174K8, 0T174K9, 0T174Z8, 0T174Z9, 0T18478, 0T18479, 0T184J8, 0T184J9, 0T184K8, 0T184K9, 0T184Z8, 0T184Z9, 0T1647A, 0T1647C, 0T164JA, 0T164JC, 0T164KA, 0T164KC, 0T164ZA, 0T164ZC, 0T1747A, 0T1747C, 0T174JA, 0T174JC, 0T174KA, 0T174KC, 0T174ZA, 0T174ZC, 0T1847A, 0T1847C, 0T184JA, 0T184JC, 0T184KA, 0T184KC, 0T184ZA, 0T184ZC, 0T1B47C, 0T1647D, 0T164JD, 0T164KD, 0T164ZD, 0T1747D, 0T174JD,0T174KD, 0T174ZD, 0T1847D, 0T184JD, 0T184KD, 0T184ZD, 0TRB47Z) }

OR

(ICD-9)

{57.71, 57.7, 57.79}

OR

(CPT)

{51580, 51590, 51575, 51596, 51585, 51595}

Step 3 - The script finds the date of radical cystectomy per each patient. If any code for a specific patient-encounter matches the codes and conforms to the logic presented above in Step 2, that encounter date is flagged as the date of the radical cystectomy. The flag and date are added as new columns to the main patient DataFrame, known as pat_dfp. Otherwise, the patient is marked as not having a radical cystectomy date. After reviewing all the patient-encounter pairs, any patient whose codes (through all their encounters) did not register a radical cystectomy is removed from pat_dfp. 

Step 4 - All encounter_id values per patient_id values that have a start_date after the end_date of radical cystectomy are deleted. Firstly, the ending dates for all the radical cystectomies are recorded after examining the encounter file. The script parses for encounters which took place after that patient’s radical cystectomy end date, and subsequently deletes these encounters.

Note: At this juncture, the procedure file will only contain encounter_ids that happened prior to, or the encounter in which, the patient received their radical cystectomy.

The following steps are to designate our three cohorts:
Cohort 1 - Patient has no bladder, has ovaries, has fallopian tubes
Cohort 2 - Patient has no bladder, has ovaries, has no fallopian tubes
Cohort 3 - Patient has no bladder, has no ovaries, and has no fallopian tubes

Step 5 - The procedure file is copied as a new cohort file for use in the upcoming steps. 

Step 6 - All patients whose encounters fall under ‘maybe ovaries/fallopian_tubes’ categories before their radical cystectomy encounter are deleted. This is determined by the set of CPT codes: {58150, 58152, 58180, 58200, 58210, 58240}. Identifying whether the encounters fall under these categories requires searching all patient-encounter pairs and marking matches to the given list of codes. Any patients who were marked as such are deleted from the pat_dfp DataFrame.

Step 7 - Each patient's radical cystectomy encounter was determined by codes in one of four code categories (ICD-10, ICD-9, CPT, HCP). The frequencies of these categories across all the radical cystectomy patients are tabulated and displayed.

Step 8 - A cohort variable is added to the pat_dfp DataFrame. Automatically, any patients who have the ICD-9 or HCP code systems are given a cohort designation of 0 (an exclusionary marking, since the focus of the analysis is on the ICD-10 and CPT codes.

Step 9 - Any patients without the cohort 0 designation are automatically assigned a cohort value of 1, unless:

If a patient has an current cohort value of 1 and code field = 58700 OR 58600 OR 58605 OR 58611 OR 58615 OR
58670 OR 58679 OR 0UT7* OR (0UT6* AND 0UT5*) OR 0U57* OR (0U56* AND 0U55*), then they are assigned a cohort value of 2, unless:

If a patient has an current cohort value of 2 and code field = 58262 OR 58291 OR 58571 OR 58573 OR 58552 OR 58554 OR 58542 OR 58544 OR 58720 OR 58661 OR 58940 OR 58943 OR (0UT0* AND 0UT1*) OR 0UT2* OR (0U50* AND 0U51*) OR 0U52*, then they are assigned a cohort value of 3. 

Note: * implies that the full code starts with the given set of characters. 

Step 9.5 - Data structures and shared data headers for Steps 10-12 are prepared. All patients with radical cystectomies are double checked and made sure they are only counted once. The script also reads in the patient, lab_result, diagnosis, vital_signs and medication_drug .csv files, which will be mined for information and data needed in the next three steps.

Step 10 - Determining patient information at the Time of Cystectomy or Up till the time of Cystectomy

Radical Cystectomy Date (Date)
Year of Birth
Race (White, Black, Asian, Other)
Ethnicity (Hispanic, Not Hispanic)
BMI – in lab_result and vital_signs file
Year of Death
(Z80.3*) Family History of Breast Cancer (Yes/No)
(Z85.3*) Personal History of Breast Cancer (Yes/No)
(Z80.41*) Family History of Ovarian Cancer (Yes/No)
(Z85.43*) Personal History of Ovarian Cancer (Yes/No)
(Z85.038*) Personal History of Colorectal Cancer (Yes/No)
(V84.09*) Lynch Syndrome (Yes/No)
(Z15.01*) Hereditary Breast and Ovarian, Li-Fraumeni (Yes/No)
(Z72.0*) Tobacco Use (Yes/No)
(I10*-I16*) Hypertensive Diseases (Yes/No)
(E08*-13*) Diabetes Mellitus (Yes/No)
(E78*) Hyperlipidemia (Yes/No)
(I20*-25*) ischemic heart disease (Yes/No)
(E66*) Overweight/Obesity (Yes/No)
(E28*) Polycystic Ovarian Disease/Syndrome (Yes/No)
(N80*) History of Endometriosis (Yes/No)
(N73*) History of PID (Yes/No)
(HS200*) History of Systemic Contraceptives (Yes/No) – in medication_drug file
(HS300*) History of estrogens (Yes/No) – in medication_drug file
(HS800*) History of Progestins (Yes/No) – in medication_drug file 
(1013911 OR 1008870) History of hysterectomy (Yes/No) – in procedure file, a CPT code. 

Step 11 - Determining patient outcomes at any time in their lifetimes (before, during, or after RC) Outcomes - all Yes/No
	If Yes -> List Year of first-time diagnosis made
(F03*) Dementia
(G31.84* or F09*) Cognitive Impairment
(G20*) Parkinson’s
(M80* or M81*) Osteoporosis
(I20*-25*) ischemic heart disease
(I21* or 22*) Myocardial Infarction
(F30*-39*) Affective Mood Disorders
(F41*) Anxiety
(F32*, F33*) Depression
(F52*) Female Sexual Dysfunction 
(I60*-63*) Stroke
(I82.4* or 82.6*) Deep Vein Thrombosis
(I26*) Pulmonary Embolism
(N81*) Genital Prolapse

Step 12 - Determining patient outcomes during or after RC Encounter – all Yes/No
If Yes -> List Year of first-time diagnosis made 
(F03*) Dementia
(G31.84* or F09*) Cognitive Impairment
(G20*) Parkinson’s
(M80* - 81*) Osteoporosis
(I20*-25*) ischemic heart disease
(I21* or I22*) Myocardial Infarction
(F30*-39*) Affective Mood Disorders
(F41*) Anxiety
(F32* F33*) Depression
(F52*) Female Sexual Dysfunction
(I60*- 63*) Stroke
(I82.4* or I82.6*) Deep Vein Thrombosis
(I26*) Pulmonary Embolism
(N81*) Genital Prolapse
(C56*) Ovarian Cancer
(C48.1* or 48.2*) Peritoneal Cancer
(C50*) Breast Cancer

Step 13 - After all the data mining, the requested variable data can be conglomerated into a single excel file.

Step 14 - The script uses the Matplotlib library to create relevant data plots post-analysis. 



