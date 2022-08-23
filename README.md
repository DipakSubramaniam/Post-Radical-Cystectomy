# Post-Radical-Cystectomy

Pranjal Agrawal
Gaurish Agrawal
Dipak Subramaniam

Data analysis of large patient/encounter records into usable format for statistical investigation

The purpose of these steps is to create a list of codes (each code is a procedure) that happened prior to the encounter, or in the same encounter, as in which the patient received radical cystectomy (removal of the bladder). 


Step 1 = drop encounter_id if no end date to encounter (encounter file). Delete these from the encounter file. 
Step 2 = Find date of radical cystectomy (RC) per patient_id (procedure file)
Date of RC is the date at which 
Code = 0TRB07Z OR 0TRB47Z
OR
Code = (0TTB0ZZ, OR 0TRB07Z, OR 0TTB3ZZ, OR 0TTB4ZZ, OR 0TTB8ZZ, OR 0TRB47Z, OR 8E0W4CZ) AND (0T16078, OR 0T16079, OR 0T160J8, 0T160J9, OR  0T160K8, 0T160K9, 0T160Z8, 0T160Z9, 0T17078, 0T17079, 0T170J8, 0T170J9, 0T170K8, 0T170K9, 0T170Z8, 0T170Z9, 0T18078, 0T18079, 0T180J8, 0T180J9, 0T180K8, 0T180K9, 0T180Z8, 0T180Z9, 0T1607A, 0T1607C, 0T160JA, 0T160JC, 0T160KA, 0T160KC, 0T160ZA, 0T160ZC, 0T1707A, 0T1707C, 0T170JA, 0T170JC, 0T170KA, 0T170KC, 0T170ZA, 0T170ZC, 0T1807A, 0T1807C, 0T180JA, 0T180JC, 0T180KA, 0T180KC, 0T180ZA, 0T180ZC, 0T1B07C, 0T1607D, 0T160JD, 0T160KD, 0T160ZD, 0T1707D, 0T170JD, 0T170KD, 0T170ZD, 0T1807D, 0T180JD, 0T180KD, 0T180ZD, 0TRB07Z, 0T16478, 0T16479, 0T164J8, 0T164J9, 0T164K8, 0T164K9, 0T164Z8, 0T164Z9, 0T17478, 0T17479, 0T174J8, 0T174J9, 0T174K8, 0T174K9, 0T174Z8, 0T174Z9, 0T18478, 0T18479, 0T184J8, 0T184J9, 0T184K8, 0T184K9, 0T184Z8, 0T184Z9, 0T1647A, 0T1647C, 0T164JA, 0T164JC, 0T164KA, 0T164KC, 0T164ZA, 0T164ZC, 0T1747A, 0T1747C, 0T174JA, 0T174JC, 0T174KA, 0T174KC, 0T174ZA, 0T174ZC, 0T1847A, 0T1847C, 0T184JA, 0T184JC, 0T184KA, 0T184KC, 0T184ZA, 0T184ZC, 0T1B47C, 0T1647D, 0T164JD, 0T164KD, 0T164ZD, 0T1747D, 0T174JD,0T174KD, 0T174ZD, 0T1847D, 0T184JD, 0T184KD, 0T184ZD, OR 0TRB47Z)
Basically each comma in this long list  is an “OR”
OR
Code = 57.71 OR 57.7 OR 57.79
OR
Code = 51580 OR 51590 OR 51575 OR 51596 OR 51585 OR 51595

Step 3 = Find associated encounter_id for that date per patient_id (procedure file), call this encounter_RC.
Step 4 = Delete all encounter_ids per patient_id that have a start date after the end date of encounter_RC. This step is 	to be done in the procedure file. 
The start and end dates of all encounters are in the encounter file

At this point, the procedure file will only contain encounter_ids that happened prior to, or the encounter in which, the patient received RC. 


The following steps are to designate our three cohorts:
Cohort 1 - Patient has no bladder, has ovaries, has fallopian tubes
Cohort 2 - Patient has no bladder, has ovaries, has no fallopian tubes
Cohort 3 - Patient has no bladder, has no ovaries, and has no fallopian tubes

Step 5 = Copy the procedure file and rename it “Cohort file”

Step 6 = Delete all patients with “maybe” ovaries and/or “maybe” fallopian tubes prior to Encounter_RC (Before deleting, text me the total  number of these patients) based off the following:

Code = 58150 OR 58152 OR 58180 OR 58200 OR 58210 OR 58240

Step 7 = Tell Pranjal the number of unique entries in variable code_systems
There will be a maximum of 4 unique entries = ICD-9, ICD-10, CPT, and HCP..Text me which of the 4 are there

Step 8 = Create variable “cohort.” Any patients with a ICD-9 or HCP.. code system get the value of “4” for the variable cohort. 

Step 9 = All other patients get the value of “1” for the variable cohort. 
The below steps all exclude cohort 4.
If patiens have Code =  58700 OR 58600 OR 58605 OR 58611 OR 58615 OR 58670 OR 58679 OR 0UT7 OR (0UT6 AND 0UT5) → value of “2” for the variable cohort. However, if patients have any of the following code:
Code = 58262 OR  58291 OR 58571 OR 58573 OR  58552 OR 58554 OR 58542 OR 58544 OR 58720 OR  58661 OR 58940 OR 58943 OR (0UT0 AND 0UT1) OR 0UT2 → value of “3” for the variable cohort. 
Things to watch out for:
If patient does not have any of the above listed codes in step 9, they are in cohort 1
If patient has cohort 2 codes only they are in cohort 2
if patient has cohort 2 and cohort 3 codes they are in cohort 3
If patient has cohort 3 codes only they are in cohort 3
Cohort 2 will be a VERY VERY small group of patients
