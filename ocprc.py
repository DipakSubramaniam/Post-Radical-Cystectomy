# Post Radical Cystectomy Processing
# Dipak Subramaniam

import time
import pandas as pd
from dask import dataframe as dd
import numpy
from collections import Counter
import warnings

# * timing
begin = time.time()

# * warning suppressor
warnings.filterwarnings('ignore')  # suppresses SettingWithCopyWarning in steps 2 & 3

# * reading in procedure file
ipfile = input('* enter the procedure csv file name: ')
pro = dd.read_csv(ipfile)  # procedure data frame
pro.to_parquet('pr.parquet', engine='pyarrow')  # storing in parquet format

# * reading in encounter file
iefile = input('* enter the encounter csv file name: ')
enc = dd.read_csv(iefile)  # encounter data frame
enc.to_parquet('en.parquet', engine='pyarrow')  # storing in parquet format

# * storing in the more efficient Apache Parquet format for easier read/write
dfp = dd.read_parquet('pr.parquet', engine='pyarrow', storage_options={"anon": True, "use_ssl": True},)
dfe = dd.read_parquet('en.parquet', engine='pyarrow', storage_options={"anon": True, "use_ssl": True},)

# * step 1 - dropping encounters if no end date present
print('* step 1 - dropping encounters without an end date in encounter csv file')
dfe_drop = dfe.dropna(subset=['end_date'])  # dropping rows if end nan values (empty spaces) exist for the end date
new_dfe = dfe_drop.compute()
old_dfe = dfe.compute()
rows_dropped = old_dfe.loc[~old_dfe.set_index(list(old_dfe.columns)).index.isin(new_dfe.set_index(list(new_dfe.columns)).index)]
encounters_dropped = rows_dropped['encounter_id'].tolist()
old_dfp = dfp.compute()
dfp_s1 = old_dfp[~old_dfp['encounter_id'].isin(encounters_dropped)]
print('* updated dropped encounters on procedure csv file')
rows_dropped_pro = old_dfp.loc[~old_dfp.set_index(list(old_dfp.columns)).index.isin(dfp_s1.set_index(list(dfp_s1.columns)).index)]
print('* step 1 complete')

# * steps 2 & 3 - finding date and encounter_id of radical cystectomy per patient_id
print('* steps 2 & 3 - finding date of radical cystectomy per patient')
icd_ten_a_codes = ['0TRB072', '0TRB47Z']  # a, b and c lists are ICD-10 codes
icd_ten_b_codes = ['0TTB0ZZ', '0TRB07Z', '0TTB3ZZ', '0TTB4ZZ', '0TTB8ZZ', '0TRB47Z', '8E0W4CZ']
icd_ten_c = '0T16078, 0T16079, 0T160J8, 0T160J9, 0T160K8, 0T160K9, 0T160Z8, 0T160Z9, 0T17078, 0T17079, 0T170J8, 0T170J9, 0T170K8, 0T170K9, 0T170Z8, 0T170Z9, 0T18078, 0T18079, 0T180J8, 0T180J9, 0T180K8, 0T180K9, 0T180Z8, 0T180Z9, 0T1607A, 0T1607C, 0T160JA, 0T160JC, 0T160KA, 0T160KC, 0T160ZA, 0T160ZC, 0T1707A, 0T1707C, 0T170JA, 0T170JC, 0T170KA, 0T170KC, 0T170ZA, 0T170ZC, 0T1807A, 0T1807C, 0T180JA, 0T180JC, 0T180KA, 0T180KC, 0T180ZA, 0T180ZC, 0T1B07C, 0T1607D, 0T160JD, 0T160KD, 0T160ZD, 0T1707D, 0T170JD, 0T170KD, 0T170ZD, 0T1807D, 0T180JD, 0T180KD, 0T180ZD, 0TRB07Z, 0T16478, 0T16479, 0T164J8, 0T164J9, 0T164K8, 0T164K9, 0T164Z8, 0T164Z9, 0T17478, 0T17479, 0T174J8, 0T174J9, 0T174K8, 0T174K9, 0T174Z8, 0T174Z9, 0T18478, 0T18479, 0T184J8, 0T184J9, 0T184K8, 0T184K9, 0T184Z8, 0T184Z9, 0T1647A, 0T1647C, 0T164JA, 0T164JC, 0T164KA, 0T164KC, 0T164ZA, 0T164ZC, 0T1747A, 0T1747C, 0T174JA, 0T174JC, 0T174KA, 0T174KC, 0T174ZA, 0T174ZC, 0T1847A, 0T1847C, 0T184JA, 0T184JC, 0T184KA, 0T184KC, 0T184ZA, 0T184ZC, 0T1B47C, 0T1647D, 0T164JD, 0T164KD, 0T164ZD, 0T1747D, 0T174JD,0T174KD, 0T174ZD, 0T1847D, 0T184JD, 0T184KD, 0T184ZD, 0TRB47Z'
icd_ten_c = icd_ten_c.replace(" ", "")  # removing all whitespaces
icd_ten_c_codes = icd_ten_c.split(',')  # parsing into list based on ','
icd_nine_codes = ['57.71', '57.7', '57.79']  # icd-9 codes for radical cystectomy
cpt_codes = ['51580', '51590', '51575', '51596', '51585', '51595']  # '51597' taken out for now, tbd, and '51570', used to be at the beginning of this list
count = 0
rc_flag = []
unique_patients = []
pat_dfp = dfp_s1  # pat_dfp = dfp.compute() # revert if needed
uniq_pat_list = pat_dfp['patient_id'].unique()
patient_rc_dict = dict.fromkeys(uniq_pat_list, 0)
patient_enc_dict = dict.fromkeys(uniq_pat_list, '')
rc_date_dict = dict.fromkeys(uniq_pat_list, '')
rc_end_date_dict = dict.fromkeys(uniq_pat_list, '')
gbpe = pat_dfp.groupby(['patient_id', 'encounter_id'])
cpe = []  # codes per encounter
print('* searching all the patient-encounter pairs and data from each pair')
for pe, grp in gbpe:  # searching all the patient-encounter pairs and each pair's data containing the other columns
    cpe = grp['code']
    tmp = []
    if set(icd_ten_a_codes).intersection(cpe):
        set_a = set(icd_ten_a_codes).intersection(cpe)
        tmp.extend(set_a)
        cond_a = True
    else:
        cond_a = False
    if set(icd_ten_b_codes).intersection(cpe):
        set_b = set(icd_ten_b_codes).intersection(cpe)
        tmp.extend(set_b)
        cond_b = True
    else:
        cond_b = False
    if set(icd_ten_c_codes).intersection(cpe):
        set_c = set(icd_ten_c_codes).intersection(cpe)
        tmp.extend(set_c)
        cond_c = True
    else:
        cond_c = False
    if set(icd_nine_codes).intersection(cpe):
        set_d = set(icd_nine_codes).intersection(cpe)
        tmp.extend(set_d)
        cond_d = True
    else:
        cond_d = False
    if set(cpt_codes).intersection(cpe):
        set_e = set(cpt_codes).intersection(cpe)
        tmp.extend(set_e)
        cond_e = True
    else:
        cond_e = False
    if ((cond_a or (cond_b and cond_c)) or cond_d or cond_e) and patient_rc_dict[pe[0]] == 0:  # all conditions met
        patient_rc_dict[pe[0]] += 1
        patient_enc_dict[pe[0]] = pe[1]
        rc_end_date_dict[pe[0]] = '0'
        date = grp['date']
        if isinstance(date.iloc[pd.Index(cpe).get_loc(tmp[0])], (int, numpy.integer)):  # if it is just an integer
            rc_date_dict[pe[0]] = date.iloc[pd.Index(cpe).get_loc(tmp[0])]
        else:  # if it is a pandas series
            tmp_list = date.iloc[pd.Index(cpe).get_loc(tmp[0])].tolist()
            rc_date_dict[pe[0]] = next(x for x in tmp_list if len(str(x)) == 8)
        rc_flag.append(1)
        count += 1
    else:
        rc_flag.append(0)
print('* counted all radical cystectomies, stored corresponding encounter_id and date values')
#  * removing x# patients who did not register an radical cystectomy but were supposed to
patient_rc_dict = {k: v for k, v in patient_rc_dict.items() if v != 0}  # patient, rc_flag
patient_enc_dict = {k: v for k, v in patient_enc_dict.items() if v != ''}  # patient, encounter_rc
rc_date_dict = {k: v for k, v in rc_date_dict.items() if v != ''}  # patient, date_rc
rc_end_date_dict = {k: v for k, v in rc_end_date_dict.items() if v != ''}  # patient, end date dummy value
pat_dfp = pat_dfp.reset_index()
flag_rc = []
print('* adding and populating rc column with radical cystectomy flags')
for ind, row in pat_dfp.iterrows():
    cur_patient = row['patient_id']
    cur_encounter = row['encounter_id']
    if cur_patient in patient_enc_dict and cur_encounter == patient_enc_dict[cur_patient]:
        flag_rc.append(1)
    else:
        flag_rc.append(0)
pat_dfp.insert(loc=9, column='rc', value=flag_rc)
print('* removed a small fraction of patients wherein their codes do not register a radical cystectomy')  # , len(patient_rc_dict)
print('* steps 2 & 3 complete')

# * step 4 - deleting all encounters per patient that have a start date after the end date of radical cystectomy
print('* step 4 - deleting all encounters per patient that have a start date after the end date of radical cystectomy')
print('* end dates being tabulated for all radical cystectomies')
indices_to_delete = []
new_dfe = new_dfe.reset_index()
for index, row in new_dfe.iterrows():
    cur_patient = row['patient_id']
    cur_start_date = row['start_date']
    cur_end_date = row['end_date']
    if (cur_patient in rc_date_dict) and (cur_start_date == int(rc_date_dict[cur_patient])):
        rc_end_date_dict[cur_patient] = str(int(cur_end_date))
print('* parsing for encounters that happened after the radical cystectomies')
count = 0
pat_dfp = pat_dfp.reset_index()
for ind, row in pat_dfp.iterrows():
    cur_patient = row['patient_id']
    cur_encounter = row['encounter_id']
    cur_start_date = row['date']
    if cur_patient in rc_date_dict:
        if cur_start_date > int(rc_end_date_dict[cur_patient]):  # current encounter start date vs rc end date
            indices_to_delete.append(1)
            count += 1
        else:
            indices_to_delete.append(0)
    else:
        indices_to_delete.append(0)
print('* identified number of rows: ', count)
pat_dfp.insert(loc=8, column='dele', value=indices_to_delete)
pd.set_option('display.max_columns', None)
pat_dfp = pat_dfp[pat_dfp.dele != 1]
print('* deleted rows containing encounters with dates after the end of the patient radical cystectomy')
print('* step 4 complete')

# * step 5 - preparing cohort file
print('* step 5 - copying procedure file into cohort file')
new_filename = 'cohort.csv'
print('* new cohort file set to be created: ', new_filename)
print('* step 5 complete')

# * step 6 - removing certain patients no under consideration
print('* step 6 - deleting patients who fall under maybe ovaries/fallopian_tubes categories before radical cystectomy')
uniq_cohort_list = pat_dfp['patient_id'].unique()
cohort_rc_dict = dict.fromkeys(uniq_cohort_list, 0)
maybe_codes = ['58150', '58152', '58180', '58200', '58210', '58240']
gbpe2 = pat_dfp.groupby(['patient_id', 'encounter_id'])
maybe_count = 0
cpe = []
for pe, grp in gbpe2:  # searching all the patient-encounter pairs and each pair's data containing the other columns
    cpe = grp['code']
    if set(cpe).intersection(maybe_codes) and cohort_rc_dict[pe[0]] == 0:
        maybe_count += 1
        cohort_rc_dict[pe[0]] += 1
res = Counter(cohort_rc_dict.values())
print("* total patients: ", len(cohort_rc_dict))
print("* patients having maybe codes: ", str(dict(res)[1]))
print("* patients not having such codes: ", str(dict(res)[0]))
maybe_list = [key for key, value in cohort_rc_dict.items() if value == 1]
print('* rows in cohort dataframe before deletion: ', pat_dfp.shape[0])
print('* deleting patients with maybe codes')
maybe_deletes = []
pat_dfp = pat_dfp.reset_index(drop=True)
for inx, row in pat_dfp.iterrows():
    if row['patient_id'] in maybe_list:
        maybe_deletes.append(1)
    else:
        maybe_deletes.append(0)
pat_dfp.insert(loc=8, column='maybe', value=maybe_deletes)
pat_dfp = pat_dfp[pat_dfp.maybe != 1]
print('* rows in cohort dataframe after deletion: ', pat_dfp.shape[0])
print('* final count of patients: ', len(pat_dfp['patient_id'].unique()))
print('* step 6 complete')

# * step 7 - code system frequencies
print('* step 7 - finding frequency of each unique code system')
# frequency = pat_dfp['code_system'].value_counts()
frq_values = pat_dfp['code_system'].value_counts().keys().tolist()
frq_counts = pat_dfp['code_system'].value_counts().tolist()
print("* unique group values in column 'code_system' : ", frq_values)
print("* frequency of values in column 'code_system' : ", frq_counts)
# print(frequency)
print('* step 7 complete')

# * steps 8 & 9 - designating cohorts
print('* steps 8 & 9 - determining patient cohorts')
cohort_uniq_patients = pat_dfp['patient_id'].unique()
cohort_last_dict = dict.fromkeys(cohort_uniq_patients, 0)
coh3cpt = ['58262', '58291', '58571', '58573', '58552', '58554', '58542', '58544', '58720', '58661', '58940', '58943']
coh3icd10a = ['0UT0']
coh3icd10b = ['0UT1']
coh3icd10c = ['0UT2']
coh2cpt = ['58700', '58600', '58605', '58611', '58615', '58670', '58679']
coh2icd10a = ['0UT5']
coh2icd10b = ['0UT6']
coh2icd10c = ['0UT7']
code_sys = ['CPT', 'HCPCS', 'ICD-10-PCS', 'ICD-9-CM']
gbpe3 = pat_dfp.groupby(['patient_id', 'encounter_id'])
cpe = []  # codes per encounter
cspe = []  # code systems per encounter
for pe, grp in gbpe3:  # searching all the patient-encounter pairs and each pair's data containing the other columns
    cpe = grp['code']
    cspe = grp['code_system']
    also_in_three = False  # variable that is toggled if codes are in cohort 3 and 2, prevents cohort 2 from claiming it
    if cohort_last_dict[pe[0]] == 0:  # checking if patient has already been accounted for
        cohort_last_dict[pe[0]] = 4  # everyone gets 4 by default, 4 = non cpt/icd-10 code patients
        # print('code systems for this pair: ', cspe)
        if set(cspe).intersection(['CPT']) or set(cspe).intersection(['ICD-10-PCS']):  # if code system is CPT or ICD-10-PCS, 1 = any CPT/ICD10 code type
            cohort_last_dict[pe[0]] = 1
            # print('codes for pair: ', cpe)
            # if codes also intersect with cohort 3 designated codes --> cohort 3
            if set(cpe).intersection(coh3cpt) or set([cd[:4] for cd in cpe]).intersection(coh3icd10c) or (set([cd[:4] for cd in cpe]).intersection(coh3icd10a) and set([cd[:4] for cd in cpe]).intersection(coh3icd10b)):
                cohort_last_dict[pe[0]] = 3
            # if codes also intersect with cohort 3 and 2 designated codes, --> cohort 3
            if (set(cpe).intersection(coh3cpt) or set([cd[:4] for cd in cpe]).intersection(coh3icd10c) or (set([cd[:4] for cd in cpe]).intersection(coh3icd10a) and set([cd[:4] for cd in cpe]).intersection(coh3icd10b))) and (set([cd[:4] for cd in cpe]).intersection(coh2cpt) or set([cd[:4] for cd in cpe]).intersection(coh2icd10c) or (set([cd[:4] for cd in cpe]).intersection(coh2icd10a) and set([cd[:4] for cd in cpe]).intersection(coh2icd10b))):
                cohort_last_dict[pe[0]] = 3
                also_in_three = True  # will prevent next if statement from being processed
            # if codes also intersect with cohort 2 designated codes --> cohort 2
            if (set(cpe).intersection(coh2cpt) or set([cd[:4] for cd in cpe]).intersection(coh2icd10c) or (set([cd[:4] for cd in cpe]).intersection(coh2icd10a) and set([cd[:4] for cd in cpe]).intersection(coh2icd10b))) and not also_in_three:
                cohort_last_dict[pe[0]] = 2

res = Counter(cohort_last_dict.values())
print("* total patients: ", len(cohort_last_dict))
print('* cohort breakdown: ', res)

# * further code for part 8 and subsequent conversion to updated cohort csv excel file
c4_list = [key for key, value in cohort_last_dict.items() if value == 4]
c3_list = [key for key, value in cohort_last_dict.items() if value == 3]
c2_list = [key for key, value in cohort_last_dict.items() if value == 2]
c1_list = [key for key, value in cohort_last_dict.items() if value == 1]
cohort_toggle = []
print('* adding cohort id values for each patient/encounter row')
pat_dfp = pat_dfp.reset_index(drop=True)
for indx, row in pat_dfp.iterrows():
    if row['patient_id'] in c3_list:
        cohort_toggle.append(3)
    elif row['patient_id'] in c2_list:
        cohort_toggle.append(2)
    else:
        cohort_toggle.append(1)  # makes the assumption that in the end there will be no cohort 4 items
pat_dfp.insert(loc=10, column='cohort', value=cohort_toggle)
column_list = ['patient_id', 'encounter_id', 'code_system', 'code', 'principal_procedure_indicator', 'date',
               'derived_by_TriNetX', 'source_id', 'rc', 'cohort']
pat_dfp[column_list].to_csv(new_filename, encoding='utf-8', index=False)  # new file
print('* new cohort file created: ', new_filename)
print('* steps 8 & 9 complete')

# * step 10 - data mining part 1
print('* step 10 - determining patient information up till or at radical cystectomy')
# * col: patient rc date, birth year, race, ethnicity, bmi, death year, {16} diagnoses, {3} medications, {1} diagnosis
# pats_uniq = pat_dfp['patient_id'].unique()
print('* step 10 complete')

# * step 11 - data mining part 2
print('* step 11 - determining patient outcomes across their lifetimes')
# * cols:
print('* step 11 complete')

# * step 12 - data mining part 3
print('* step 12 - determining patient outcomes during or after radical cystectomy')
# * cols:
print('* step 12 complete')

# test print toggles
pd.set_option('display.max_columns', None)  # set in order to display all columns in dataframe when test printing
pd.set_option('display.max_rows', None)  # set in order to display all rows in dataframe when test printing

# * timing
print("* %s seconds" % (time.time() - begin))
