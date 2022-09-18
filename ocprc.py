# Post Radical Cystectomy Processing Script
# Ovarian Cancer After Radical Cystectomy For Bladder Cancer
# Code by: Dipak Subramaniam
# In coordination with: Pranjal Agrawal
# Version: 0.1
# Start date: August 16, 2022

import time
import pandas as pd
from dask import dataframe as dd
import numpy
from collections import Counter
import warnings
import matplotlib.pyplot as plt

# * timing
begin = time.time()

# * warning suppressor
warnings.filterwarnings('ignore')  # suppresses SettingWithCopyWarning in steps 2 & 3
# * reading procedure and encounter files
print('* reading in procedure and encounter csv files')
dfp = dd.read_csv('procedure.csv')
dfe = dd.read_csv('encounter.csv')
old_dfp = dfp.compute()
patient_tracker = []
print('* ~ num patients: ', len(old_dfp['patient_id'].unique()))
patient_tracker.append(len(old_dfp['patient_id'].unique()))

# * step 1 - dropping encounters if no end date present
print('* step 1 - dropping encounters without an end date in encounter csv file')
dfe_drop = dfe.dropna(subset=['end_date'])  # dropping rows if end nan values (empty spaces) exist for the end date
new_dfe = dfe_drop.compute()
old_dfe = dfe.compute()
rows_dropped = old_dfe.loc[~old_dfe.set_index(list(old_dfe.columns)).index.isin(new_dfe.set_index(list(new_dfe.columns)).index)]
encounters_dropped = rows_dropped['encounter_id'].tolist()
# old_dfp = dfp.compute()
dfp_s1 = old_dfp[~old_dfp['encounter_id'].isin(encounters_dropped)]
print('* updated dropped encounters on procedure csv file')
rows_dropped_pro = old_dfp.loc[~old_dfp.set_index(list(old_dfp.columns)).index.isin(dfp_s1.set_index(list(dfp_s1.columns)).index)]
print('* ~ num patients: ', len(dfp_s1['patient_id'].unique()))
patient_tracker.append(len(dfp_s1['patient_id'].unique()))
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
pat_dfp = dfp_s1
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
date_rc = []
flag_rc = []
dict_rc_count = pat_dfp.to_dict('records')

print('* adding and populating rc column with radical cystectomy flags')
for row in dict_rc_count:
    cur_patient = row['patient_id']
    cur_encounter = row['encounter_id']
    if cur_patient in rc_date_dict:
        date_rc.append(rc_date_dict[cur_patient])
    else:
        date_rc.append('no__date')
    if cur_patient in patient_enc_dict and cur_encounter == patient_enc_dict[cur_patient]:
        flag_rc.append(1)
    else:
        flag_rc.append(0)

pat_dfp.insert(loc=9, column='rc', value=flag_rc)
pat_dfp.insert(loc=10, column='rc_date', value=date_rc)
print('* removed a small fraction of patients wherein their codes do not register a radical cystectomy')
print('* ~ num patients: ', len(pat_dfp['patient_id'].unique()))
patient_tracker.append(len(pat_dfp['patient_id'].unique()))
patient_tracker.append(len(pat_dfp['patient_id'].unique()))
print('* steps 2 & 3 complete')

# * step 4 - deleting all encounters per patient that have a start date after the end date of radical cystectomy
print('* step 4 - deleting all encounters per patient that have a start date after the end date of radical cystectomy')
print('* end dates being tabulated for all radical cystectomies')
indices_to_delete = []
new_dfe = new_dfe.reset_index()
dict_end_enc = new_dfe.to_dict('records')
for row in dict_end_enc:
    cur_patient = row['patient_id']
    cur_start_date = row['start_date']
    cur_end_date = row['end_date']
    if (cur_patient in rc_date_dict) and (cur_start_date == int(rc_date_dict[cur_patient])):
        rc_end_date_dict[cur_patient] = str(int(cur_end_date))

print('* parsing for encounters that happened after the radical cystectomies')
count = 0
pat_dfp = pat_dfp.reset_index()
dict_del_after = pat_dfp.to_dict('records')
for row in dict_del_after:
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

pat_dfp.insert(loc=11, column='dele', value=indices_to_delete)
pd.set_option('display.max_columns', None)
pat_dfp = pat_dfp[pat_dfp.dele != 1]
print('* deleted rows containing encounters with dates after the end of the patient radical cystectomy')
print('* ~ num patients: ', len(pat_dfp['patient_id'].unique()))
patient_tracker.append(len(pat_dfp['patient_id'].unique()))
print('* step 4 complete')

# * step 5 - preparing cohort file
print('* step 5 - copying procedure file into cohort file')
new_filename = 'cohort.csv'
print('* new cohort file set to be created: ', new_filename)
print('* ~ num patients: ', len(pat_dfp['patient_id'].unique()))
patient_tracker.append(len(pat_dfp['patient_id'].unique()))
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
maybe_list = [key for key, value in cohort_rc_dict.items() if value == 1]

print('* identifying patients with maybe codes')
maybe_deletes = []
pat_dfp = pat_dfp.reset_index(drop=True)
dict_maybe = pat_dfp.to_dict('records')
for row in dict_maybe:
    if row['patient_id'] in maybe_list:
        maybe_deletes.append(1)
    else:
        maybe_deletes.append(0)
pat_dfp.insert(loc=12, column='maybe', value=maybe_deletes)
pat_dfp = pat_dfp[pat_dfp.maybe != 1]
print('* deleting patients with maybe codes')
print('* ~ num patients: ', len(pat_dfp['patient_id'].unique()))
patient_tracker.append(len(pat_dfp['patient_id'].unique()))
print('* step 6 complete')

# * step 7 - code system frequencies
print('* step 7 - finding frequency of each unique code system')
# frequency = pat_dfp['code_system'].value_counts()
frq_values = pat_dfp['code_system'].value_counts().keys().tolist()
frq_counts = pat_dfp['code_system'].value_counts().tolist()
print("* unique values in column 'code_system' : ", frq_values)
print("* frequency of values in column 'code_system' : ", frq_counts)
print('* ~ num patients: ', len(pat_dfp['patient_id'].unique()))
patient_tracker.append(len(pat_dfp['patient_id'].unique()))
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
gbpe3 = pat_dfp.groupby(['patient_id', 'encounter_id'])
cpe = []  # codes per encounter
cspe = []  # code systems per encounter

print('* searching all patient-encounter pairs for cohort classification')
for pe, grp in gbpe3:  # searching all the patient-encounter pairs and each pair's data containing the other columns
    cpe = grp['code']
    cspe = grp['code_system']
    also_in_three = False  # variable that is toggled if codes are in cohort 3 and 2, prevents cohort 2 from claiming it
    if cohort_last_dict[pe[0]] == 0:  # checking if patient has already been accounted for
        if 'HCPCS' in cspe.tolist() or 'ICD-9-CM' in cspe.tolist():
            cohort_last_dict[pe[0]] = 4  # 4 = non cpt/icd-10 code patients
        else:
            cohort_last_dict[pe[0]] = 1  # everyone else gets a 1
            if set(cpe).intersection(coh3cpt) or set([cd[:4] for cd in cpe]).intersection(coh3icd10c) or (
                    set([cd[:4] for cd in cpe]).intersection(coh3icd10a) and set([cd[:4] for cd in cpe]).intersection(
                    coh3icd10b)):
                cohort_last_dict[pe[0]] = 3
            # if codes also intersect with cohort 3 and 2 designated codes, --> cohort 3
            if (set(cpe).intersection(coh3cpt) or set([cd[:4] for cd in cpe]).intersection(coh3icd10c) or (
                    set([cd[:4] for cd in cpe]).intersection(coh3icd10a) and set([cd[:4] for cd in cpe]).intersection(
                    coh3icd10b))) and (
                    set([cd[:4] for cd in cpe]).intersection(coh2cpt) or set([cd[:4] for cd in cpe]).intersection(
                    coh2icd10c) or (set([cd[:4] for cd in cpe]).intersection(coh2icd10a) and set(
                    [cd[:4] for cd in cpe]).intersection(coh2icd10b))):
                cohort_last_dict[pe[0]] = 3
                also_in_three = True  # will prevent next if statement from being processed
            # if codes also intersect with cohort 2 designated codes --> cohort 2
            if (set(cpe).intersection(coh2cpt) or set([cd[:4] for cd in cpe]).intersection(coh2icd10c) or (
                    set([cd[:4] for cd in cpe]).intersection(coh2icd10a) and set([cd[:4] for cd in cpe]).intersection(
                    coh2icd10b))) and not also_in_three:
                cohort_last_dict[pe[0]] = 2
patient_tracker.append(len(pat_dfp['patient_id'].unique()))

# * cohort analysis
res = Counter(cohort_last_dict.values())
print("* total patients: ", len(cohort_last_dict))
print('* cohort breakdown: ', res)

# * further code for part 8 and subsequent conversion to updated cohort csv excel file
c4_list = [key for key, value in cohort_last_dict.items() if value == 4]
c3_list = [key for key, value in cohort_last_dict.items() if value == 3]
c2_list = [key for key, value in cohort_last_dict.items() if value == 2]
c1_list = [key for key, value in cohort_last_dict.items() if value == 1]
print("* total patients who have icd-9/hcpcs (cohort 4): ", len(c4_list))
print("* total patients in cohort 1: ", len(c1_list))
print("* total patients in cohort 2: ", len(c2_list))
print("* total patients in cohort 3: ", len(c3_list))
cohort_toggle = []

print('* adding cohort id values for each patient/encounter row')
pat_dfp = pat_dfp.reset_index(drop=True)
dict_coh = pat_dfp.to_dict('records')
for row in dict_coh:
    if row['patient_id'] in c3_list:
        cohort_toggle.append(3)
    elif row['patient_id'] in c2_list:
        cohort_toggle.append(2)
    elif row['patient_id'] in c1_list:
        cohort_toggle.append(1)
    else:
        cohort_toggle.append(4) # assumption that cohort 4 will be excluded

pat_dfp.insert(loc=13, column='cohort', value=cohort_toggle)
pat_dfp = pat_dfp[pat_dfp.cohort != 4]
column_list = ['patient_id', 'encounter_id', 'code_system', 'code', 'principal_procedure_indicator', 'date',
               'derived_by_TriNetX', 'source_id', 'rc', 'rc_date', 'cohort']
pat_dfp[column_list].to_csv(new_filename, encoding='utf-8', index=False)  # new file
print('* new cohort file created: ', new_filename)
print('* ~ num patients: ', len(pat_dfp['patient_id'].unique()))
patient_tracker.append(len(pat_dfp['patient_id'].unique()))
print('* steps 8 & 9 complete')

# * step 9.5 - data mining pre-work
print('* step 9.5 - preparing data structures')
column10 = ['pat_id', 'rc_date', 'cohort', 'encs_per_patient', 'yob', 'race', 'ethnicity', 'bmi', 'yod', 'fh_breast', 'ph_breast',
            'fh_ovarian', 'ph_ovarian', 'ph_colorectal', 'lynch_syndrome', 'hereditary_breast&ovarian_lf', 'tobacco',
            'hypertensive', 'diabetes', 'hyperlipidemia', 'ischemic_hd', 'overweight/obese', 'polycystic_od/s',
            'endometriosis', 'pid', 'systemic_contraceptives', 'estrogens', 'progestins', 'hysterectomy']
column11 = ['pat_id', 'rc_date', 'cohort', 'encs_per_patient', 'dementia', 'cog_impair', 'parkinsons', 'osteoperosis', 'ischemic_hd',
            'myocard_infarc', 'aff_mood_disorders', 'anxiety', 'depression', 'fem_sex_dysf', 'stroke',
            'deep_vein_thromb', 'pulm_emb', 'genital_prolapse']
column12 = ['pat_id', 'rc_date', 'cohort', 'encs_per_patient', 'dementia', 'cog_impair', 'parkinsons', 'osteoperosis', 'ischemic_hd',
            'myocard_infarc', 'aff_mood_disorders', 'anxiety', 'depression', 'fem_sex_dysf', 'stroke',
            'deep_vein_thromb', 'pulm_emb', 'genital_prolapse', 'ovarian_cancer', 'peritoneal_cancer', 'breast_cancer']
data_set_10 = {column: [] for column in column10}
data_set_11 = {column: [] for column in column11}
data_set_12 = {column: [] for column in column12}
cohort_uniq2 = pat_dfp['patient_id'].unique()
cohort_patient_flags = dict.fromkeys(cohort_uniq_patients, 0)
s10_11_12 = 0
coh_dict = pat_dfp.to_dict('records')

# * dataframe where rc is flagged and actually has a date
for row in coh_dict:
    if row['rc'] == 1 and row['rc_date'] != 'no__date' and cohort_patient_flags[row['patient_id']] == 0:
        cohort_patient_flags[row['patient_id']] += 1  # flag ensuring each valid patient with an rc is counted only once
        data_set_10['pat_id'].append(row['patient_id'])
        data_set_11['pat_id'].append(row['patient_id'])
        data_set_12['pat_id'].append(row['patient_id'])
        data_set_10['rc_date'].append(row['rc_date'])
        data_set_11['rc_date'].append(row['rc_date'])
        data_set_12['rc_date'].append(row['rc_date'])
        data_set_10['cohort'].append(row['cohort'])
        data_set_11['cohort'].append(row['cohort'])
        data_set_12['cohort'].append(row['cohort'])
        s10_11_12 += 1

# * setting the rest of the data to empty lists of each lists' size
special_cols = ['pat_id', 'rc_date', 'cohort', 'encs_per_patient']
for key in data_set_10:
    if key not in special_cols:
        data_set_10[key] = [None] * s10_11_12
    if key == 'encs_per_patient':
        data_set_10[key] = [0] * s10_11_12
for key in data_set_11:
    if key not in special_cols:
        data_set_11[key] = ['N'] * s10_11_12
    if key == 'encs_per_patient':
        data_set_11[key] = [0] * s10_11_12
for key in data_set_12:
    if key not in special_cols:
        data_set_12[key] = ['N'] * s10_11_12
    if key == 'encs_per_patient':
        data_set_12[key] = [0] * s10_11_12

# * reading in patient, lab_result, diagnosis and medication_drug .csv files
print('* reading patient, vital_signs, diagnosis and medication_drug csv files')
df_pat = pd.read_csv('patient.csv')
dict_pat = df_pat.to_dict('list')
df_vit = pd.read_csv('vitals_signs.csv')
dict_vit = df_vit.to_dict('records')
df_dia = pd.read_csv('diagnosis.csv')
dict_dia = df_dia.to_dict('records')
df_med = pd.read_csv('medication_drug.csv')
dict_med = df_med.to_dict('records')
dict_pro = pat_dfp.to_dict('records')

# * reading in encounters / each patient from encounter .csv file
print('* reading encounter csv file to determine number of encounters per patient')
for row in dict_end_enc:
    if row['patient_id'] in data_set_10['pat_id']:
        pos = data_set_10['pat_id'].index(row['patient_id'])  # positions between sets 10, 11 and 12 are the same
        data_set_10['encs_per_patient'][pos] += 1
        data_set_11['encs_per_patient'][pos] += 1
        data_set_12['encs_per_patient'][pos] += 1

print('* ~ num patients: ', s10_11_12)
patient_tracker.append(s10_11_12)
print('* step 9.5 complete')

# * step 10 - data mining part 1
print('* step 10 - determining patient information up till or at radical cystectomy')
# * col: patient, rc date, birth year, race, ethnicity, bmi, death year, {16} diagnoses, {3} medications, {1} diagnosis
# * for each patient in data set 10, grabbing corresponding patient csv data
print('* for each patient in data set 10, grabbing corresponding patient csv data')
elm_ind = 0
for elm in data_set_10['pat_id']:
    pos = dict_pat['patient_id'].index(elm)
    data_set_10['yob'][elm_ind] = dict_pat['year_of_birth'][pos]
    data_set_10['race'][elm_ind] = dict_pat['race'][pos]
    data_set_10['ethnicity'][elm_ind] = dict_pat['ethnicity'][pos]
    data_set_10['yod'][elm_ind] = dict_pat['month_year_death'][pos]
    elm_ind += 1

# * reformatting birth and death values
data_set_10['yob'] = [str(e) for e in data_set_10['yob']]
data_set_10['yob'] = ['n/a' if v in ['nan'] else v for v in data_set_10['yob']]  # replacing 'nan' with 'n/a'
data_set_10['yob'] = [v[:-2] if '.0' in v else v for v in data_set_10['yob']]  # deleting '.0'
data_set_10['yod'] = [str(e) for e in data_set_10['yod']]
data_set_10['yod'] = ['living' if v in ['nan'] else v for v in data_set_10['yod']]  # replacing 'nan' with 'living'
data_set_10['yod'] = [v[:-4] if '.0' in v else v for v in data_set_10['yod']]  # deleting '##.0'

# * for each patient in data set 10, grabbing corresponding vitals_signs csv data
print('* for each patient in data set 10, grabbing corresponding vitals_signs csv data')
for row in dict_vit:
    if row['units_of_measure'] == 'kg/m2' and row['patient_id'] in data_set_10['pat_id']:
        pos = data_set_10['pat_id'].index(row['patient_id'])
        data_set_10['bmi'][pos] = row['value']  # str(int(row['value']))
data_set_10['bmi'] = ['None' if v is None else v for v in data_set_10['bmi']]  # replacing None with 'n/a'
data_set_10['bmi'] = [str(e) for e in data_set_10['bmi']]
data_set_10['bmi'] = ['n/a' if v in ['nan'] else v for v in data_set_10['bmi']]  # replacing 'nan' with 'n/a'
print('* total patients (', elm_ind, ') - valid bmi values (', data_set_10['bmi'].count('None'), ') = (', (elm_ind-data_set_10['bmi'].count('None')), ')')

# * for each patient in data set 10, grabbing corresponding diagnosis csv data
print('* for each patient in data set 10, grabbing corresponding diagnosis csv data')
# * at the time of cystectomy or up till the time of cystectomy
for row in dict_dia:
    if row['patient_id'] in data_set_10['pat_id'] and row['date'] <= data_set_10['rc_date'][data_set_10['pat_id'].index(row['patient_id'])]:
        pos = data_set_10['pat_id'].index(row['patient_id'])
        if row['code'][:5] == 'Z80.3':
            data_set_10['fh_breast'][pos] = 'Y'
        if row['code'][:5] == 'Z85.3':
            data_set_10['ph_breast'][pos] = 'Y'
        if row['code'][:6] == 'Z80.41':
            data_set_10['fh_ovarian'][pos] = 'Y'
        if row['code'][:6] == 'Z85.43':
            data_set_10['ph_ovarian'][pos] = 'Y'
        if row['code'][:7] == 'Z85.038':
            data_set_10['ph_colorectal'][pos] = 'Y'
        if row['code'][:6] == 'V84.09':
            data_set_10['lynch_syndrome'][pos] = 'Y'
        if row['code'][:6] == 'Z15.01':
            data_set_10['hereditary_breast&ovarian_lf'][pos] = 'Y'
        if row['code'][:5] == 'Z72.0':
            data_set_10['tobacco'][pos] = 'Y'
        if row['code'][:3] in ['I{}'.format(i) for i in range(10, 17)]:
            data_set_10['hypertensive'][pos] = 'Y'
        if row['code'][:3] in ['E{}'.format(i) for i in ["%.2d" % i for i in range(8, 14)]]:
            data_set_10['diabetes'][pos] = 'Y'
        if row['code'][:3] == 'E78':
            data_set_10['hyperlipidemia'][pos] = 'Y'
        if row['code'][:3] in ['I{}'.format(i) for i in range(20, 26)]:
            data_set_10['ischemic_hd'][pos] = 'Y'
        if row['code'][:3] == 'E66':
            data_set_10['overweight/obese'][pos] = 'Y'
        if row['code'][:3] == 'E28':
            data_set_10['polycystic_od/s'][pos] = 'Y'
        if row['code'][:3] == 'N80':
            data_set_10['endometriosis'][pos] = 'Y'
        if row['code'][:3] == 'N73':
            data_set_10['pid'][pos] = 'Y'

# * for each patient in data set 10, grabbing corresponding procedure csv (bmi) data
print('* for each patient in data set 10, grabbing corresponding procedure csv (bmi) data')
for row in dict_pro:
    if row['patient_id'] in data_set_10['pat_id'] and row['date'] <= data_set_10['rc_date'][data_set_10['pat_id'].index(row['patient_id'])]:
        pos = data_set_10['pat_id'].index(row['patient_id'])
        if row['code'] in ['1013911', '1008870'] or row['code'] in [1013911, 1008870]:
            data_set_10['hysterectomy'][pos] = 'Y'

# * for each patient in data set 10, grabbing corresponding medication_drug csv data
print('* for each patient in data set 10, grabbing corresponding medication_drug csv data')
for row in dict_med:
    if row['patient_id'] in data_set_10['pat_id'] and row['start_date'] <= data_set_10['rc_date'][data_set_10['pat_id'].index(row['patient_id'])]:
        pos = data_set_10['pat_id'].index(row['patient_id'])
        if len(str(row['code'])) >= 5:
            if str(row['code'])[:5] == 'HS200':
                data_set_10['systemic_contraceptives'][pos] = 'Y'
            if str(row['code'])[:5] == 'HS300':
                data_set_10['estrogens'][pos] = 'Y'
            if str(row['code'])[:5] == 'HS800':
                data_set_10['progestins'][pos] = 'Y'

# * reformatting Y/N values for last 20 columns in data set 10
for key in data_set_10:
    if key in column10[-20:]:  # targeting all the yes/no variables
        data_set_10[key] = ['N' if v is None else v for v in data_set_10[key]]  # replacing None with 'N'
patient_tracker.append(s10_11_12)
print('* step 10 complete')

# * step 11 - data mining part 2
print('* step 11 - determining patient outcomes across their lifetimes')
# * col: patient, rc date, (14) diagnoses
# * for each patient in data set 11, grabbing corresponding diagnosis csv data
print('* for each patient in data set 11, grabbing corresponding diagnosis csv data')
# * at any point of life
first_11 = dict.fromkeys(data_set_11['pat_id'], [])
for key in first_11:
    first_11[key] = [0] * (len(column11) - len(special_cols))
for row in dict_dia:
    if row['patient_id'] in data_set_11['pat_id']:
        pos = data_set_11['pat_id'].index(row['patient_id'])
        if row['code'][:3] == 'F03' and first_11[row['patient_id']][0] == 0:
            data_set_11['dementia'][pos] = str(row['date'])[:4]
            first_11[row['patient_id']][0] += 1
        if (row['code'][:6] == 'G31.84' or row['code'][:3] == 'F09') and first_11[row['patient_id']][1] == 0:
            data_set_11['cog_impair'][pos] = str(row['date'])[:4]
            first_11[row['patient_id']][1] += 1
        if row['code'][:3] == 'G20' and first_11[row['patient_id']][2] == 0:
            data_set_11['parkinsons'][pos] = str(row['date'])[:4]
            first_11[row['patient_id']][2] += 1
        if row['code'][:3] in ['M80', 'M81'] and first_11[row['patient_id']][3] == 0:
            data_set_11['osteoperosis'][pos] = str(row['date'])[:4]
            first_11[row['patient_id']][3] += 1
        if row['code'][:3] in ['I{}'.format(i) for i in range(20, 26)] and first_11[row['patient_id']][4] == 0:
            data_set_11['ischemic_hd'][pos] = str(row['date'])[:4]
            first_11[row['patient_id']][4] += 1
        if row['code'][:3] in ['I21', 'I22'] and first_11[row['patient_id']][5] == 0:
            data_set_11['myocard_infarc'][pos] = str(row['date'])[:4]
            first_11[row['patient_id']][5] += 1
        if row['code'][:3] in ['F{}'.format(i) for i in range(30, 40)] and first_11[row['patient_id']][6] == 0:
            data_set_11['aff_mood_disorders'][pos] = str(row['date'])[:4]
            first_11[row['patient_id']][6] += 1
        if row['code'][:3] == 'F41' and first_11[row['patient_id']][7] == 0:
            data_set_11['anxiety'][pos] = str(row['date'])[:4]
            first_11[row['patient_id']][7] += 1
        if row['code'][:3] in ['F32', 'F33'] and first_11[row['patient_id']][8] == 0:
            data_set_11['depression'][pos] = str(row['date'])[:4]
            first_11[row['patient_id']][8] += 1
        if row['code'][:3] == 'F52' and first_11[row['patient_id']][9] == 0:
            data_set_11['fem_sex_dysf'][pos] = str(row['date'])[:4]
            first_11[row['patient_id']][9] += 1
        if row['code'][:3] in ['I{}'.format(i) for i in range(60, 64)] and first_11[row['patient_id']][10] == 0:
            data_set_11['stroke'][pos] = str(row['date'])[:4]
            first_11[row['patient_id']][10] += 1
        if row['code'][:5] in ['I82.4', 'I82.6'] and first_11[row['patient_id']][11] == 0:
            data_set_11['deep_vein_thromb'][pos] = str(row['date'])[:4]
            first_11[row['patient_id']][11] += 1
        if row['code'][:3] == 'I26' and first_11[row['patient_id']][12] == 0:
            data_set_11['pulm_emb'][pos] = str(row['date'])[:4]
            first_11[row['patient_id']][12] += 1
        if row['code'][:3] == 'N81' and first_11[row['patient_id']][13] == 0:
            data_set_11['genital_prolapse'][pos] = str(row['date'])[:4]
            first_11[row['patient_id']][13] += 1
patient_tracker.append(s10_11_12)
print('* step 11 complete')

# * step 12 - data mining part 3
print('* step 12 - determining patient outcomes during or after radical cystectomy')
# * col: patient, rc date, (17) diagnoses
# * for each patient in data set 12, grabbing corresponding diagnosis csv data
print('* for each patient in data set 12, grabbing corresponding diagnosis csv data')
# * at the time of cystectomy or after the time of cystectomy
first_12 = dict.fromkeys(data_set_12['pat_id'], [])
for key in first_12:
    first_12[key] = [0] * (len(column12) - len(special_cols))
for row in dict_dia:
    if row['patient_id'] in data_set_12['pat_id'] and row['date'] >= data_set_12['rc_date'][data_set_12['pat_id'].index(row['patient_id'])]:
        pos = data_set_12['pat_id'].index(row['patient_id'])
        if row['code'][:3] == 'F03' and first_12[row['patient_id']][0] == 0:
            data_set_12['dementia'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][0] += 1
        if (row['code'][:6] == 'G31.84' or row['code'][:3] == 'F09') and first_12[row['patient_id']][1] == 0:
            data_set_12['cog_impair'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][1] += 1
        if row['code'][:3] == 'G20' and first_12[row['patient_id']][2] == 0:
            data_set_12['parkinsons'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][2] += 1
        if row['code'][:3] in ['M80', 'M81'] and first_12[row['patient_id']][3] == 0:
            data_set_12['osteoperosis'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][3] += 1
        if row['code'][:3] in ['I{}'.format(i) for i in range(20, 26)] and first_12[row['patient_id']][4] == 0:
            data_set_12['ischemic_hd'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][4] += 1
        if row['code'][:3] in ['I21', 'I22'] and first_12[row['patient_id']][5] == 0:
            data_set_12['myocard_infarc'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][5] += 1
        if row['code'][:3] in ['F{}'.format(i) for i in range(30, 40)] and first_12[row['patient_id']][6] == 0:
            data_set_12['aff_mood_disorders'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][6] += 1
        if row['code'][:3] == 'F41' and first_12[row['patient_id']][7] == 0:
            data_set_12['anxiety'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][7] += 1
        if row['code'][:3] in ['F32', 'F33'] and first_12[row['patient_id']][8] == 0:
            data_set_12['depression'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][8] += 1
        if row['code'][:3] == 'F52' and first_12[row['patient_id']][9] == 0:
            data_set_12['fem_sex_dysf'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][9] += 1
        if row['code'][:3] in ['I{}'.format(i) for i in range(60, 64)] and first_12[row['patient_id']][10] == 0:
            data_set_12['stroke'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][10] += 1
        if row['code'][:5] in ['I82.4', 'I82.6'] and first_12[row['patient_id']][11] == 0:
            data_set_12['deep_vein_thromb'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][11] += 1
        if row['code'][:3] == 'I26' and first_12[row['patient_id']][12] == 0:
            data_set_12['pulm_emb'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][12] += 1
        if row['code'][:3] == 'N81' and first_12[row['patient_id']][13] == 0:
            data_set_12['genital_prolapse'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][13] += 1
        if row['code'][:3] == 'C56' and first_12[row['patient_id']][14] == 0:
            data_set_12['ovarian_cancer'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][14] += 1
        if row['code'][:5] in ['C48.1', 'C48.2'] and first_12[row['patient_id']][15] == 0:
            data_set_12['peritoneal_cancer'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][15] += 1
        if row['code'][:3] == 'C50' and first_12[row['patient_id']][16] == 0:
            data_set_12['breast_cancer'][pos] = str(row['date'])[:4]
            first_12[row['patient_id']][16] += 1
patient_tracker.append(s10_11_12)
print('* step 12 complete')

print('* step 13 - gathering requested variable data into single excel file')
outfile = 'data_set.xlsx'
df10 = pd.DataFrame.from_dict(data_set_10)
df11 = pd.DataFrame.from_dict(data_set_11)
df12 = pd.DataFrame.from_dict(data_set_12)
print('* writing to ', outfile)
with pd.ExcelWriter(outfile) as writer:
    df10.to_excel(writer, sheet_name='Before Or On Encounter_RC')
    df11.to_excel(writer, sheet_name='At Any Time In Life')
    df12.to_excel(writer, sheet_name='During Or After Encounter_RC')
print('* step 13 complete')

print('* step 14 - creating frequency/time plot')
# * histogram for number of rcs over time periods
b_list = [1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2022]
freq_rc, bins, bars = plt.hist([int(str(y)[:4]) for y in df10['rc_date'].tolist()], bins=b_list, edgecolor="orange", color="blue", rwidth=0.95)
plt.xlabel('RC Years')
plt.xticks(b_list, rotation=45)
plt.ylabel('Number of RCs')
plt.title('Frequency of Radical Cystectomies (RC) over Time')
for bar in bars:
    x = (bar._x0 + bar._x1) / 2
    y = bar._y1 + 0.05
    plt.text(x, y, int(bar._y1), ha='center')
plt.tight_layout()
plt.savefig('rc_plot.jpg')
plt.show()

# * chart of number of patients at each step
fig = plt.figure()
ax = fig.add_subplot(111)
x = ['init', '1', '2', '3', '4', '5', '6', '7', '8', '9', '9.5', '10', '11', '12']
y = patient_tracker
plt.xlabel("Steps")
plt.ylabel("Number of Patients")
plt.title("Number of Patients at Each Step")
ax.plot(x, y, '-.')
val = 0
for index in range(len(x)):
    if val == 0:
        ax.text(x[index], y[index] + 0.1, str(y[index]) + " - " + str(x[index]), size=10)
    else:
        if y[index] != y[index-1]:
            ax.text(x[index], y[index] + 0.1, str(y[index]) + " - step " + str(x[index]), size=10)
    val += 1
fig.savefig('num_patients_plot.jpg')
plt.show()
print('* step 14 complete')

print('* step 15 - propensity matching')

print('* step 15 complete')

# test print toggles
pd.set_option('display.max_columns', None)  # set in order to display all columns in dataframe when test printing
pd.set_option('display.max_rows', None)  # set in order to display all rows in dataframe when test printing

# * timing
print("* %s seconds" % (time.time() - begin))
