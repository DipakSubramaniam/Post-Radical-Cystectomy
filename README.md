# Post-Radical-Cystectomy
Data analysis of large patient/encounter records into usable format for statistical investigation

start - reading in procedure and encounter csv files
step 1 - dropping encounters without an end date in encounter csv file
steps 2 & 3 - finding date of radical cystectomy per patient
step 4 - deleting all encounters per patient that have a start date after the end date of radical cystectomy
step 5 - copying procedure file into cohort file
step 6 - deleting patients who fall under maybe ovaries/fallopian_tubes categories before radical cystectomy
step 7 - finding frequency of each unique code system
steps 8 & 9 - determining patient cohorts
step 9.5 - preparing data structures and shared data headers
step 10 - determining patient information up till or at radical cystectomy
step 11 - determining patient outcomes across their lifetimes
step 12 - determining patient outcomes during or after radical cystectomy
step 13 - gathering requested variable data into single excel file
step 14 - creating data plots
step 15 - propensity matching (TBD)

