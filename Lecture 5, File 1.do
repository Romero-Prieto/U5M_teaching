///File 1: Data preparation/// 
local      pATh         = "/Users/lshjr3/Documents/MRC_U5M_Workshop"         /*Adjust path*/

tempfile   tempBR                                                            /*All tempfiles must be declared before using them.*/
tempfile   tempIR
tempfile   temp
tempfile   mAStEr

local      country      = "GM NG ML"                                         /*DHS country identifier.*/
local      DHS_GM       = "81 61"                                            /*DHS file identifier, The Gambia.*/
local      DHS_NG       = "7B 6A 53 4B 21"                                   /*DHS file identifier, Nigeria.*/
local      DHS_ML       = "7A 6A 53 41 32 01"                                /*DHS file identifier, Mali.*/
local      vARsBR       = "caseid v001 v002 v003 v005 v008 v016 v011 v012 v023 v024 v025 bidx b3 b4 b6 b17 v101 v102" /*Variables from Birth Recode.*/
local      vARsIR       = "caseid v001 v002 v003 v005 v008 v016 v011 v012 v023 v024 v025 v018 v101 v102 vcal_1" /*Variables from Individual Recode.*/

clear
clear mata
generate   survey       = ""
save      `mAStEr', replace                                                  /*Tempfile used to append information from each DHS.*/

use       `vARsBR' in 1/1 using "`pATh'/MLBR7AFL.DTA", clear                 /*Most updated DHS. This file is created to make room for all variales.*/
drop       in 1                                                              /*Makes an empty file, just with variable names.*/
save      `tempBR', replace                                                  /*Saves the file as a temporary file.*/

use       `vARsIR' in 1/1 using "`pATh'/MLIR7AFL.DTA", clear                 /*Follows the same approach with the IR file.*/
drop       in 1
save      `tempIR', replace

foreach co of local country {                                                /*Code is repeated for every country.*/                                      
	foreach survey of local DHS_`co' {                                       /*Code is repeated for every survey.*/
		cls
		use        caseid v0* b* using "`pATh'/`co'BR`survey'FL.DTA", clear  /*Imports all variables starting by v0 and b.*/
		duplicates drop                                                      /*DHS samples do not have duplicates. This is only to show the appropriate place for this cleaning step.*/
		merge 1:1  caseid v002 v012 bidx using `tempBR', nogenerate noreport nolabel keep(master) /*Makes room for the variables that are in the list but not in this specific sample.*/
		keep      `vARsBR'                                                   /*Keeps just the variables of the list.*/
		generate   mother         = 1                                        /*Creates an identifier of women in the BR. They should be mothers or mothers to be.*/
		save      `temp', replace                                            /*Saves date as a temporary file to be used later.*/
		
		use        caseid v* b* using "`pATh'/`co'IR`survey'FL.DTA", clear   /*Imports all variables starting by v and b. IR includes all women aged 15–49, regardless of maternity.*/ 
		duplicates drop                                                      /*Makes a hypothetical cleaning.*/
		merge 1:1  caseid v002 v012 using `tempIR', nogenerate noreport nolabel keep(master) /*Brings the name of thos variables not in the sample.*/
		keep      `vARsIR'                                                   /*Keeps just the variables of the list.*/
		merge 1:m  caseid v001 v002 v003 v008 v016 v008 v012 v023 v024 v025 using `temp', nogenerate /*Combines the IR and the BR into a single file.*/
		
		sort       caseid bidx                                               /*Assumes a natural sort. The code of the mother and the code of each birth.*/
		bysort     caseid: generate   k = _n                                 /*Creates an index. One number per record of a woman.*/
		bysort     caseid: generate   K = _N                                 /*Creates an index. Total number of records per woman.*/

		generate   cluster        = v001                                     /*DHS ids for cluster, household, and respondent.*/
		generate   HH             = v002
		generate   respondent     = v003
		generate   survey         = "`co'`survey'"                           /*Generates a indentifier for each survey.*/
		
		generate   DOB            = mdy(v011 - floor((v011 - 1)/12)*12,1,floor((v011 - 1)/12) + 1900) /*Generates the date of birth of the woman/mother.*/
		recode     v016        (. = 1)                                       /*Exact date only available for the most recent surveys. If not avialable, the first of the month is assumed.*/
		generate   interview      = mdy(v008 - floor((v008 - 1)/12)*12,v016,floor((v008 - 1)/12) + 1900) /*Generates the date of the interview.*/
		
		/*Date of Birth - Birth Histories*/
		sort       b3 b17
		if b17[1] == . {
			/*The exact date is only avaliable for the most recent surveys. If not reported, a random day of the month is assumed. Identifies the limits of this random date.*/
			generate   B_min          = min(mdy(b3 - floor((b3 - 1)/12)*12,1,floor((b3 - 1)/12) + 1900),interview)   if b3    != . /*The first day of the reported month is the lower limit.*/
			generate   B_max          = min(mdy(b3 + 1 - floor(b3/12)*12,1,floor(b3/12) + 1900),interview)           if b3    != . /*The first day of the following month is the upper limit.*/
			}
		else {
			generate   B_min          = min(mdy(b3 - floor((b3 - 1)/12)*12,b17,floor((b3 - 1)/12) + 1900),interview) if b3    != . /*Calculates the date of birth.*/
			replace    b17            = b17 - 1                                                                      if B_min == . & b17 != ./*Dates on the 31th of months with actually 30.*/
			replace    B_min          = min(mdy(b3 - floor((b3 - 1)/12)*12,b17,floor((b3 - 1)/12) + 1900),interview) if B_min == . & b3  != ./*Calculates again, if incorrect date.*/ 
			generate   B_max          = B_min                                /*Because exact dates are available.*/
			}
		
		format     %tdDD/NN/CCYY interview B_* DOB		                     /*Gives date format to the date variables.*/
		rename     b4 sex                                                    /*Identifies the sex and age of the child.*/
		rename     v012 age
		generate   W              = v005/1000000                             /*Identifies the sampling weights, rounded and *10^6. Useful when decimals are not available, not the case.*/
		
		/*Ages at death (in days) - Birth Histories*/
		generate   D_min          = b6 - 100                                        if b6    != .   & b6 <= 200 /*The report is in days within the first month of life.*/
		generate   D_max          = D_min + 1                                       if b6    != .   & b6 <= 200 /*Assumes a plausible maximum of one additional day.*/
		replace    D_min          = 0                                               if b6    == 198 | b6 == 199 /*If days were reported but a number was not provided (rare).*/ 
		replace    D_max          = 365.25/12                                       if b6    == 198 | b6 == 199 /*Max and min bound the first month of life.*/	
		replace    D_min          = (b6 - 200)*365.25/12                            if b6    != .   & b6 >= 200 & b6  < 300 /*The report is 2-24 months.*/
		replace    D_max          = (b6 - 200 + 1)*365.25/12                        if b6    != .   & b6 >= 200 & b6  < 300 /*Assumes a maximum of one additional month.*/
		replace    D_min          = 0                                               if b6    == 298 | b6 == 299 /*If months were reported but a number was not provided (rare).*/ 
		replace    D_max          = 24*365.25/12                                    if b6    == 298 | b6 == 299 /*Max and min bound the first 2 years of life.*/
		replace    D_min          = (b6 - 300)*365.25                               if b6    != .   & b6 >= 300 /*The report is in years after the second birthday.*/
		replace    D_max          = (b6 - 300 + 1)*365.25                           if b6    != .   & b6 >= 300 /*Assumes a plausible maximum of one additional year.*/
		replace    D_min          = 0                                               if b6    == 398 | b6 == 399 /*If years were reported but a number was not provided (rare).*/ 
		replace    D_max          = max(year(interview) - year(B_min),0)*365.25     if b6    == 398 | b6 == 399 /*Age at death could be from 0 to the age at interview.*/
		
		sort       caseid bidx                                                                                  /*Assumes a natural sort.*/
		generate   temp           = 1                                               if B_min != .               /*Creates a temp variable to indicate a child is born.*/
		bysort     caseid: egen       Born  = sum(temp)                                                         /*Creates a Summary of children ever born.*/ 
		replace    Born           = .                                               if k     != 1               /*Constraints the information to be one total per woman.*/
		recode     temp        (1 = .)                                              if b6    != .               /*Constratints the temp variable to indicate a child is alive.*/
		bysort     caseid: egen       Alive = sum(temp)                                                         /*Creates a Summary of children alive.*/ 
		replace    Alive          = .                                               if k     != 1               /*Constraints the information to be one total per woman.*/
		
		/*Pregnancy Calendar, if available*/
		generate   CAL            = "C" + vcal_1                                    if k     == 1   | k  == .   /*Information related to the pregnancy calendar.*/
		generate   row            = v018                                            if k     == 1   | k  == .   /*Useful information to calculate perinatal mortality.*/
		
		keep       survey cluster HH respondent age sex interview B_* D_* W bidx caseid DOB mother Born Alive CAL row /*Retains the relevent variables for U5M.*/
		order      survey cluster HH respondent age sex interview B_* D_* W bidx caseid DOB mother Born Alive CAL row /*Arranges the dataset in this particular order.*/
		sort       survey cluster caseid bidx                                /*Sorts the observation following this particular structure.*/
		append     using `mAStEr'                                            /*Appends the dataset to other surveys previously processed.*/
		save      `mAStEr', replace                                          /*Saves the progress as a temporary file.*/
		}
	}
save     "`pATh'/Results/DHSdata.dta", replace


///Quality Assessemnts///
local      pATh         = "/Users/lshjr3/Documents/MRC_U5M_Workshop"         /*Adjust path*/
use      "`pATh'/Results/DHSdata.dta", clear                                 /*Loads DHS prepared data.*/

bysort survey: egen       date = min(interview)                              /*Creates an anlytical threshold for each survey.*/
generate   sample       = 1 if B_max >= mdy(month(date),day(date),year(date) - 5) /*Creates an analytical (arbitrary) period of 5 years.*/
generate   female       = 1 if sample == 1 & sex == 2                        /*Creates an indicator for female births within the analytical window.*/
generate   male         = 1 if sample == 1 & sex == 1                        /*Creates an indicator for male births within the analytical window.*/
collapse  (sum) female male [aw = W], by(survey)
generate   SRB          = male/female


local      pATh         = "/Users/lshjr3/Documents/MRC_U5M_Workshop"         /*Adjust path*/
use      "`pATh'/Results/DHSdata.dta", clear                                 /*Loads DHS prepared data.*/
bysort survey: egen       date = min(interview)                              /*Creates an anlytical threshold for each survey.*/
generate   sample       = 1 if B_max + (D_min + D_max)/2 >= mdy(month(date),day(date),year(date) - 5) /*Creates an analytical (arbitrary) period of 5 years.*/

local      lISt         = "GM61 GM81 ML01 ML32 ML41 ML53 ML6A ML7A NG21 NG4B NG53 NG6A NG7B" /*Declares the surveys to be analysed.*/
foreach survey of local lISt {
	generate   `survey' = sample if survey == "`survey'"
	}

drop if    D_min       == .   	
collapse (sum) `lISt' [aw = W], by(D_min)                                    /*Gets the number of deaths by age, for each survey.*/
generate   lABel        = string(floor(D_min))           + "d" if D_min  < 365.25*2/12 /*Creates age labels.*/
replace    lABel        = string(floor(D_min/365.25*12)) + "m" if lABel == "" & D_min <= 365.25*5
replace    lABel        = string(floor(D_min/365.25))    + "y" if lABel == ""
order      lABel
twoway    (line ML01 D_min if D_min < 28) (line ML32 D_min if D_min < 28)    /*Plots the first two surveys of Mali.*/
twoway    (line ML6A D_min if D_min < 28) (line ML7A D_min if D_min < 28)    /*Plots the last two surveys of Mali.*/
