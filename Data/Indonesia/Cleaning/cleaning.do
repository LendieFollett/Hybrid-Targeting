* PROJECT: Community revealed preferences and the proxy means test
* BY: Lendie Follett and Heath Henderson
* DATE: 8/17/21
* LOCATION: /Users/hendersonhl/Documents/Articles/Hybrid-Targeting/Data/Indonesia/Cleaning
* PURPOSE: To clean Alatas et al. data

*** Set up
clear all
cd "/Users/hendersonhl/Documents/Articles/Hybrid-Targeting/Data/Indonesia/Cleaning"
global data /Users/hendersonhl/Documents/Articles/Hybrid-Targeting/Data/Indonesia/Alatas Files/codeddata
global table12 /Users/hendersonhl/Documents/Articles/Hybrid-Targeting/Data/Indonesia/Alatas Files/Output/Tables/Table12
global raw /Users/hendersonhl/Documents/Articles/Hybrid-Targeting/Data/Indonesia/Alatas Files/Data/baseline (english)

*** Get X variables
* Note: Alatas et al. use these variables in Table 11
use "$data/analysis01.dta"
areg RANK RANKCONSUMPTION pcfloor tfloor twall toilet water lighting troof ///
fcook house credit hhsize hhsize2 hhage hhage2 hhmale hhmarried hhmalemarr ///
hhsector1 hhsector2 hhsector3 formal informal hheduc2 hheduc3 hheduc4 age04 ///
eschild jschild sschild higheduc2 higheduc3 higheduc4 depratio se2_115 ///
se2_117 se2_13_15 se2_11_12 se2_14 se2_114 se2_18 se2_17 se2_124 se2_125 ///
se2_126_127 se2_119 se2_122_123 se5_110 se5_11_12, ///
absorb(hhid_who_rank) cluster(hhea) // Replicate column 2 of Table 11
order hhid pcfloor tfloor twall toilet water lighting troof fcook house credit ///
hhsize hhsize2 hhage hhage2 hhmale hhmarried hhmalemarr hhsector1 hhsector2 ///
hhsector3 formal informal hheduc2 hheduc3 hheduc4 age04 age14 ///
eschild jschild sschild higheduc2 higheduc3 higheduc4 depratio se2_115 se2_117 ///
se2_13_15 se2_11_12 se2_14 se2_114 se2_18 se2_17 se2_124 se2_125 se2_126_127 ///
se2_119 se2_122_123 se5_110 se5_11_12 
keep hhid-se5_11_12
collapse pcfloor-se5_11_12, by(hhid)
sort hhid
save xvars.dta, replace

* Merge in elite connectedness variable
use "$data/mistargeting_CORRECTED.dta", clear
keep hhid_suseti elite_connectedness
rename hhid_suseti hhid
merge 1:1 hhid using xvars.dta
keep if _merge==3
drop _merge
order hhid pcfloor-hhsector3 elite_connectedness
rename elite_connectedness connected
sort hhid
save xvars.dta, replace

*** Get ranking variable
* Note: Alatas et al. use these variables in Table 12. The variable TREATRANK_2
* is the percentile ranks from the community meetings, even for hybrid
* communities. The variable TREATRANK_1 adds to this variable the relevant
* PMT rankings, in some cases replacing the community rankings done in hybrid
* communities. The variable RANK_exp, which is used in the regression below, 
* is just a re-expression of TREATRANK_1. 
clear all
use "$table12/Table12data.dta"
areg RANK_exp LOGCONSUMPTION loghhsize sharekids hh_eleorless elite_connectedness ///
ETHICMINORITY RELIGIOUSMINORITY WIDOW disabled death sick loss tobacco_alcohol ///
total_savings share_banksavings share_debt connected fam_outside ///
relig_participate projects_work projects_money  if COMMUNITY==1, ///
absorb(kecagroup) cluster(hhea) // Replicate column 5 of Table 12
keep hhid_suseti TREATRANK_2
rename hhid_suseti hhid
rename TREATRANK_2 rank
sort hhid
save ranks.dta, replace

*** Get food consumption data
* Note: For food expenditure, there are only a very small number of missing
* values in the various food categories. We thus follow Alatas et al. and 
* implicitly set those missing values to zero by using rowtotal.
clear all
use "$raw/English_hh_ks1.dta"
foreach var of varlist ks01_1a-ks01_14b {
    replace `var'=. if `var'>999990 
}
egen food = rowtotal(ks01_1a-ks01_13b) // Food expenditure, excluding alcohol and tobacco
replace food=. if food==0  // Three obs. oddly report no food consumption
gen nonstaple = 1 - (ks01_1a + ks01_1b + ks01_2)/food // Non-staple consumption
gen cereals=1 if (ks01_1a>0 & ks01_1a!=.) | (ks01_1b>0 & ks01_1b!=.) // Dietary diversity
replace cereals=0 if ks01_1a==0 & ks01_1b==0
gen tubers=1 if ks01_2>0 & ks01_2!=.
replace tubers=0 if ks01_2==0
gen vegetables=1 if ks01_6>0 & ks01_6!=.
replace vegetables=0 if ks01_6==0
gen fruits=1 if ks01_8>0 & ks01_8!=.
replace fruits=0 if ks01_8==0
gen meats=1 if ks01_4>0 & ks01_4!=.
replace meats=0 if ks01_4==0
gen eggs=1 if ks01_5a>0 & ks01_5a!=.
replace eggs=0 if ks01_5a==0
gen seafood=1 if (ks01_3a>0 & ks01_3a!=.) | (ks01_3b>0 & ks01_3b!=.)
replace seafood=0 if ks01_3a==0 & ks01_3b==0
gen legumes=1 if ks01_7>0 & ks01_7!=.
replace legumes=0 if ks01_7==0
gen milk=1 if ks01_5b>0 & ks01_5b!=.
replace milk=0 if ks01_5b==0
gen oils=1 if ks01_9>0 & ks01_9!=.
replace oils=0 if ks01_9==0
gen sweets=1 if ks01_10>0 & ks01_10!=.
replace sweets=0 if ks01_10==0
gen spices=1 if ks01_11>0 & ks01_11!=.
replace spices=0 if ks01_11==0
gen diversity=cereals + tubers + vegetables + fruits + meats + eggs + seafood ///
+ legumes + milk + oils + sweets + spices
keep hhid food nonstaple diversity
sort hhid
save food.dta, replace

*** Get treatment variables
* Note: Alatas et al. use these variables in Table 2, among others
clear all
use "$data/mistargeting_CORRECTED.dta"
order hhid_suseti prop kab kec desa hhea CONSUMPTION PMT COMMUNITY HYBRID ///
ELITE POOREST10 DAYMEETING RTS 
keep hhid_suseti-RTS
rename hhid_suseti hhid
sort hhid
merge 1:1 hhid using xvars.dta
drop _merge
merge 1:1 hhid using ranks.dta
drop _merge
merge 1:1 hhid using food.dta
drop _merge

*** Verify select data
* Note: The X variables and ranking variable are taken directly from datasets
* used in regressions that have been replicated.
codebook hhid // Total no. of households is 5756
codebook hhea // Total no. of villages is 640

* Compare with Table 1 in Alatas et al.
egen tag = tag(hhea)  // Flag one observation per village
tab COMMUNITY if tag==1
tab HYBRID if tag==1
tab PMT if tag==1
drop tag

* Compare with Table 2 in Alatas et al.
sum CONSUMPTION RTS

*** Miscellaneous cleaning
replace food = (food/hhsize)/1000  // In per capita terms (in thousand Rp)
gen nonfood = 1 - food/CONSUMPTION
order hhid-hhea PMT-RTS rank CONSUMPTION food nonfood nonstaple diversity
replace rank=. if PMT==1  // Only keep community ranking
gen hhsize_ae = (hhsize - (1 - 0.93)*age14)^0.85  // See Olken (2005) for details
drop age14
gen consumption_ae = (CONSUMPTION*hhsize)/hhsize_ae
order hhid-CONSUMPTION consumption_ae food-hhsize2 hhsize_ae

* Renaming and labeling
rename prop province
rename kab district
rename kec subdistrict
rename desa village
rename hhea ea
rename PMT pmt
rename COMMUNITY community
rename HYBRID hybrid
rename ELITE elite
rename POOREST10 poorest10
rename DAYMEETING daymeeting
rename RTS treated
rename CONSUMPTION consumption
rename se2_115 ac
rename se2_117 computer
rename se2_13_15 radio
rename se2_11_12 tv
rename se2_14 dvd
rename se2_114 satellite
rename se2_18 gas
rename se2_17 refrigerator
rename se2_124 bicycle
rename se2_125 motorcycle
rename se2_126_127 auto
rename se2_119 hp
rename se2_122_123 jewelry
rename se5_110 chicken
rename se5_11_12 cow
label var hhid "Household identifier"
label var province "Province code"
label var ea "Enumeration area"
label var elite "Elite subtreatment"
label var poorest10 "10 poorest subtreatment"
label var daymeeting "Day meeting subtreatment"
label var treated "Household received cash"
label var rank "Percentile rank in subvillage"
label var consumption_ae "Adult equivalent consumption (in thousand Rp)"
label var food "Per capita food expenditure (in thousand Rp)"
label var nonfood "Proportion of expenditure on non-food items"
label var nonstaple "Proportion of food consumption from non-staples"
label var diversity "Dietary diversity"
label var pcfloor "Household floor area per capita"
label var tfloor "Not earth floor"
label var twall "Brick or cement wall"
label var toilet "Private toilet"
label var water "Clean drinking water"
label var lighting "PLN electricity"
label var troof "Concrete or corrugated roof"
label var fcook "Cooks with firewood"
label var house "Owns house privately"
label var credit "Has received credit"
label var hhsize "Household size"
label var hhsize2 "Household size squared"
label var hhsize_ae "Adult equivalent household size"
label var hhage "Age of head of household"
label var hhage2 "Age of head of household squared"
label var hhmale "Head of household is male"
label var hhmarried "Head of household is married"
label var hhmalemarr "Head of household is male and married"
label var hhsector1 "Head of household works in agriculture sector"
label var hhsector2 "Head of household works in industry sector"
label var hhsector3 "Head of household works in service sector"
label var connected "Connected to elite households"
label var formal "Head of household works in formal sector"
label var informal "Head of household works in informal sector"
label var hheduc2 "Education attainment of HH head is elementary school"
label var hheduc3 "Education attainment of HH head is junior school"
label var hheduc4 "Education attainment of HH head is senior high school or higher"
label var age04 "Number of children 0-4"
label var eschild "Number of children in elementary school"
label var jschild "Number of children in junior high school"
label var sschild "Number of children in senior high school"
label var higheduc2 "Highest education attainment within HH is elementary school"
label var higheduc3 "Highest education attainment within HH is junior school"
label var higheduc4 "Highest education attainment within HH is senior high or higher"
label var depratio "Total dependency ratio"
label var ac "AC"
label var computer "Computer"
label var radio "Radio/cassette player"
label var tv "TV"
label var dvd "DVD player"
label var satellite "Satellite dish"
label var gas "Gas burner"
label var refrigerator "Refrigerator"
label var bicycle "Bicycle"
label var motorcycle "Motorcycle"
label var auto "Car/minibus/truck"
label var hp "HP"
label var jewelry "Jewelry"
label var chicken "Chicken"
label var cow "Caribou/cow"

*** Save datasets

* Save final dataset w/ missing values
save "alatas(missing).dta", replace
outsheet using "alatas(missing).csv", comma nolabel replace

* Drop incomplete observations
drop if depratio==.
drop if pcfloor==.
drop if toilet==.
drop if jewelry==.
drop if hhage==.
drop if water==.
drop if chicken==.
drop if rank==. & (community==1 | hybrid==1) 
drop if diversity==.
drop if nonfood==.

* Save final dataset w/o missing values
save "alatas.dta", replace
outsheet using "alatas.csv", comma nolabel replace

* Save variable names and labels
describe, replace clear
keep name varlab
rename name Name
rename varlab Definition
outsheet using "variables.csv", comma nolabel replace
clear all





