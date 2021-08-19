* PROJECT: Community revealed preferences and the proxy means test
* BY: Lendie Follett and Heath Henderson
* DATE: 8/17/21
* LOCATION: /Users/hendersonhl/Documents/Articles/Hybrid-Targeting/Data/Burkina Faso/Cleaning
* PURPOSE: To clean the Hillebrecht et al. data

*** Set up
clear all
cd "/Users/hendersonhl/Documents/Articles/Hybrid-Targeting/Data/Burkina Faso/Cleaning"
global data /Users/hendersonhl/Documents/Articles/Hybrid-Targeting/Data/Burkina Faso/Hillebrecht Files/data
use "$data/data_hhsurvey_and_ranking.dta"

*** Replicate Hillebrecht's PMT
*** Note: The results for Hillebrecht's PMT are presented in the spreadsheet
*** "TAB_A3_weights_econometric_pmt.csv" in the "results" folder.

* Generate indicator for expenditure variable
* Note: This is the outcome variable for the PMT. It's not clear why 
* Hillebrecht et al. dichotomize the expenditure variable.
bysort year trg_unit: egen quota = sum(elig_fin) // elig_fin is the treatment indicator
sort year trg_unit expd identifier  // identifier is the numeric ID variable
by year trg_unit: gen rank = _n  // genderates rank based on expd
gen elig_expend = (rank <= quota)  

* Run PMT regressions
* Note: These results are indentical to those in the spreadsheet.
tab trg_unit, gen(trg)
global xvars hhm_level2	hhm_level3 hhm_level4 hhh_literate hhh_no_primagr ///
ethn_min wat_gd wat_cha san_alt roof_gd wall_gd floor_gd toil_gd nmb_rooms ///	
hhm_d_cart hhm_d_plow hhm_d_bike hhm_d_mbike hhm_d_car hhm_d_radio hhm_d_tv ///
hhm_d_fridge hhm_d_kitchen hhm_d_horse_donkey hhm_d_goat_sheep hhm_d_chicken ///
hhm_d_bullock hhm_d_pig hha_age1660	hha_age60 hhh_married hhh_nwidow	hhh_male
reg elig_expend $xvars trg1-trg58 if year==2008 & rural==1, vce(cluster hh_id) 
reg elig_expend $xvars trg1-trg58 if year==2008 & rural==0, vce(cluster hh_id) 
drop trg1-trg62

*** Data cleaning

* Keep select variables
order hh_id year trg_unit rural elig_fin infA infB infC expd $xvars
keep hh_id-hhh_male

* Renaming
rename hh_id hhid
rename trg_unit community
rename elig_fin treated
rename infA informant1
rename infB informant2
rename infC informant3
rename expd consumption
rename hhm_level2 primary	
rename hhm_level3 secondary
rename hhm_level4 tertiary
rename hhh_literate literacy
rename hhh_no_primagr agriculture
rename ethn_min minority
rename wat_gd well
rename wat_cha water
rename san_alt waste
rename roof_gd roof
rename wall_gd walls
rename floor_gd floor
rename toil_gd toilet
rename nmb_rooms rooms
rename hhm_d_cart cart
rename hhm_d_plow plow
rename hhm_d_bike bike  
rename hhm_d_mbike motorbike
rename hhm_d_car car
rename hhm_d_radio radio
rename hhm_d_tv tv
rename hhm_d_fridge fridge
rename hhm_d_kitchen kitchen
rename hhm_d_horse_donkey horse
rename hhm_d_goat_sheep goat
rename hhm_d_chicken chicken
rename hhm_d_bullock bullock 
rename hhm_d_pig pig
rename hha_age1660 age1660
rename hha_age60 age60
rename hhh_married married
rename hhh_nwidow widow
rename hhh_male male

* Flip select indicators
replace agriculture = (agriculture==0)
replace widow = (widow==0)

* Labels
label var hhid "Household identifier"
label var community "Community identifier"
label var rural "HH lives in rural area"
label var treated "HH eligible for discount"
label var informant1 "Rank given by informant 1"
label var informant2 "Rank given by informant 2"
label var informant3 "Rank given by informant 3"
label var consumption "Monthly per capita consumption (CFA)"
label var primary "Any HH member with primary education"
label var secondary "Any HH member with secondary education"
label var tertiary "Any HH member with tertiary education"
label var literacy "Household head literate"
label var agriculture "HH head occupation is agricultural"
label var minority "HH is ethnic minority"
label var well "HH uses running water or good wells"
label var water "Drinking water is changed at least every other day"
label var waste "Wastewater by cesspool gutters or septic tank"
label var roof "Roof is concrete, metal sheets, or tile"
label var walls "Walls are not mud or straw"
label var floor "Floor is made of cement"
label var toilet "Toilet not in open field"
label var rooms "Number of rooms"
label var cart "HH owns at least one cart"
label var plow "HH owns at least one plow"
label var bike "HH ownse at least one bike"
label var motorbike "HH owns at least one motorbike"
label var car "HH owns at least one car"
label var radio "HH owns at least one radio"
label var tv "HH owns at least one TV"
label var fridge "HH owns at least one fridge"
label var kitchen "HH has kitchen"
label var horse "HH owns at least one horse or donkey"
label var goat "HH owns at least one goat or sheep"
label var chicken "HH owns at least one chicken"
label var bullock "HH owns at least one bullock"
label var pig "HH owns at least one pig"
label var age1660 "Share of HH members between 16 and 60 years"
label var age60 "Share of HH members above 60 years"
label var married "HH head is married"
label var widow "HH head is widowed"
label var male "HH head is male"

* Re-order variables
order hhid-consumption rooms floor walls roof toilet well water waste male ///
married widow age1660 age60 minority primary secondary tertiary literacy ///
agriculture cart-pig

*** Save datasets

* Save final dataset w/ missing values
save "hillebrecht(missing).dta", replace
outsheet using "hillebrecht(missing).csv", comma nolabel replace

* Save final dataset w/o missing values
drop if walls == .a | roof == .a | toilet == .a
drop if informant1 == . | informant2 == . | informant3 == . 
save "hillebrecht.dta", replace
outsheet using "hillebrecht.csv", comma nolabel replace
clear all

* Save variable names and labels
describe, replace clear
keep name varlab
rename name Name
rename varlab Definition
outsheet using "variables.csv", comma nolabel replace
clear all





