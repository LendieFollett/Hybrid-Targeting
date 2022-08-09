* PROJECT: A hybrid approach to targeting social assistance
* BY: Lendie Follett and Heath Henderson
* DATE: 8/9/22
* LOCATION: /Users/hendersonhl/Documents/Articles/Hybrid-Targeting/Data/Indonesia/Cleaning
* PURPOSE: To clean household-level ranking data

* Set file paths
cd "/Users/hendersonhl/Documents/Articles/Hybrid-Targeting/Data/Indonesia/Cleaning"
global codeddata "/Users/hendersonhl/Documents/Articles/Hybrid-Targeting/Data/Indonesia/Alatas Files/codeddata"
global baseline "/Users/hendersonhl/Documents/Articles/Hybrid-Targeting/Data/Indonesia/Alatas Files/Data/baseline"

* Open data from household-level ranking exercise
use "$baseline/hh_cr04a.dta", clear

* Renaming, dropping, and labeling
rename cr04a_line rank 
rename hhid_cr04 hhid_ranked
order hhea hhid hhid_ranked rank 
keep hhea-rank
label var hhea "Identifier for community"
label var hhid "Identifier for household that is ranking"
label var hhid_ranked "Identifier for household being ranked"
label var rank "Rank of household"
sort hhea hhid rank

* Remove village head rankings
merge m:1 hhid using "$codeddata/intermediate_data/rthead.dta"
drop if RTHEAD==1
keep hhea-rank
rename hhid hhid_ranker
sort hhea hhid_ranker rank

* Export
save "hh_rankings.dta", replace
outsheet using "hh_rankings.csv", comma nolabel replace
clear all
