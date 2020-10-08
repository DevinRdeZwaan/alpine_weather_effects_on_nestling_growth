#####################################################################################
# R Code for:

# de Zwaan, D.R., Drake, A., Greenwood, J.L., and Martin, K. (2020). Timing and intensity of weather events shape nestling development strategies in three alpine breeding songbirds. Frontiers in Ecology and Evolution. doi: 10.3389/fevo.2020.570034

#####################################################################################

### Project goals:
Across taxa, offspring size traits are linked to survival and life-time fitness. Inclement weather can be a major constraint on offspring growth and parental care. Despite the adaptive benefits of larger offspring, we have a limited understanding of the effects of severe environmental conditions across developmental stages and how coping strategies differ among species. We assessed the influence of inclement weather on offspring size and mass traits within populations of three alpine breeding songbirds in British Columbia: 1) horned lark (Eremophila alpestris), 2) dark-eyed junco (Junco hyemalis), and 3) savannah sparrow (Passerculus sandwichensis). Specifically, we investigated at which stages during early-life development offspring are most vulnerable to inclement weather and whether thresholds exist in the developmental response to severe weather events. 

### Version: R 3.6.3

# Table of Contents:

# 1) '01_de Zwaan et al_2020_ Sliding window code'

### Description
Code for sliding window approach to identify key weather variables and time windows. Also, code to extract weather variables for the selected time windows from the raw weather data for further analysis.

# 2) '02_de Zwaan et al_2020_ Model fitting and selection code'

### Description
Code for General Additive Models and model selection of the candidate weather variables. 

