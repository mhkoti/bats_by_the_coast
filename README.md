# bats_by_the_coast
Material and code for an article: Lundberg et al. 2024 - Contrasting seasonal distribution patterns could reduce interspecific competition in two boreal aerial hawking bats

This repository includes a pre-processed dataset and all necessary code for replicating the analyses in Material and code for an article: Lundberg et al. 2024 - Contrasting seasonal distribution patterns could reduce interspecific competition in two boreal aerial hawking bats

In the data file 'data.csv' there are following columns:

site: identifier for the recording site ie. the exact location where the device was located

x: E coordinate (ETRS89 / TM35FIN(E,N) - Finland - EPSG:3067)

y: N coordinate (ETRS89 / TM35FIN(E,N) - Finland - EPSG:3067)

coast: euclidean distance to the baltic coast (km)

north: euclidean distance to the 60 degN latitude (km). All sites are located north from 60 degN

year: the year (2019 or 2020)

rec.period: ordinal identifying the two-night-long recording period through the years 2019 and 2020. Each year has unique numbers.

rec.period2: ordinal identifying the two-night-long recording period per season. Corresponding dates in 2019 and 2020 share numbers.

night: starting date of a night

species: Eptesicus nilssonii ('EPTNIL') or Pipistrellus nathusii ('PIPNAT')

obs: number of 10 second time windows during which the species have been recorded

obs.min: number of 1 minute time windows during which the species have been recorded

This dataset is compiled from the recordings made during a citizen sience project, in which Finnish upper secondary students and bat researches deployed Audiomoths for two nights once in every two weeks during the summer and early autumn 2019 and 2020. Detailed description of the data collection is available at Lundberg et al. 2021-Next-generation ultrasonic recorders facilitate effective bat activity and distribution monitoring by citizen scientists (https://doi.org/10.1002/ecs2.3866).

During pre-processing of the data, raw bat observations have been cleaned to remove misidentifications and non-focal taxa. The remaining observations have been allocated into recording periods (two consecutive nights during which Audiomoth were scheduled to record) grouped by recording night, year, species and observation site.

The R script 'bbtcoast_script.R' includes all necessary code for fitting the models presented in the article. The user only needs to store 'bbtcoast_data.csv' in the working directory and run the script. The script was developed and tested with R version 4.2.0, tidyverse 1.3.1, spaMM_4.1.20, DHARMa 0.4.5.
