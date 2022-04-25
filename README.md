# AppModelling

This repository contains the underlying model code for the manuscript "The effect of notification window length on the epidemiological impact of COVID-19 mobile contact tracing applications", written in Matlab version 2021a.

To generate the underlying plots For Figures 2, and S3-S5, simply run 'NotificationWindow.m'. To generate the plots for Figure S1, uncomment relevant lines of code so that the base case individual is detected via a PCR test rather than an LFT.

The script 'NotificationWindowSimulation.m' runs the model as an explicit simulation (rather than calculating R via integrals), and is used to explore the impact of heterogeneity on results. This can be run to Figures S6 and S7, with values of 'HetMarker' designating what heterogeneities are simulated (HetMarker = 1 , heterogeneous contact rates, HetMarker = 2, heterogeneous infectious periods). To generate these figures, 'Rwindow', 'Reduction2' and 'Reduction5' from NotificationWindow.m  must be saved as 'Rwindowtrue', 'Reduction2true', and 'Reduction5true' respectively.

Test probability profiles included in this repository were obtained directly from [Hellewell et al. (2020)](https://cmmid.github.io/topics/covid19/pcr-positivity-over-time.html). 
