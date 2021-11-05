# AppModelling

This repository contains the underlying model code for the manuscript "The effect of notification window length on the epidemiologicalimpact of COVID-19 mobile contact tracing applications", written in Matlab version 2021a.

To generate the underlying plots For Figures 2, S2, and S3, simply run 'NotificationWindow.m'. To generate the plots for Figure S1, uncomment relevant lines of code so that the base case individual is detected via a PCR test rather than an LFT.

The script 'NotificationWindowSimulation.m' runs the model as an explicit simulation (rather than calculating R via integrals). This can be run to generate analogous figures to Figures 2B and 2C.

Test probability profiles included in this repository were obtained directly from [Hellewell et al. (2020)](https://cmmid.github.io/topics/covid19/pcr-positivity-over-time.html). 
