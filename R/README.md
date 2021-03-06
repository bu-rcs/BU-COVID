# R scripts used in the pipeline

A few scripts in the pipeline are written in R. They are used to set-up input network files and to post-process the results.

* *pre-process_classroom_data.R* - Read input files describing classroom schedule and construct classroom networks
* *pre-process_house_data.R* - Read input files with housing information and construct building, floor, room, and houshold networks
* *plot_no_intervention_results.R* - Read files located in the results folders and display percentage of populations affected by the virus if no social distancing or other intervention is implemented.
* *plot_compare_interventions.R* - Compare results of various interventions
