# ISO_Open_Clusters
Repo for fitting stellar parameters of stars in open clusters using isochrones analysis

## Files and Directories Structure:
 - iso_input: photometric data for open clusters (in CSV format)
 - main.py & iso_code.py: Python script for running Isochrones
 - plots: Directory to save Isochrones analysis plots, the plots are organized in the format "cluster_name/Gaia_edr3_source_id"
 - posteriors: Directory to save Isochrones posteriors CSV data, organized in the format "cluster_name/Gaia_edr3_source_id_take2.csv"
 - records: a file for keeping records of Isochrones in case the job fails at some point. In the records, a CSV file is created for each
process named by the base_name, and records the loop index, dr3_source_id, and running time of that particular star. If the job for that process failed, start a new job by running Isochrones starting from the end of the index column.
