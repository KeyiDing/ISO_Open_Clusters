# ISO_Open_Clusters
Repo for fitting stellar parameters of stars in open clusters using [Isochrones](https://github.com/timothydmorton/isochrones) analysis on [Advanced Research Computing at Hopkins](https://www.arch.jhu.edu/).

## Files and Directories Structure:
 - iso_input: photometric data for open clusters (in CSV format)
 - main.py & iso_code.py: Python script for running Isochrones
 - main_mp.py: run isochrones through multiprocessing
 - plots: Directory to save Isochrones analysis plots, the plots are organized in the format "cluster_name/Gaia_edr3_source_id"
 - posteriors: Directory to save Isochrones posteriors CSV data, organized in the format "cluster_name/Gaia_edr3_source_id_take2.csv"
 - records: a file for keeping records of Isochrones in case the job fails at some point. In the records, a CSV file is created for each
process named by the base_name, and records the loop index, dr3_source_id, and running time of that particular star. If the job for that process failed, start a new job by running Isochrones starting from the end of the index column.

## Usage
1. In main_mp.py, change the input photometry file path, name of the open cluster, and the number of processes based on input.
2. In iso_code.py, specify extinctionV and parallax in section 1, and the Prior distribution in section 2.
3. In job, set "ntasks-per-node" to the number of processes, and enter your email under "mail-user".

## Note on Multiprocessing
 - Since the code does not require communication between processes, distributed memory is used to run the code on Rockfish, with each process having 16 GB of memory. Based on our past experiments, Isochrones usually takes less than 13 GB of memory, but feel free to increase it if needed.
 - In the job submission file, the number of CPU cores per node is specified through "ntasks-per-node", therefore it's essential to make sure that the number of processes in the Multiprocessing pool is less or equal to "ntasks-per-node". The maximum number of processes per node is 48.
 - If one of the process has an issue, the entire pool of processes would be stopped. In this case, manually specify the starting and ending point of the running loop for each process, and resubmit the job.
