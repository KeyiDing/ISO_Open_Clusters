# ISO_Open_Clusters
Repo for fitting stellar parameters of stars in open clusters using [Isochrones](https://github.com/timothydmorton/isochrones) analysis on [Advanced Research Computing at Hopkins](https://www.arch.jhu.edu/).

## Files and Directories Structure:
 - iso_input: photometric data for open clusters (in CSV format)
 - iso_code.py: Python script for running Isochrones
 - main_array.py: Python script for SLURM Job array
 - job_array: SLURM job submission file using array
 - main_mp.py: run isochrones through multiprocessing
 - job_mp: SLURM job submission file using multiprocessing
 - plots: Directory to save Isochrones analysis plots, the plots are organized in the format "cluster_name/Gaia_edr3_source_id"
 - posteriors: Directory to save Isochrones posteriors CSV data, organized in the format "cluster_name/Gaia_edr3_source_id_take2.csv"
 - nearby_cluster_av.nb: Mathematica file for numerically computing the extinctions of nearby clusters
 - records: a file for keeping records of Isochrones in case the job fails at some point. In the records, a CSV file is created for each
process named by the base_name, and records the loop index, dr3_source_id, and running time of that particular star. If the job for that process failed, start a new job by running Isochrones starting from the end of the index column

## Usage

### Method 1: SLURM Job Array Parallelism
1. in main_array.py, change the input photometry file path, name of the open cluster, and the length of the job array based on input.
2. In iso_code.py, specify extinctionV and parallax in section 1, and the Prior distribution in section 2.
3. In job_array, set "array" from 0 to length of job array -1. Also change the job name and output file directory.


### Method 2: Python Multiprocessing
1. In main_mp.py, change the input photometry file path, name of the open cluster, and the number of processes based on input.
2. In iso_code.py, specify extinctionV and parallax in section 1, and the Prior distribution in section 2.
3. In job_mp, set "ntasks-per-node" to the number of processes, and enter your email under "mail-user".

## Note on Python Multiprocessing
 - :warning: **Based on my experiment, in the multiprocessing method, the program won't terminate and no error code will be reported if one process files, and it is difficult to keep track of the error message for each process. Therefore the SLURM Job Array method is more encouraged to run the analysis in a large scale.** :warning:
 - Since the code does not require communication between processes, distributed memory is used to run the code on Rockfish, with each process having 16 GB of memory. Based on our past experiments, Isochrones usually takes less than 13 GB of memory, but feel free to increase it if needed.
 - In the job submission file, the number of CPU cores per node is specified through "ntasks-per-node", therefore it's essential to make sure that the number of processes in the Multiprocessing pool is less or equal to "ntasks-per-node". Theoretically the maximum number of processes per node is 48, but based on my experiment, a multiprocessing pool with approximately 40 processes will lead to IO mutex lock conflicts. I'm still investigating the potential cause.
 - If one of the process has an issue, the entire pool of processes will be stopped. In this case, manually specify the starting and ending point of the running loop for each process, and resubmit the job.
