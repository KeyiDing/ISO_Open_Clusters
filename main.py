import iso_code
import pandas as pd
import os.path
import sys
import time

def main():
    #name of the open cluster
    name = "Praesepe_plx_bound"
    #input photometry data
    data_input = pd.read_csv("iso_input/Praesepe_nonbinary.csv", dtype={'dr2_source_id': int, 'dr3_source_id': int})
    #specify the number of processes to run Isochrones in parallel, change it based on your need
    nprocess=10

    #clean up records folder before each new run
    filelist = [f for f in os.listdir("records") if f.endswith(".csv")]
    for f in filelist:
        os.remove(os.path.join("records", f))

    #read the index of current thread
    ind=4
    base = "chain{}".format(ind)
    record = pd.DataFrame(data={'index': [], 'source_id': [], 'time': []})

    for i in range(139,156):
        iso_code.run_isochrones(data_input.iloc[[i]],name,base)



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

