import iso_code
import pandas as pd
import os.path
import sys
import time

def main():
    #name of the open cluster
    name = "Praesepe"
    #input photometry data
    data_input = pd.read_csv("iso_input/Praesepe_nonbinary.csv", dtype={'dr2_source_id': int, 'dr3_source_id': int})
    #specify the number of processes to run Isochrones in parallel, change it based on your need
    nprocess=20

    #clean up records folder before each new run
    filelist = [f for f in os.listdir("records") if f.endswith(".csv")]
    for f in filelist:
        os.remove(os.path.join("records", f))

    #read the index of current thread
    ind = int(sys.argv[1])
    base = "chain{}".format(ind)

    #run isochrones
    iso_code.iso_process(nprocess,data_input,ind,name,base)



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

