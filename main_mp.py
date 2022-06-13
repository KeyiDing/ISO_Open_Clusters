import iso_code
import pandas as pd
import os.path
import sys
from multiprocessing import Pool
import time

def iso_process(tot,df,ind,name,base):
    #evenly divide the whole input photometry dataframe into n sub-dataframes
    start = int(ind*(len(df)/tot))
    end = int((ind+1)*(len(df)/tot))

    #create a "record" dataframe to record running time for each process
    #in each entry of the record dataframe contains the index of a star in the input photometry file, the source id, and the isochrones running time for that star
    record = pd.DataFrame(data={'index': [], 'source_id': [], 'time': []})

    #run isochrones and record corrosponding info
    for i in range(start,end):
        start_time = time.time()
        iso_code.run_isochrones(df.iloc[[i]],name,base)
        end_time = time.time()
        print(i)
        record.loc[record.shape[0]] = [str(i), int(row['dr3_source_id'].iloc[0]), end_time - start_time]
        record.to_csv("./records/{}.csv".format(base), index=False)

if __name__ == '__main__':
    #name of the open cluster
    name = "Blanco1"
    #input photometry data
    data_input = pd.read_csv("iso_input/Blanco1_nonbinary.csv", dtype={'dr2_source_id': int, 'dr3_source_id': int})
    #specify the number of processes to run Isochrones in parallel, change it based on your need
    nprocess=3

    #clean up records folder before each new run
    filelist = [f for f in os.listdir("records") if f.endswith(".csv")]
    for f in filelist:
        os.remove(os.path.join("records", f))

    #create a multiprocessing pool of size nprocess
    with Pool(processes=nprocess) as pool:
        for i in range(nprocess):
            #format Multinest base name for each process
            base = "chain"+str(i)
            #apply each process in the pool
            #since our code doesn't require communication between processes, and there's no further steps after each process finishes,
            #asynchronization is not really needed, but is generally a good practice to avoid race
            pool.apply_async(iso_process, args=(nprocess,data_input,i,name,base,))
        pool.close()
        pool.join()
