import iso_code
import pandas as pd
import os.path
import sys
from multiprocessing import Pool
from multiprocessing import Process

def iso_process(tot,df,ind,name,base):
    start = int(ind*(len(df)/tot))
    end = int((ind+1)*(len(df)/tot))
    for i in range(start,end):
        iso_code.run_isochrones(i,df.iloc[[i]], name, base)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    name = "alphaPer"
    data_input = pd.read_csv("iso_input/alphaPer_nonbinary.csv", dtype={'dr2_source_id': int, 'dr3_source_id': int,'j_msigcom':float,'h_msigcom':float,'k_msigcom':float})
    tot=5
    with Pool(processes=tot) as pool:
        for i in range(tot):
            base = "chain"+str(i)
            pool.apply_async(iso_process, args=(tot,data_input,i,name,base,))
        pool.close()
        pool.join()
