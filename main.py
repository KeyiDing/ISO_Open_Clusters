import iso_code
import pandas as pd
import os.path
import sys
import time

def main():
    #name of the open cluster
    name = "M67_pht_fluxmag"
    #input photometry data
    data_input = pd.read_csv("iso_input/M67_pht_fluxmag.csv", dtype={'dr2_source_id': int, 'dr3_source_id': int})

    base = "chain0"

    for i in range(0,len(data_input)):
        iso_code.run_isochrones(data_input.iloc[[i]],name,base)
        with open("./records/record.txt","w") as out:
            out.write(str(i))



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

