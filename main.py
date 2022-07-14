import iso_code
import pandas as pd
import os.path
import sys
import time

def main():
    #name of the open cluster
    name = "Hyades_binary_rejected"
    #input photometry data
    data_input = pd.read_csv("iso_input/binary_rejected/Hyades_binary_rejected.csv", dtype={'dr2_source_id': int, 'dr3_source_id': int})

    base = "chain0"

    for i in range(200,len(data_input)):
        iso_code.run_isochrones(data_input.iloc[[i]],name,base)
        with open("./records/record.txt","w") as out:
            out.write(str(i))



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

