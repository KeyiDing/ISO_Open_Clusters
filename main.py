import iso_code
import pandas as pd
import os.path
import sys

def main():
    name="alphaPer"
    #I am reding the input data as a pandas dataframe because it handles large files much better
    #and saves memory. Also, it is faster. I set up the data type for source id to read the entire number as integer
    #and avoid truncating it
    data_input = pd.read_csv("iso_input/alphaPer_nonbinary.csv", dtype={'dr2_source_id': int, 'dr3_source_id': int})
    # data_input = data_input[(data_input['Cluster'] == 'NGC2682 ')]
    name = name.strip()
    #Although standard pratice is to vectorize usage, in our case we do need to use a for loop.
    #Just a for loop to run it for all the stars in my file list
    len = data_input.shape[0]
    base = "chain1"

    for i in range(209,len):
        #The iloc function will get the row we are interested in running and keep the data structure
        iso_code.run_isochrones(data_input.iloc[[i]],name,base)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

