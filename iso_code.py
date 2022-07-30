import numpy as np
import pandas as pd
from isochrones.mist import MIST_Isochrone
from isochrones.mist.bc import MISTBolometricCorrectionGrid
from isochrones import get_ichrone, SingleStarModel
from isochrones.priors import FlatPrior, PowerLawPrior,GaussianPrior
import matplotlib.pyplot as plt
from astropy.io import ascii as at
from astropy.coordinates import SkyCoord
from astropy import units as u
import os
import time
import sys

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
        run_isochrones(df.iloc[[i]],name,base)
        end_time = time.time()
        record.loc[record.shape[0]] = [str(i), int(df.iloc[[i]]['dr3_source_id'].iloc[0]), end_time - start_time]
        record.to_csv("./records/{}.csv".format(base), index=False)


# row: the phtometry of current star we are working on
# name: name of the open cluster
# base: Multinest base name of the current process
def run_isochrones(row,name,base):
    # section 1: collect input physical and photometry data from the input dataframe

    #the first set of data are the known physical parameters of the star
    #in our code, we use the mean parallax of the open cluster as the parallax of each star
    #change the value of extinctionV and parallax based on the open clusters
    
    # Nearby Clusters
    # # get Celestial Coordinates
    # r=row["ra2000"].iloc[0]
    # d=row["dec2000"].iloc[0]
    # plx=row["parallax"].iloc[0]
    # c_icrs = SkyCoord(ra=r * u.degree, dec=d * u.degree, frame='icrs')
    # # convert coordinates into galactic, and calculate extinction
    # extinctionV = 0.7*(1-np.exp(-(1/plx)*np.sin(c_icrs.galactic.b)/0.125))*0.125/np.sin(c_icrs.galactic.b)
    # params_iso = {'parallax':(row["parallax"].iloc[0], row["parallax_error"].iloc[0])}
    
    # NGC0188
    params_iso = {'parallax':(0.5593,0.0011),'AV':(0.33,0.033)}

    # The second set of data are the actual stellar magnitudes (or brightness) in
    # each pass band. This is done in two ways. First you tell isochrones which
    # bands it should look for, what is being done in the bands array, and then
    # you actually tell isochrones the iloc and uncertainties of each pass band. This
    # is done through a dictionary, which I am calling mags_iso below.

    # count how many Photometry are available
    count = 0
    #always use Gaia_G
    bands = ['Gaia_G_DR2Rev']
    mags_iso = {'Gaia_G_DR2Rev':(row['gmag'].iloc[0],0.0003)}
    
    if row['FUVmag'].isna().iloc[0]==False and row['e_FUVmag'].isna().iloc[0]==False:
        count += 1
        bands.append('GALEX_FUV')
        mags_iso['GALEX_FUV'] = (row['FUVmag'].iloc[0],row['e_FUVmag'].iloc[0])
        
    if row['NUVmag'].isna().iloc[0]==False and row['e_NUVmag'].isna().iloc[0]==False:
        count += 1
        bands.append('GALEX_NUV')
        mags_iso['GALEX_NUV'] = (row['NUVmag'].iloc[0],row['e_NUVmag'].iloc[0])
    
        
    # a flag variable that checks if SkyMapper U band exists
    #if only SkyMapper u-band data available, then use SkyMapper
    skymapper_flag = 0
    if row['u_psf'].isna().iloc[0]==False and row['e_u_psf'].isna().iloc[0]==False:
        skymapper_flag = 1
        bands.append('SkyMapper_u')
        mags_iso['SkyMapper_u'] = (row['u_psf'].iloc[0], row['e_u_psf'].iloc[0])

    if skymapper_flag == 1:
        count += 1
        if row['v_psf'].isna().iloc[0]==False and row['e_v_psf'].isna().iloc[0]==False:
            bands.append('SkyMapper_v')
            mags_iso['SkyMapper_v'] = (row['v_psf'].iloc[0], row['e_v_psf'].iloc[0])
        if row['g_psf'].isna().iloc[0]==False and row['e_g_psf'].isna().iloc[0]==False:
            bands.append('SkyMapper_g')
            mags_iso['SkyMapper_g'] = (row['g_psf'].iloc[0], row['e_g_psf'].iloc[0])
        if row['r_psf'].isna().iloc[0]==False and row['e_r_psf'].isna().iloc[0]==False:
            bands.append('SkyMapper_r')
            mags_iso['SkyMapper_r'] = (row['r_psf'].iloc[0], row['e_r_psf'].iloc[0])
        if row['i_psf'].isna().iloc[0]==False and row['e_i_psf'].isna().iloc[0]==False:
            bands.append('SkyMapper_i')
            mags_iso['SkyMapper_i'] = (row['i_psf'].iloc[0], row['e_i_psf'].iloc[0])
        if row['z_psf'].isna().iloc[0]==False and row['e_z_psf'].isna().iloc[0]==False:
            bands.append('SkyMapper_z')
            mags_iso['SkyMapper_z'] = (row['z_psf'].iloc[0], row['e_z_psf'].iloc[0])

    sdss_flag = 0
    #  if only SDSS u-band data available, then use SDSS
    if row['psfMag_uSDSS'].isna().iloc[0]==False and row['psfMagErr_uSDSS'].isna().iloc[0]==False:
        sdss_flag = 1
        bands.append('SDSS_u')
        mags_iso['SDSS_u'] = (row['psfMag_uSDSS'].iloc[0], row['psfMagErr_uSDSS'].iloc[0])

    if sdss_flag == 1:
        count += 1
        if row['psfMag_gSDSS'].isna().iloc[0]==False and row['psfMagErr_gSDSS'].isna().iloc[0]==False:
            bands.append('SDSS_g')
            mags_iso['SDSS_g'] = (row['psfMag_gSDSS'].iloc[0], row['psfMagErr_gSDSS'].iloc[0])
        if row['psfMag_rSDSS'].isna().iloc[0]==False and row['psfMagErr_rSDSS'].isna().iloc[0]==False:
            bands.append('SDSS_r')
            mags_iso['SDSS_r'] = (row['psfMag_rSDSS'].iloc[0], row['psfMagErr_rSDSS'].iloc[0])
        if row['psfMag_iSDSS'].isna().iloc[0]==False and row['psfMagErr_iSDSS'].isna().iloc[0]==False:
            bands.append('SDSS_i')
            mags_iso['SDSS_i'] = (row['psfMag_iSDSS'].iloc[0], row['psfMagErr_iSDSS'].iloc[0])
        if row['psfMag_zSDSS'].isna().iloc[0]==False and row['psfMagErr_zSDSS'].isna().iloc[0]==False:
            bands.append('SDSS_z')
            mags_iso['SDSS_z'] = (row['psfMag_zSDSS'].iloc[0], row['psfMagErr_zSDSS'].iloc[0])

    # if neither SkyMapper or SDSS u-band data is available, then use the survey with the most available bands
    if sdss_flag == 0 and skymapper_flag == 0:
        if row['v_psf'].isna().iloc[0]==False and row['e_v_psf'].isna().iloc[0]==False:
            count += 1
            bands.append('SkyMapper_v')
            mags_iso['SkyMapper_v'] = (row['v_psf'].iloc[0], row['e_v_psf'].iloc[0])
        if row['g_psf'].isna().iloc[0]==False and row['e_g_psf'].isna().iloc[0]==False:
            count += 1
            bands.append('SkyMapper_g')
            mags_iso['SkyMapper_g'] = (row['g_psf'].iloc[0], row['e_g_psf'].iloc[0])
        if row['r_psf'].isna().iloc[0]==False and row['e_r_psf'].isna().iloc[0]==False:
            count += 1
            bands.append('SkyMapper_r')
            mags_iso['SkyMapper_r'] = (row['r_psf'].iloc[0], row['e_r_psf'].iloc[0])
        if row['i_psf'].isna().iloc[0]==False and row['e_i_psf'].isna().iloc[0]==False:
            count += 1
            bands.append('SkyMapper_i')
            mags_iso['SkyMapper_i'] = (row['i_psf'].iloc[0], row['e_i_psf'].iloc[0])
        if row['z_psf'].isna().iloc[0]==False and row['e_z_psf'].isna().iloc[0]==False:
            count += 1
            bands.append('SkyMapper_z')
            mags_iso['SkyMapper_z'] = (row['z_psf'].iloc[0], row['e_z_psf'].iloc[0])

        if row['psfMag_gSDSS'].isna().iloc[0]==False and row['psfMagErr_gSDSS'].isna().iloc[0]==False:
            count += 1
            bands.append('SDSS_g')
            mags_iso['SDSS_g'] = (row['psfMag_gSDSS'].iloc[0], row['psfMagErr_gSDSS'].iloc[0])
        if row['psfMag_rSDSS'].isna().iloc[0]==False and row['psfMagErr_rSDSS'].isna().iloc[0]==False:
            count += 1
            bands.append('SDSS_r')
            mags_iso['SDSS_r'] = (row['psfMag_rSDSS'].iloc[0], row['psfMagErr_rSDSS'].iloc[0])
        if row['psfMag_iSDSS'].isna().iloc[0]==False and row['psfMagErr_iSDSS'].isna().iloc[0]==False:
            count += 1
            bands.append('SDSS_i')
            mags_iso['SDSS_i'] = (row['psfMag_iSDSS'].iloc[0], row['psfMagErr_iSDSS'].iloc[0])
        if row['psfMag_zSDSS'].isna().iloc[0]==False and row['psfMagErr_zSDSS'].isna().iloc[0]==False:
            count += 1
            bands.append('SDSS_z')
            mags_iso['SDSS_z'] = (row['psfMag_zSDSS'].iloc[0], row['psfMagErr_zSDSS'].iloc[0])

    # always use 2Mass if available
    if row['j_m'].isna().iloc[0] ==  False and row['j_msigcom'].isna().iloc[0] ==  False:
        count += 1
        bands.append('2MASS_J')
        mags_iso['2MASS_J'] = (row['j_m'].iloc[0],row['j_msigcom'].iloc[0])
    if row['h_m'].isna().iloc[0] ==  False and row['h_msigcom'].isna().iloc[0] ==  False:
        count += 1
        bands.append('2MASS_H')
        mags_iso['2MASS_H'] = (row['h_m'].iloc[0],row['h_msigcom'].iloc[0])
    if row['k_m'].isna().iloc[0] ==  False and row['k_msigcom'].isna().iloc[0] ==  False:
        count += 1
        bands.append('2MASS_Ks')
        mags_iso['2MASS_Ks'] = (row['k_m'].iloc[0],row['k_msigcom'].iloc[0])

    # always use Wise if available
    if row['w1mpro'].isna().iloc[0] ==  False and row['w1sigmpro'].isna().iloc[0] ==  False:
        count += 1
        bands.append('WISE_W1')
        mags_iso['WISE_W1'] = (row['w1mpro'].iloc[0],row['w1sigmpro'].iloc[0])
    if row['w2mpro'].isna().iloc[0] ==  False and row['w2sigmpro'].isna().iloc[0] ==  False:
        count += 1
        bands.append('WISE_W2')
        mags_iso['WISE_W2'] = (row['w2mpro'].iloc[0],row['w2sigmpro'].iloc[0])
        

    # always use Pan-Starrs if available
    if row['gMeanPSFMag'].isna().iloc[0] == False and row['gMeanPSFMagErr'].isna().iloc[0] == False and row['gMeanPSFMag'].iloc[0]>14.5:
        count += 1
        bands.append("PS_g")
        mags_iso["PS_g"] = (row["gMeanPSFMag"].iloc[0],row["gMeanPSFMagErr"].iloc[0])
    if row['rMeanPSFMag'].isna().iloc[0] == False and row['rMeanPSFMagErr'].isna().iloc[0] == False and row['rMeanPSFMag'].iloc[0]>15:
        count += 1
        bands.append("PS_r")
        mags_iso["PS_r"] = (row["rMeanPSFMag"].iloc[0],row["rMeanPSFMagErr"].iloc[0])
    if row['iMeanPSFMag'].isna().iloc[0] == False and row['iMeanPSFMagErr'].isna().iloc[0] == False and row['iMeanPSFMag'].iloc[0]>15:
        count += 1
        bands.append("PS_i")
        mags_iso["PS_i"] = (row["iMeanPSFMag"].iloc[0],row["iMeanPSFMagErr"].iloc[0])
    if row['zMeanPSFMag'].isna().iloc[0] == False and row['zMeanPSFMagErr'].isna().iloc[0] == False and row['zMeanPSFMag'].iloc[0]>14:
        count += 1
        bands.append("PS_z")
        mags_iso["PS_z"] = (row["zMeanPSFMag"].iloc[0],row["zMeanPSFMagErr"].iloc[0])
    if row['yMeanPSFMag'].isna().iloc[0] == False and row['yMeanPSFMagErr'].isna().iloc[0] == False and row['yMeanPSFMag'].iloc[0]>13:
        count += 1
        bands.append("PS_y")
        mags_iso["PS_y"] = (row["yMeanPSFMag"].iloc[0],row["yMeanPSFMagErr"].iloc[0])
        

    # if nothing other than Gaia DR2 G-band data are available by this step,
    # then add Gaia DR2 G_BP and G_RP data with appropriate uncertainties
    if count == 0:
        err = 0
        gmag = row["gmag"].iloc[0]
        if gmag < 13:
            err = 0.002
        elif gmag < 17:
            err = 0.01
        else:
            err = 0.2
        if row['bpmag'].isna().iloc[0] == False:
            bands.append('Gaia_BP_DR2Rev')
            mags_iso['Gaia_BP_DR2Rev'] = (row['bpmag'].iloc[0],err)
        if row['rpmag'].isna().iloc[0] == False:
            bands.append('Gaia_RP_DR2Rev')
            mags_iso['Gaia_RP_DR2Rev'] = (row['rpmag'].iloc[0],err)
    print(params_iso)
    print(mags_iso)


    # Section 2: get the grid of models and specify prior distribution

    # get_ichrone and SingleStarModel commands gets the grid of models from the
    # available grids with the specific pass bands you requested (get_ichrone)
    # and SingleStarModel create an initial stellar model based on your parameters
    # and the grid of isochrones from get_ichrone.
    mist = get_ichrone('mist', basic=False, bands=bands)
    model1 = SingleStarModel(mist, **params_iso, **mags_iso)


    # specify Prior distribution and bound to do Bayesian Statistics
    
    # Nearby clusters
    # model1.set_prior(AV=FlatPrior((0, 2*extinctionV)))
    # model1._bounds['AV'] = (0, 2*extinctionV)
    # # model1.set_prior(AV=GaussianPrior(mean=extinctionV, sigma=0.1 * extinctionV))
    # MeanDistance = 1000 / params_iso['parallax'][0]
    # HighDistance = 1000 / (params_iso['parallax'][0]-3*params_iso['parallax'][1])
    # LowDistance = 1000 / (params_iso['parallax'][0]+3*params_iso['parallax'][1])
    # model1._bounds["distance"] = (LowDistance, HighDistance)
    # model1._bounds['age'] = (8, np.log10(13.721e9))
    # model1._bounds['feh'] = (-1,0.5)

    # M67
    # LowerExtinction = np.max([0,extinctionV-3*0.10*extinctionV])
    # UpperExtinction = extinctionV+3*0.10*extinctionV
    # model1.set_prior(AV=FlatPrior((LowerExtinction, UpperExtinction)))
    # model1._bounds['AV'] = (LowerExtinction, UpperExtinction)
    # model1._bounds['distance'] = (813,873)
    # model1._bounds['age'] = (np.log10(1.0e9),np.log10(13.721e9))
    
    # NGC0188
    extinctionV = 0.33
    LowerExtinction = np.max([0,extinctionV-3*0.10*extinctionV])
    UpperExtinction = extinctionV+3*0.10*extinctionV
    model1.set_prior(AV=FlatPrior((LowerExtinction, UpperExtinction)))
    model1._bounds['AV'] = (LowerExtinction, UpperExtinction)
    MeanDistance = 1000 / params_iso['parallax'][0]
    HighDistance = 1000 / (params_iso['parallax'][0]-3*params_iso['parallax'][1])
    LowDistance = 1000 / (params_iso['parallax'][0]+3*params_iso['parallax'][1])
    model1._bounds["distance"] = (LowDistance, HighDistance)
    model1._bounds['age'] = (np.log10(1.0e9),np.log10(13.721e9))
    model1._bounds['feh'] = (-1,0.5)

    # Section 3: runs and saves the results.
    # create plots and posteriors folder
    if os.path.isdir('./plots/{f}_plots'.format(f=name)) == False:
        os.mkdir('./plots/{f}_plots'.format(f=name))
    if os.path.isdir("./plots/{n}_plots/{id}".format(n=name,id=int(row['dr3_source_id'].iloc[0]))) == False:
        os.mkdir("./plots/{n}_plots/{id}".format(n=name,id=int(row['dr3_source_id'].iloc[0])))
    if os.path.isdir('./posteriors/{f}_posteriors'.format(f=name)) == False:
        os.mkdir('./posteriors/{f}_posteriors'.format(f=name))
    

    # fit parameters to the model
    model1.fit(refit=True,n_live_points=1000,evidence_tolerance=0.5,max_iter=100000,basename=base)
    if len(model1.derived_samples)<8 or len(model1.derived_samples)<len(bands):
        return

    # save the posterior, physical plot, and observed plot
    model1.derived_samples.to_csv("./posteriors/{f}_posteriors/{id}_take2.csv".format(f=name,id=int(row['dr3_source_id'].iloc[0])), index_label='index')
    plot1 = model1.corner_observed()
    plt.savefig("./plots/{f}_plots/{id1}/corner_{id2}.png".format(f=name,id1=int(row['dr3_source_id'].iloc[0]),id2=int(row['dr3_source_id'].iloc[0])))
    plot2 = model1.corner_physical()
    plt.savefig("./plots/{f}_plots/{id1}/physical_{id2}.png".format(f=name,id1=int(row['dr3_source_id'].iloc[0]),id2=int(row['dr3_source_id'].iloc[0])))

    # close the plots to prevent memory leak
    plt.clf()
    plt.close('all')



# 1. update extinctionV based on https://stilism.obspm.fr/ for alphaPer (likelihood and prior), and NearBy clusters
#  - update AV upper bound for not converging stars
# 2. For M67, replace SDSS with PAN-STARRS to check for the effect on UV band
# 3. for Blanco1, compare photometry use for the "structure"
# 4. develop a automatic way to identify binaries and variable stars based on Simbad
# 5. Prepare for report and paper: Latex for AAS, compare with spectrum (Galah, Lamost) and compare (dispersion...)
# LAMOST: http://www.lamost.org/dr7/v2.0/catalogue