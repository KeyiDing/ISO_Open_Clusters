import numpy as np
import pandas as pd
from isochrones.mist import MIST_Isochrone
from isochrones.mist.bc import MISTBolometricCorrectionGrid
from isochrones import get_ichrone, SingleStarModel
from isochrones.priors import FlatPrior, PowerLawPrior
import matplotlib.pyplot as plt
from astropy.io import ascii as at
import os


def run_isochrones(row,name,base):
    #row is now the name of the row we are working on
    """
    here there are two alternatives:
    Keep working with the same data structure you had before (i. e. an Astropy Table)
    or use the pandas data structure. I will exemplify both

    First keeping the Tables structure:To do this you must transform the pandas
    structure into an astropy Table: (see packages imported above)
    """
    # row = tb.from_pandas(row)
    """
    If you choose this your code will keep the same structure but all the [0]]
    should be changed to 0:
    """
    # params_iso = {'parallax':(row['parallax'][0], row['parallax_error'][0])}

    """
    the second alternative is to keep the pandas dataframe structue. In this
    option be aware that getting a single value from a pandas dataframe is a bit different.
    (notice that I am using the .values[0] command because I also know there is
    only one row in my data!)
    example:
    """
    params_iso = {'parallax':(row['parallax'].values[0], row['parallax_error'].values[0])}

    """
    Please change the remaining of the code accoring to your choice. Memory wise
    pandas is generally better, but because we are dealing in this section with a single row it
    should not matter
    """
    #parallax cannot be negative
    if row['parallax'].values[0]<0:
        return

    """
    The second set of data are the actual stellar magnitudes (or brightness) in
    each pass band. This is done in two ways. First you tell isochrones which
    bands it should look for, what is being done in the bands array, and then
    you actually tell isochrones the values and uncertainties of each pass band. This
    is done through a dictionary, which I am calling mags_iso below.
    """
    # count how many Photometry are available
    count = 0

    bands = ['Gaia_G_DR2Rev']
    mags_iso = {'Gaia_G_DR2Rev':(row['gmag'].values[0],0.0003)}
    """
    In this example all the stars have Gaia G, 2MASS J, H, and K, and Wise W1, W2, and W3
    pass band data. However some stars data for some other passbands. The if(s) below
    append the bands and dictionary with data for other pass bands when available. See the
    data in the file appended.
    """
    # if GALEX FUV/NUV are available, then always include them
    if row['FUVmag'].isna().values[0]==False and row['e_FUVmag'].isna().values[0]==False:
        count += 1
        bands.append('GALEX_FUV')
        mags_iso['GALEX_FUV'] = (row['FUVmag'].values[0],row['e_FUVmag'].values[0])

    if row['NUVmag'].isna().values[0]==False and row['e_NUVmag'].isna().values[0]==False:
        count += 1
        bands.append('GALEX_NUV')
        mags_iso['GALEX_NUV'] = (row['NUVmag'].values[0],row['e_NUVmag'].values[0])
        
    # a flag variable that checks if SkyMapper U band exists
    #if only SkyMapper u-band data available, then use SkyMapper
    skymapper_flag = 0
    if row['u_psf'].isna().values[0]==False and row['e_u_psf'].isna().values[0]==False:
        skymapper_flag = 1
        bands.append('SkyMapper_u')
        mags_iso['SkyMapper_u'] = (row['u_psf'].values[0], row['e_u_psf'].values[0])

    if skymapper_flag == 1:
        count += 1
        if row['v_psf'].isna().values[0]==False and row['e_v_psf'].isna().values[0]==False:
            bands.append('SkyMapper_v')
            mags_iso['SkyMapper_v'] = (row['v_psf'].values[0], row['e_v_psf'].values[0])
        if row['g_psf'].isna().values[0]==False and row['e_g_psf'].isna().values[0]==False:
            bands.append('SkyMapper_g')
            mags_iso['SkyMapper_g'] = (row['g_psf'].values[0], row['e_g_psf'].values[0])
        if row['r_psf'].isna().values[0]==False and row['e_r_psf'].isna().values[0]==False:
            bands.append('SkyMapper_r')
            mags_iso['SkyMapper_r'] = (row['r_psf'].values[0], row['e_r_psf'].values[0])
        if row['i_psf'].isna().values[0]==False and row['e_i_psf'].isna().values[0]==False:
            bands.append('SkyMapper_i')
            mags_iso['SkyMapper_i'] = (row['i_psf'].values[0], row['e_i_psf'].values[0])
        if row['z_psf'].isna().values[0]==False and row['e_z_psf'].isna().values[0]==False:
            bands.append('SkyMapper_z')
            mags_iso['SkyMapper_z'] = (row['z_psf'].values[0], row['e_z_psf'].values[0])

    sdss_flag = 0
    #  if only SDSS u-band data available, then use SDSS
    if row['psfMag_uSDSS'].isna().values[0]==False and row['psfMagErr_uSDSS'].isna().values[0]==False:
        sdss_flag = 1
        bands.append('SDSS_u')
        mags_iso['SDSS_u'] = (row['psfMag_uSDSS'].values[0], row['psfMagErr_uSDSS'].values[0])

    if sdss_flag == 1:
        count += 1
        if row['psfMag_gSDSS'].isna().values[0]==False and row['psfMagErr_gSDSS'].isna().values[0]==False:
            bands.append('SDSS_g')
            mags_iso['SDSS_g'] = (row['psfMag_gSDSS'].values[0], row['psfMagErr_gSDSS'].values[0])
        if row['psfMag_rSDSS'].isna().values[0]==False and row['psfMagErr_rSDSS'].isna().values[0]==False:
            bands.append('SDSS_r')
            mags_iso['SDSS_r'] = (row['psfMag_rSDSS'].values[0], row['psfMagErr_rSDSS'].values[0])
        if row['psfMag_iSDSS'].isna().values[0]==False and row['psfMagErr_iSDSS'].isna().values[0]==False:
            bands.append('SDSS_i')
            mags_iso['SDSS_i'] = (row['psfMag_iSDSS'].values[0], row['psfMagErr_iSDSS'].values[0])
        if row['psfMag_zSDSS'].isna().values[0]==False and row['psfMagErr_zSDSS'].isna().values[0]==False:
            bands.append('SDSS_z')
            mags_iso['SDSS_z'] = (row['psfMag_zSDSS'].values[0], row['psfMagErr_zSDSS'].values[0])

    # if neither SkyMapper or SDSS u-band data is available, then use the survey with the most available bands
    if sdss_flag == 0 and skymapper_flag == 0:
        if row['v_psf'].isna().values[0]==False and row['e_v_psf'].isna().values[0]==False:
            count += 1
            bands.append('SkyMapper_v')
            mags_iso['SkyMapper_v'] = (row['v_psf'].values[0], row['e_v_psf'].values[0])
        if row['g_psf'].isna().values[0]==False and row['e_g_psf'].isna().values[0]==False:
            count += 1
            bands.append('SkyMapper_g')
            mags_iso['SkyMapper_g'] = (row['g_psf'].values[0], row['e_g_psf'].values[0])
        if row['r_psf'].isna().values[0]==False and row['e_r_psf'].isna().values[0]==False:
            count += 1
            bands.append('SkyMapper_r')
            mags_iso['SkyMapper_r'] = (row['r_psf'].values[0], row['e_r_psf'].values[0])
        if row['i_psf'].isna().values[0]==False and row['e_i_psf'].isna().values[0]==False:
            count += 1
            bands.append('SkyMapper_i')
            mags_iso['SkyMapper_i'] = (row['i_psf'].values[0], row['e_i_psf'].values[0])
        if row['z_psf'].isna().values[0]==False and row['e_z_psf'].isna().values[0]==False:
            count += 1
            bands.append('SkyMapper_z')
            mags_iso['SkyMapper_z'] = (row['z_psf'].values[0], row['e_z_psf'].values[0])

        if row['psfMag_gSDSS'].isna().values[0]==False and row['psfMagErr_gSDSS'].isna().values[0]==False:
            count += 1
            bands.append('SDSS_g')
            mags_iso['SDSS_g'] = (row['psfMag_gSDSS'].values[0], row['psfMagErr_gSDSS'].values[0])
        if row['psfMag_rSDSS'].isna().values[0]==False and row['psfMagErr_rSDSS'].isna().values[0]==False:
            count += 1
            bands.append('SDSS_r')
            mags_iso['SDSS_r'] = (row['psfMag_rSDSS'].values[0], row['psfMagErr_rSDSS'].values[0])
        if row['psfMag_iSDSS'].isna().values[0]==False and row['psfMagErr_iSDSS'].isna().values[0]==False:
            count += 1
            bands.append('SDSS_i')
            mags_iso['SDSS_i'] = (row['psfMag_iSDSS'].values[0], row['psfMagErr_iSDSS'].values[0])
        if row['psfMag_zSDSS'].isna().values[0]==False and row['psfMagErr_zSDSS'].isna().values[0]==False:
            count += 1
            bands.append('SDSS_z')
            mags_iso['SDSS_z'] = (row['psfMag_zSDSS'].values[0], row['psfMagErr_zSDSS'].values[0])

    if row['j_m'].isna().values[0] ==  False and row['j_msigcom'].isna().values[0] ==  False:
        count += 1
        bands.append('2MASS_J')
        mags_iso['2MASS_J'] = (row['j_m'].values[0],row['j_msigcom'].values[0])

    if row['h_m'].isna().values[0] ==  False and row['h_msigcom'].isna().values[0] ==  False:
        count += 1
        bands.append('2MASS_H')
        mags_iso['2MASS_H'] = (row['h_m'].values[0],row['h_msigcom'].values[0])

    if row['k_m'].isna().values[0] ==  False and row['k_msigcom'].isna().values[0] ==  False:
        count += 1
        bands.append('2MASS_Ks')
        mags_iso['2MASS_Ks'] = (row['k_m'].values[0],row['k_msigcom'].values[0])

    if row['w1mpro'].isna().values[0] ==  False and row['w1sigmpro'].isna().values[0] ==  False:
        count += 1
        bands.append('WISE_W1')
        mags_iso['WISE_W1'] = (row['w1mpro'].values[0],row['w1sigmpro'].values[0])
    if row['w2mpro'].isna().values[0] ==  False and row['w2sigmpro'].isna().values[0] ==  False:
        count += 1
        bands.append('WISE_W2')
        mags_iso['WISE_W2'] = (row['w2mpro'].values[0],row['w2sigmpro'].values[0])


    # if nothing other than Gaia DR2 G-band data are available by this step,
    # then add Gaia DR2 G_BP and G_RP data with appropriate uncertainties
    if count == 0:
        if row['bpmag'].isna().values[0] == False:
            bands.append('Gaia_BP_DR2Rev')
            mags_iso['Gaia_BP_DR2Rev'] = (row['bpmag'].values[0],row['dr2_bp_flux_over_error'].values[0])
        if row['rpmag'].isna().values[0] == False:
            bands.append('Gaia_RP_DR2Rev')
            mags_iso['Gaia_RP_DR2Rev'] = (row['rpmag'].values[0],row['dr2_rp_flux_over_error'].values[0])
    print(params_iso)
    print(mags_iso)
    """
    get_ichrone and SingleStarModel commands gets the grid of models from the
    available grids with the specific pass bands you requested (get_ichrone)
    and SingleStarModel create an initial stellar model based on your parameters
    and the grid of isochrones from get_ichrone.
    """
    mist = get_ichrone('mist', basic=False, bands=bands)
    model1 = SingleStarModel(mist, **params_iso, **mags_iso)
    """
    Below we give the code the necessary information about these stars that it
    needs to run the Bayesian statistics.
    """
    #You set the information about possible composition.
    # model1.set_prior(feh=FlatPrior((-2,0)), AV=PowerLawPrior(alpha=-2., bounds=(0.0001,1.4)))
    #specific for M67 (experimental, subect to change):
    model1._bounds['distance'] = (800, 900)
    model1._bounds['age'] = (np.log10(1.0e9), np.log10(13.721e9))
    model1._bounds['AV'] = (0.00, 0.40)
    model1.set_prior(AV=FlatPrior((0.08, 0.28)))

    """
    you bound your grid to certain distances. The lower and upper limits of
    distances. Do not worry about the names of the columns (r_lo_photogeo)
    because it is just telling how the distance was calculated in the
    particular study they came from.
    """
    # model1._bounds['distance'] = (mags['r_lo_photogeo'].values[0]]-20,mags['r_hi_photogeo'].values[0]]+20)
    """
    We could bound the age and mass, but we will not do that in this example.
    """
    #model1._bounds['age'] = (np.log10(1.0e9),np.log10(13.721e9))
    #model1._bounds['mass'] = (0.1,2)
    """
    below are limits to the extinction and composition. We only bound the
    extinction, for which we have data.
    """
    #model1._bounds['AV'] = (0,.3) #
    #model1._bounds['feh'] = (params['feh'][idx].values[0] - params['err_feh'][idx].values[0], params['feh'][idx].values[0] + params['err_feh'][idx].values[0])
    #Runs and saves the results.
    if os.path.isdir("./plots/{n}_plots/{id}".format(n=name,id=int(row['dr3_source_id'].values[0]))) == False:
        os.mkdir("./plots/{n}_plots/{id}".format(n=name,id=int(row['dr3_source_id'].values[0])))
    model1.fit(refit=True,n_live_points=1000,evidence_tolerance=0.5,max_iter=75000,basename=base)
    if len(model1.derived_samples)<8 or len(model1.derived_samples)<len(bands):
        model1 = SingleStarModel(mist, **params_iso, **mags_iso)
        model1.set_prior(feh=FlatPrior((-2, 0)), AV=PowerLawPrior(alpha=-2., bounds=(0.0001, 1.0)))
        model1.fit(refit=True, n_live_points=1000, evidence_tolerance=0.5, max_iter=150000)
        if len(model1.derived_samples)<8 or len(model1.derived_samples)<len(bands):
            model1.fit(refit=True, n_live_points=1000, evidence_tolerance=0.5, max_iter=225000)
            model1.set_prior(feh=FlatPrior((-2, 0)), AV=PowerLawPrior(alpha=-2., bounds=(0.0001, 1.0)))
            model1.fit(refit=True, n_live_points=1000, evidence_tolerance=0.5, max_iter=150000)
            if len(model1.derived_samples) < 8 or len(model1.derived_samples)<len(bands):
                return
    if os.path.isdir('./posteriors/{f}_posteriors'.format(f=name)) == False:
        os.mkdir('./posteriors/{f}_posteriors'.format(f=name))
    if os.path.isdir('./plots/{f}_plots'.format(f=name)) == False:
        os.mkdir('./plots/{f}_plots'.format(f=name))
    model1.derived_samples.to_csv("./posteriors/{f}_posteriors/{id}_take2.csv".format(f=name,id=int(row['dr3_source_id'].values[0])), index_label='index')
    plot1 = model1.corner_observed()
    plt.savefig("./plots/{f}_plots/{id1}/corner_{id2}.png".format(f=name,id1=int(row['dr3_source_id'].values[0]),id2=int(row['dr3_source_id'].values[0])))
    plot2 = model1.corner_physical()
    plt.savefig("./plots/{f}_plots/{id1}/physical_{id2}.png".format(f=name,id1=int(row['dr3_source_id'].values[0]),id2=int(row['dr3_source_id'].values[0])))
    """
    Here you might have one of the causes of your memory leak. matplotlib
    uses a lot of memory and you should close all open instances every time
    to avoid problems. This is true regardless of the code oyu are writing.
    """
    plt.clf()# cleans your plot1
    plt.close('all') #closes everything
    #Now when it leaves this function all local variables will be deleted!
    #but keep in mind that matplotlib might keep plots open in the memory so we need
    #to close them