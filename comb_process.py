import splat
import splat.plot as splot
import splat.photometry as sphot
import splat.empirical as spem
from splat.empirical import estimateDistance, typeToColor, typeToMag
import splat.database as spdb

# other useful imports
import matplotlib.pyplot as plt
import numpy as np
import scipy
import pandas as pd
pd.set_option('display.max_columns', None)
import astropy.units as u
from astropy.io import fits
from astropy.utils.data import download_file

# Python standard library
import time
from datetime import date
import warnings
import math
import statistics
from itertools import combinations


# Third-party software
import numpy as np
import numpy.ma as ma

# Astropy
from astropy import coordinates as coord
from astropy import units as u
from astropy.table import Table, vstack, hstack, join
from astropy.coordinates import SkyCoord
from astropy.io import fits

# Astroquery. This tutorial requires 0.3.5 or greater.
import astroquery
from astroquery.simbad import Simbad
from astroquery.vo_conesearch import conf, conesearch, vos_catalog
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch

# Set up matplotlib
import matplotlib.pyplot as plt

import glob
import getpass
import os
import sys
os.chdir('/Users/cdann/Python/splat/resources/Spectra/SPEX-PRISM/')

data_col_header = {
    'UT_DATE' : ['DATE_OBS'],
    'UT_TIME' : ['TIME_OBS'],
    'OBSERVER' : ['OBSERVER'],
    'RA' : ['RA', 'RAJ2000', 'RA_ICRS'],
    'DEC' : ['DEC', 'DEJ2000', 'DE_ICRS'],
    'Hour_Angle' : ['HA'],
    'Position_Angle' : ['POSANGLE'],
    'Grating' : ['GRAT'],
    'Slit_Width' : ['SLIT'],
#    'program_info' :  I do not know what he wants with this yet
    'NAME' : ['OBJECT'],
    'EXP_time' : ['ITIME'],
    'Airmass' : ['AIRMASS'],
}
data_col_names = ['File_Name', 'UT_DATE', 'UT_TIME', 'OBSERVER', 'RA', 'DEC', 'Hour_Angle', 'Position_Angle', 'Grating', 
                  'Slit_Width', 'tell_list', 'program_info', 'NAME', 'OBSERVED_DESIGNATION', 'S_N', 'EXP_time', 'Airmass', 'RED_date', 'RED_person']
data_non_head = ['File_Name', 'tell_list', 'program_info', 'OBSERVED_DESIGNATION',
                  'S_N', 'RED_date', 'RED_person']
RA_list = {'2MASS': '2MASS_RAJ2000', 'SIMBAD' : 'SIMBAD_ra', 'GAIA' : 'GAIA_ra', 'WISE' : 'ALLWISE_RAJ2000', 'SDSS' : 'SDSS_RAdeg',
          'PANSTARRS': 'PANSTARRS_RAJ2000'}
DEC_list = {'2MASS': '2MASS_DEJ2000', 'SIMBAD' : 'SIMBAD_dec', 'GAIA' : 'GAIA_dec', 'WISE' : 'ALLWISE_DEJ2000', 'SDSS' : 'SDSS_DEdeg',
          'PANSTARRS': 'PANSTARRS_DEJ2000'}
cats_list1 = ['2MASS', 'SIMBAD']
cats_list2 = ['2MASS', 'SIMBAD', 'GAIA', 'ALLWISE', 'SDSS', 'II/349/ps1']

def data_sheet(path):
    #The purpose of this function is to recieve a path to a folder filled with spex_prism.fits files and return a csv file containing
    #any relevant information from within the fits file.
    assert os.path.isdir(path) == True, 'No such path exists.  Try a different path.'
    
    #Search for spex_prisim files
    files = glob.glob(path+'/spex*.fits')  #Remove _prism_ to also be able to include sxd_merged
    files = [f for f in files if 'tellspec' not in f]
    files = [f for f in files if 'modVega' not in f]
#    return(files)

    
    if len(files) == 0: print('Cannot find any spex*.fits data files in {}'.format(path))
    else:
        name_list = []
        for x in range(len(files)):
            if 'notelluric' in files[x]:
                tell_search = files[x].replace('_notelluric', '')
#                if any(y == tell_search for y in files):
#                    continue
                if any(y != tell_search for y in files):
                    name_list += [files[x]]
            else:
                name_list += [files[x]]
                
                
    spreadsheet = pd.DataFrame()
    for c in list(data_col_names): spreadsheet[c] = ['']*len(files)
    for i in range(len(files)):
        temp_list = []
        header = fits.open(files[i])[0].header
        for c in list(data_col_header.keys()): 
            if all(data_col_header[c][n] not in header for n in range(len(data_col_header[c]))):
                print("No column name found in header")
                print(data_col_header[c])
            for n in range(len(data_col_header[c])):
                if data_col_header[c][n] in header:
                    spreadsheet.loc[i,c] = header[data_col_header[c][n]]
            
        for x in data_non_head:
            if x == 'File_Name':
                spreadsheet.loc[i,x] = files[i].split('\\')[1]
            if x == 'tell_list':
                if 'notelluric' in files[i]:
                    spreadsheet.loc[i,x] = 'Combined Spectrum'
                else:
                    spreadsheet.loc[i,x] = 'Telluric Corrected'
            if x == 'program_info':
                spreadsheet.loc[i,x] = 'NaN'
            if x == 'OBSERVED_DESIGNATION':
                coords = splat.properCoordinates(header['RA']+ ' ' + header['DEC'])
                spreadsheet.loc[i,x] = splat.coordinateToDesignation(coords)
            if x == 'S_N':
                sp = splat.Spectrum(file = files[i])
                spreadsheet.loc[i,x] = sp.snr
            if x == 'RED_date':
                today = date.today()
                spreadsheet.loc[i,x] = today.strftime("%Y%m%d")
            if x == 'RED_person':
                spreadsheet.loc[i,x] = getpass.getuser()
                
    for i in range(len(spreadsheet)):
        spreadsheet['UT_DATE'][i] = spreadsheet['UT_DATE'][i].replace('-', '')
        
    
    date_taken = "/spectra_{UT}.csv"
    date_taken = date_taken.format(UT = spreadsheet['UT_DATE'][0].replace('-', ''))
    #Set index = False to get rid of the numbered column at the very beginning.
    spreadsheet.to_csv(path + date_taken, index=False)

    return(spreadsheet)

    #Search for SXD merged files
    files = glob.glob(path+'/spex-sxd-merged_*.fits')
    files = [f for f in files if 'tellspec' not in f]

    if len(files) == 0: print('Cannot find any spex-sxd-merged*.fits data files in {}'.format(path))
    else:
        return(files)
    
    return(path)


def info_tab(data, catalogs_given, up_coord_flag):
    #This function will have two separate ways of working based on the up_coord_flag.  If the flag is set to True, then the function will search
    # through the given catalogs list and determine the closest possible coordinates of the object in the OBSERVED_DESIGNATION dataframe.  Once
    # those coordinates are found, a dataframe is returned that gives the original OBSERVED_DESIGNATION, the most likely coordinates, and a column
    # that flags for issues.  If the flag is set to False, then the function will look for RA and DEC coordinates in the dataframe and use them to
    # search through a list of given catalogs to retrieve any relevent information.  The output should start with the same columns as it before it
    # went through the function followed by catalog data.  The second pass will also have to create a csv file.
    
    assert type(data)==pd.DataFrame, 'The data must be in the form of a pandas DataFrame'
    assert any(data.columns == 'OBSERVED_DESIGNATION')== True, 'The data must have a column that contains the OBSERVED_DESIGNATION'
    assert type(catalogs_given)== list, 'The catalogs must be in the form of a list'
    assert up_coord_flag == True or up_coord_flag == False, 'The up_coord_flag must be either True or False'
    
    #First, depending whether or not up_coord_flaf is True or False, the setup will be different.  Start both with a fresh dataframe.
    # If up_coord_falg is True, then only keep the OBSERVED_DESIGNATION column and create an empty RA and DEC column.  
    # If it is False, then keep the OBSERVED_DESIGNATION column along with RA and DEC columns.
    init_col = data[['OBSERVED_DESIGNATION', 'UT_DATE']]
    d = {'RA': [], 'DEC': [], 'flag': []}
    init_df = pd.DataFrame(data=d)
    for x in range(len(init_col)):
        if up_coord_flag == True:
            c = splat.properCoordinates(init_col['OBSERVED_DESIGNATION'][x])
            init_df.loc[x, 'RA'] = c.ra.value
            init_df.loc[x, 'DEC'] = c.dec.value
        if up_coord_flag == False:
            assert any(data.columns == 'RA')== True, 'The data must have a column that contains the RA'
            assert any(data.columns == 'DEC')== True, 'The data must have a column that contains the DEC'
            init_df = init_df.append({'RA': data['RA'][x], 'DEC': data['DEC'][x]}, ignore_index=True)
    init_df = pd.concat([init_col, init_df], axis = 1)
    init_df['flag'] = init_df['flag'].replace([init_df['flag'][0]], '')
#    return(init_df)

    #Now that the setup is finished, both dataframes will have a different purpose to serve.  First, start by searching the first catalog 
    #in the list at 30 arcseconds.  If up_coord_flag is True, then find the closest coordinates for each source and replace the RA and DEC
    # with the new coordinates.  If the angDist is greater than 15 arcseconds, add a flag in the in the flag column. If no objects are found
    # with 30 arcseconds, keep the old coordinates and move on to the next catalog and repeat.  If up_coord_flag is False, then run through the 
    # catalogs list and add all relevent information into the dataframe.
    rev_tab = pd.DataFrame()
    rev_tab = pd.concat([rev_tab, init_df], axis = 1)
    for x in range(len(catalogs_given)):
#        print(rev_tab)
        init_tab = spdb.queryXMatch(rev_tab, radius= 30 *u.arcsec, catalog=catalogs_given[x], drop_repeats = True)
        flag_val = catalogs_given[x]+' search flag'
        if catalogs_given[x] == 'GAIA':
            for c in range(len(list(init_tab.columns))):
                if '-EDR3' in init_tab.columns[c]:
                    temp_col = init_tab.columns[c].replace("-EDR3", "")
                    init_tab.rename(columns = {init_tab.columns[c] : temp_col}, inplace = True)                   
        
        for y in range(len(init_tab)):
            if catalogs_given[x] == 'II/349/ps1':
                if init_tab['II/349/PS1_angDist'][y] > 15:
                    init_tab['flag'][y] = init_tab['flag'][y] + flag_val
            else:
                if init_tab[catalogs_given[x]+'_angDist'][y] > 15:
                    init_tab['flag'][y] = init_tab['flag'][y] + flag_val      
        if up_coord_flag == True:
            temp_tab = init_tab.drop(['RA', 'DEC'], axis = 1)
            temp_tab.drop(temp_tab.columns.difference(['OBSERVED_DESIGNATION', 'flag', RA_list[catalogs_given[x]],
                                                               DEC_list[catalogs_given[x]], 'UT_DATE']), 1, inplace=True)
            temp_tab = temp_tab.rename(columns={RA_list[catalogs_given[x]]: "RA", DEC_list[catalogs_given[x]]: "DEC"})
            for z in range(len(temp_tab)):
                if str(init_tab[catalogs_given[x]+'_angDist'][z]) == 'nan':
                    temp_tab['RA'].iloc[z] = init_tab['RA'].iloc[z]
                    temp_tab['DEC'].iloc[z] = init_tab['DEC'].iloc[z]
            rev_tab = temp_tab
        if up_coord_flag == False:
            if catalogs_given[x] == 'II/349/ps1':
                for i in range(len(list(init_tab.columns))):
                    if 'II/349/PS1' in init_tab.columns[i]:
                        temp_col = init_tab.columns[i].replace("II/349/PS1", "PANSTARRS")
                        init_tab.rename(columns = {init_tab.columns[i] : temp_col}, inplace = True)
            rev_tab = pd.concat([rev_tab, init_tab], axis = 1)
            rev_tab = rev_tab.loc[:,~rev_tab.columns.duplicated()]
            
    if up_coord_flag == False:
        path = sys.argv[1]
        date_taken = "/sources_{UT}.csv"
        date_taken = date_taken.format(UT = data['UT_DATE'][0])
        #Set index = False to get rid of the numbered column at the very beginning.
        print(path + date_taken)
        rev_tab.to_csv(path + date_taken, index=False)
    
    return(rev_tab)

basefolder = '/Volumes/splat/data/spex/'


def readsxd(file,output='',name='',**kwargs):
    funit = u.erg/u.s/u.cm/u.cm/u.Angstrom
    wunit = u.micron
    xrngs = [[1.95,2.47],[1.43,1.81],[1.1,1.5],[0.92,1.2],[0.83,1.02],[0.735,0.93],[0.7,0.8]]
    with fits.open(file, **kwargs) as hdulist:
        header = hdulist[0].header
        meta = {'header': header}
        data = hdulist[0].data
    spec = []
    orders = header['ORDERS'].split(',')
    for i in range(len(data[:,0,0])):
        sp = splat.Spectrum(wave=data[i,0,:]*wunit,flux=data[i,1,:]*funit,noise=data[i,2,:]*funit,header=header,instrument='SPEX-SXD',name='{} order {}'.format(name,orders[i]))
        sp.trim(xrngs[i])
        spec.append(sp)
        
    if output=='multispec': return spec
    elif output=='1dspec':
        spc = spec[0]
        for s in spec[1:]: spc = splat.stitch(spc,s,scale=False)
        spc.name=name
        return spc
        
    else: return spec

def plotmultispec(sparr,output='',ncol=2):
    nrow = numpy.ceil(len(sparr)/ncol)
    plt.figure(figsize=[ncol*6,nrow*4])
    sparr.reverse()
    for i,sp in enumerate(sparr):
        plt.subplot(nrow,ncol,i+1)
        plt.plot(sp.wave.value,sp.flux.value,'k-')
        plt.legend([sp.name],fontsize=14)
        plt.plot(sp.wave.value,sp.noise.value,'k-',alpha=0.3)
        plt.ylim([0,numpy.nanquantile(sp.flux.value,0.95)*1.2])
        plt.xlabel('Wavelength ($\mu$m)',fontsize=14)
        plt.ylabel('F$_{\lambda}$'+' ({})'.format(sp.flux.unit),fontsize=14)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
    sparr.reverse()
    if output!='':
        plt.tight_layout()
        try:
            plt.savefig(output)
        except:
            print('Could not save to {}'.format(output))
    return


def process(folder,smth=50):
    if os.path.isdir(folder) == False: raise ValueError('Cannot find folder {}'.format(folder))

# plot individual SXD orders
    files = glob.glob(folder+'/spex-sxd_*.fits')
    files = [f for f in files if 'tellspec' not in f]
    files = [f for f in files if 'modVega' not in f]

    if len(files) == 0: print('Cannot find any spex-sxd*.fits data files in {}'.format(folder))
    else:
        names = [f.split('_')[-2] for f in files]
        dates = [(f.split('_')[-1]).split('.')[0] for f in files]
        for i in range(len(files)): 
            sps = readsxd(files[i],name=names[i])
            plotmultispec(sps,output=folder+'/plot_orders_{}_{}_sxd.pdf'.format(names[i],dates[i]))

# classify merged SXD files
    files = glob.glob(folder+'/spex-sxd-merged_*.fits')
    files = [f for f in files if 'tellspec' not in f]
    files = [f for f in files if 'modVega' not in f]

    if len(files) == 0: print('Cannot find any spex-sxd-merged*.fits data files in {}'.format(folder))
    else:
        names = [f.split('_')[-2] for f in files]
        dates = [(f.split('_')[-1]).split('.')[0] for f in files]
        for i,m in enumerate(files): 
            sp = splat.Spectrum(file=m,instrument='SPEX-SXD',name=names[i])
            sp.smooth(smth)
            sp.trim([0.7,2.45])
            splat.classifyByStandard(sp,plot=True,method='kirkpatrick',telluric=True,output=folder+'/plot_classify_{}_{}_sxd.pdf'.format(names[i],dates[i]))

# classify prism files
    files = glob.glob(folder+'/spex_prism_*.fits')
    files = [f for f in files if 'tellspec' not in f]
    files = [f for f in files if 'modVega' not in f]

    if len(files) == 0: print('Cannot find any spex_prism*.fits data files in {}'.format(folder))
    else:
        names = [f.split('_')[-2] for f in files]
        dates = [(f.split('_')[-1]).split('.')[0] for f in files]
        for i,m in enumerate(files): 
            sp = splat.Spectrum(file=m,instrument='SPEX-PRISM',name=names[i])
            splat.classifyByStandard(sp,plot=True,method='kirkpatrick',telluric=True,output=folder+'/plot_classify_{}_{}_prism.pdf'.format(names[i],dates[i]))

    return

# external function call
if __name__ == '__main__':
	if len(sys.argv) > 1: 
		process(sys.argv[1])
		data_sheet(sys.argv[1])
		first = info_tab(data_sheet(path = sys.argv[1]), catalogs_given = cats_list1, up_coord_flag = True)
		info_tab(data = first, catalogs_given = cats_list2, up_coord_flag = False)
        