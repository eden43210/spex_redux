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
import shutil
import re

data_col_header = {
    'UT_DATE' : ['DATE_OBS'],
    'UT_TIME' : ['TIME_OBS'],
    'OBSERVER' : ['OBSERVER'],
    'RA' : ['RA', 'RAJ2000', 'RA_ICRS'],
    'DEC' : ['DEC', 'DEJ2000', 'DE_ICRS'],
    'HOUR_ANGLE' : ['HA'],
    'POSITION_ANGLE' : ['POSANGLE'],
    'GRATING' : ['GRAT'],
    'SLIT_WIDTH' : ['SLIT'],
#    'program_info' :  I do not know what he wants with this yet
    'NAME' : ['OBJECT'],
    'EXP_TIME' : ['ITIME'],
    'AIRMASS' : ['AIRMASS'],
}
data_col_names = ['DATA_FILE', 'UT_DATE', 'UT_TIME', 'OBSERVER', 'RA', 'DEC', 'HOUR_ANGLE', 'POSITION_ANGLE', 'GRATING', 
                  'SLIT_WIDTH', 'PROGRAM_ID', 'NAME', 'OBSERVED_DESIGNATION', 'S_N', 'EXP_TIME', 'AIRMASS', 'REDUCTION_DATE', 'REDUCTION_PERSON']
data_non_head = ['DATA_FILE', 'PROGRAM_ID', 'OBSERVED_DESIGNATION',
                  'S_N', 'REDUCTION_DATE', 'REDUCTION_PERSON']
RA_list = {'2MASS': '2MASS_RAJ2000', 'SIMBAD' : 'SIMBAD_ra', 'GAIA' : 'GAIA_ra', 'WISE' : 'ALLWISE_RAJ2000', 'SDSS' : 'SDSS_RAdeg',
          'PANSTARRS': 'PANSTARRS_RAJ2000'}
DEC_list = {'2MASS': '2MASS_DEJ2000', 'SIMBAD' : 'SIMBAD_dec', 'GAIA' : 'GAIA_dec', 'WISE' : 'ALLWISE_DEJ2000', 'SDSS' : 'SDSS_DEdeg',
          'PANSTARRS': 'PANSTARRS_DEJ2000'}
cats_list1 = ['2MASS', 'SIMBAD']
cats_list2 = ['2MASS', 'SIMBAD', 'GAIA', 'ALLWISE', 'SDSS', 'II/349/ps1']
file_type_prefix = ['prism', 'sxd', 'lxd', 'com']
destination_folders = {'prism' : '/data/SpeX-reduced/calibrated/prism/',
                      'sxd' : '/data/SpeX-reduced/calibrated/sxd/',
                      'com' : '/data/SpeX-reduced/combspec/',
                      'csv' : '/data/SpeX-reduced/databases/',
                      'plot' : '/data/SpeX-reduced/review_plots/'}

#file_type_prefix = {'spex_prism': 'prism' , 'spex_sxd': 'sxd', 'spex_lxd': 'lxd', 'com_spec': 'notelluric'}

def search_files(path, keyword, file_type = '.fits'):
    #The purpose of this function is to search through a given folder for files that contain a chosen keyname.  IT also checks to see what type of file is available  
    assert os.path.isdir(path) == True, 'No such path exists.  Try a different path.'
    file_list = os.listdir(path)
    files = []
    for y in range(len(file_list)):
        if keyword in file_list[y] and file_type in file_list[y]:
            files += [path + file_list[y]]

    files = [f for f in files if 'tellspec' not in f]
    files = [f for f in files if 'modVega' not in f]
    return(files)

def copy_files(start_path, file_name, file_type, destination):
    #The purpose of this function is to identify a group of files that are to be copied and pasted into a specifed folder
    assert os.path.isdir(start_path) == True, 'No such path exists.  Try a different path.'
    #Start by searching through the chosen folder to return a list of all available files with that name.
    files = search_files(start_path, file_name, file_type)
    #Go through the list of files and copy them to their new directory.
    for a in range(len(files)):
        copy_name = files[a].split('proc/')[1]
        if file_name == 'com':
            com_date = [int(i) for i in files[a].split('/') if i.isdigit()]
            dir = os.path.join(destination, str(com_date[0]))
            if not os.path.exists(dir):
                os.mkdir(dir)
            shutil.copy2(files[a], dir + '/' + copy_name)
        else:
            shutil.copy2(files[a], destination + copy_name)

def data_sheet(path, file_type = 'spex_prism'):
    #The purpose of this function is to recieve a path to a folder filled with spex_prism.fits files and return a csv file containing
    #any relevant information from within the fits file.
    assert os.path.isdir(path) == True, 'No such path exists.  Try a different path.'
#    assert any(file_type_prefix[x] == file_type for x in range(len(file_type_prefix)))== True, 'Given file type does not exist'
    
    #Search for spex_prisim files
    files = search_files(path, file_type, '.fits')
    
    if len(files) == 0: print('Cannot find any {}*.fits data files in {}'.format(file_type, path))
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
            if x == 'DATA_FILE':
                spreadsheet.loc[i,x] = files[i].split('proc')[1]
            if x == 'PROGRAM_ID':
                spreadsheet.loc[i,x] = 'NaN'
            if x == 'OBSERVED_DESIGNATION':
                coords = splat.properCoordinates(header['RA']+ ' ' + header['DEC'])
                spreadsheet.loc[i,x] = splat.coordinateToDesignation(coords)
            if x == 'S_N':
                sp = splat.Spectrum(file = files[i])
                spreadsheet.loc[i,x] = sp.snr
            if x == 'REDUCTION_DATE':
                today = date.today()
                spreadsheet.loc[i,x] = today.strftime("%Y%m%d")
            if x == 'REDUCTION_PERSON':
                spreadsheet.loc[i,x] = getpass.getuser()
                
    for i in range(len(spreadsheet)):
        spreadsheet['UT_DATE'][i] = spreadsheet['UT_DATE'][i].replace('-', '')
        
    if file_type == 'sxd':
        spreadsheet.insert(1, 'MERGED_DATA_FILE', spreadsheet['DATA_FILE'])
        for z in range(len(spreadsheet['MERGED_DATA_FILE'])):
            if 'sxd_' in spreadsheet['MERGED_DATA_FILE'][z]:
                spreadsheet['MERGED_DATA_FILE'] = spreadsheet['MERGED_DATA_FILE'].replace(spreadsheet['MERGED_DATA_FILE'][z], '')
            if 'sxd-merged' in spreadsheet['DATA_FILE'][z]:
                spreadsheet['DATA_FILE'] = spreadsheet['DATA_FILE'].replace(spreadsheet['DATA_FILE'][z], '')
    
    if len(spreadsheet.index) == 0:
        print("No objects found")
        return(spreadsheet)
        
    date_taken = "/spectra_{UT}_{FT}.csv" #FT stands for File_Type and UT stands for Universal Time    
#    if file_type == 'com_spec':
#        date_taken = "/spectra_{UT}_notelluric.csv"
#    else:
#        date_taken = "/spectra_{UT}.csv"
#    print(file_type_prefix[file_type])
    date_taken = date_taken.format(UT = spreadsheet['UT_DATE'][0], FT = file_type)
    #Set index = False to get rid of the numbered column at the very beginning.
    spreadsheet.to_csv(path + date_taken, index=False)

    return(spreadsheet)


def info_tab(data, catalogs_given, up_coord_flag, file_type):
    #This function will have two separate ways of working based on the up_coord_flag.  If the flag is set to True, then the function will search
    # through the given catalogs list and determine the closest possible coordinates of the object in the OBSERVED_DESIGNATION dataframe.  Once
    # those coordinates are found, a dataframe is returned that gives the original OBSERVED_DESIGNATION, the most likely coordinates, and a column
    # that flags for issues.  If the flag is set to False, then the function will look for RA and DEC coordinates in the dataframe and use them to
    # search through a list of given catalogs to retrieve any relevent information.  The output should start with the same columns as it before it
    # went through the function followed by catalog data.  The second pass will also have to create a csv file.
    
    if len(data.index) == 0:
        return(data)
    
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
            #Note: On my laptop, the data says GAIA-EDR3 but on guac the data says GAIA-DR2.  Not sure why.
            for c in range(len(list(init_tab.columns))):
                if '-EDR3' in init_tab.columns[c]:
                    temp_col = init_tab.columns[c].replace("-EDR3", "")
                    init_tab.rename(columns = {init_tab.columns[c] : temp_col}, inplace = True)
                if '-DR2' in init_tab.columns[c]:
                    temp_col = init_tab.columns[c].replace("-DR2", "")
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
        date_taken = "/sources_{UT}_{FT}.csv" #FT stands for File_Type and UT stands for Universal Time           
        
#        if sys.argv[2] == 'com_spec':
#            print("YES")
#            date_taken = "/sources_{UT}_notelluric.csv"
#        else:
#            date_taken = "/sources_{UT}.csv"
        date_taken = date_taken.format(UT = data['UT_DATE'][0], FT = file_type)
        #Set index = False to get rid of the numbered column at the very beginning.
        rev_tab.to_csv(pathway + date_taken, index=False)
    
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
    nrow = np.ceil(len(sparr)/ncol)
    plt.figure(figsize=[ncol*6,nrow*4])
    sparr.reverse()
    for i,sp in enumerate(sparr):
        plt.subplot(nrow,ncol,i+1)
        plt.plot(sp.wave.value,sp.flux.value,'k-')
        plt.legend([sp.name],fontsize=14)
        plt.plot(sp.wave.value,sp.noise.value,'k-',alpha=0.3)
        plt.ylim([0,np.nanquantile(sp.flux.value,0.95)*1.2])
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


def process(folder,file_type = 'spex', smth=50):
    if os.path.isdir(folder) == False: raise ValueError('Cannot find folder {}'.format(folder))
        
    files = search_files(folder, file_type, '.fits')
        
        
    if file_type == 'sxd':
# plot individual SXD orders
        files_ind = []
        for x in files:
            if 'sxd_' in x:
                files_ind += [x]
        if len(files_ind) == 0: print('Cannot find any spex-sxd*.fits data files in {}'.format(folder))
        else:
            names = [f.split('_')[-2] for f in files_ind]
            dates = [(f.split('_')[-1]).split('.')[0] for f in files_ind]
            for i in range(len(files_ind)): 
                sps = readsxd(files_ind[i],name=names[i])
                plotmultispec(sps,output=folder+'/plot-orders_spex-sxd_{}_{}.pdf'.format(names[i],dates[i]))
                plt.savefig(folder+'/plot-orders_spex-sxd_{}_{}.png'.format(names[i],dates[i]))

# classify merged SXD files
        files_merged = []
        for x in files:
            if 'sxd-' in x:
                files_merged += [x]
        if len(files_merged) == 0: print('Cannot find any spex-sxd-merged*.fits data files in {}'.format(folder))
        else:
            names = [f.split('_')[-2] for f in files_merged]
            dates = [(f.split('_')[-1]).split('.')[0] for f in files_merged]
            for i,m in enumerate(files_merged): 
                sp = splat.Spectrum(file=m,instrument='SPEX-SXD',name=names[i])
                sp.smooth(smth)
                sp.trim([0.7,2.45])
                splat.classifyByStandard(sp,plot=True,method='kirkpatrick',telluric=True,output=folder+'/plot-classify_spex-sxd_{}_{}.pdf'.format(names[i],dates[i]))
                plt.savefig(folder+'/plot-classify_spex-sxd_{}_{}.png'.format(names[i],dates[i]))
                sp.plot(file=folder +'plot_spex-sxd_{}_{}.pdf'.format(names[i], dates[i]),telluric=True)
                plt.savefig(folder+'/plot_spex-sxd_{}_{}.png'.format(names[i],dates[i]))

    if file_type == 'prism':
# classify prism files
        if len(files) == 0: print('Cannot find any {}*.fits data files in {}'.format(file_type, folder))
        else:
            names = [f.split('_')[-2] for f in files]
            dates = [(f.split('_')[-1]).split('.')[0] for f in files]
            for i,m in enumerate(files): 
                sp = splat.Spectrum(file=m,instrument=file_type,name=names[i])
                splat.classifyByStandard(sp,plot=True,method='kirkpatrick',telluric=True,output=folder+'/plot-classify_spex-prism_{}_{}.pdf'.format(names[i],dates[i]))
                plt.savefig(folder+'/plot-classify_spex-prism_{}_{}.png'.format(names[i],dates[i]))
                sp.plot(file=folder + 'plot_spex-prism_{}_{}.pdf'.format(names[i],dates[i]),telluric=True)
                plt.savefig(folder+'/plot_spex-prism_{}_{}.png'.format(names[i],dates[i]))

    return

def check_files(files, keyword):
    #This is a revised, simplified version of the previous check_files function.  Just use re.search to determine if the file name has 'spex-prism' or 'spex-sxd'
    # in the name.  If it doesn't, use re.search to see if the individual words 'spex' and 'prism' or 'spex' and 'sxd' are in the file name.  If so, change it
    # to the proper notation, if not, print error and ask for user input.
    obj_list = list(files)
    
    for y in range(len(files)):
        fits_name = files[y].split("proc/")[1]
        file_dir = files[y].split("proc/")[0]
        split_file = fits_name.split("_")
        if re.search('spex-' + keyword, split_file[0]) is None:
            if split_file[0] == 'spex' and split_file[1] == keyword:
                prop_name = 'spex-' + keyword + '_' + split_file[2] + '_' + split_file[3]
                try:
                    os.rename(files[y], file_dir + 'proc/' + prop_name)
                    obj_list[y] = file_dir + 'proc/' + prop_name
                    continue
                except FileExistsError:
                    print("There was a naming issue with the following file:" + fits_name)
                    print("It appears that there is another file with the same file name in the folder. Please check that they are named properly.") 
                    print("For now, '-v2' will be added to the end of the object name.")
                    prop_name = 'spex-' + keyword + '_' + split_file[2] + '-v2' + '_' + split_file[3]
                    os.rename(files[y], file_dir + 'proc/' + prop_name)
                    obj_list[y] = file_dir + 'proc/' + prop_name
                    continue
            else:
                print("There was a naming issue with the following file:" + fits_name)
                print("Please carefully check the file for any mistakes and post the correct name below")
                prop_name = input("Correct name:")
                print("Attempting file check with the following file name: " + prop_name)
                os.rename(files[y], file_dir + 'proc/' + prop_name)
                obj_list[y] = file_dir + 'proc/' + prop_name
                continue
        if len(split_file) != 3:
            print("There was a naming issue with the following file:" + fits_name)
            print("Please carefully check the file for any mistakes and post the correct name below")
            prop_name = input("Correct name:")
            print("Attempting file check with the following file name: " + prop_name)
            os.rename(files[y], file_dir + 'proc/' + prop_name)
            obj_list[y] = file_dir + 'proc/' + prop_name
            continue
        
    for z in range(len(obj_list)):
        new_fits_name = obj_list[z].split("proc/")[1]
        new_file_dir = obj_list[z].split("proc/")[0]
        new_split_file = new_fits_name.split("_")
        corr_coords = check_coord(new_split_file[1])
        print("Reference file: " + new_fits_name)
        corr_date = check_date(new_split_file[2])
        prop_name = new_split_file[0] + '_' + corr_coords + '_' + corr_date
        os.rename(obj_list[z], new_file_dir + 'proc/' + prop_name)           

def check_coord(text):
    #This is a simpler version of the original check_coords function.  The sole purpose of this function is to check if the name has 9 characters and follows the 
    # following format: [0-9]{4}(-|+)[0-9]{4}.  If it does, add a J at the beginning.  Ignore anything else
    
    if type(re.search(r"[0-9]{4}[-+][0-9]{4}", text)) is re.Match and len(text) == 9:
        text = re.sub(r'([0-9]{4}[-+])', r'J\1', text, 1)
    if type(re.search(r"[0-9]{4}[-+][0-9]{4}", text)) is re.Match and len(text) == 12 and type(re.search(r"[-v][0-9]", text)) is re.Match:
        text = re.sub(r'([0-9]{4}[-+])', r'J\1', text, 1)
    return(text)

def check_date(text):
    num_text = re.findall("\d", text)
#    print(num_text)
    if len(num_text) == 6 and num_text[0] != 2:
        text = re.sub(r'([0-9]{6})', r'20\1', text)
        return(text)
    if len(num_text) != 8:
        print("Date incorrectly written for the following file: " + text)
        print("Please enter the correct date with the following format: YYYYMMDD")
        text_date = input("Correct coordinates:")
        text = text_date + '.fits'
    return(text)


# external function call

pathway = os.getcwd() + '/'
if __name__ == '__main__':
	if len(sys.argv) == 2:
		print(sys.argv[1])
		pathway = sys.argv[1]
	for x in list(file_type_prefix):
		if x == 'com':
			copy_files(pathway, 'com', '.fits', destination_folders['com'])
			continue
		if x == 'lxd':
			continue               #At the moment, the main focus of this code will revolve around prism and sxd data.  Don't waste time on lxd until later
#		print(pathway)
		file_check = search_files(pathway, x, '.fits')
		check_files(file_check, x)
		copy_files(pathway, x, '.fits', destination_folders[x])
		process(pathway, file_type = x)
		data_sheet(pathway, file_type = x)
#		print(data_sheet(path = pathway, file_type = x))  
		first = info_tab(data_sheet(path = pathway, file_type = x), catalogs_given = cats_list1, up_coord_flag = True, file_type = x)
		info_tab(data = first, catalogs_given = cats_list2, up_coord_flag = False, file_type = x)
	source = glob.glob(pathway + '/sources*.csv')
	data = []
	for x in range(len(source)):
		data += [pd.read_csv(source[x])]
	if len(data) == 2:
		output = data[1].append(data[0])
	if len(data) == 1:
		output = data[0]
	date_taken = source[0].split('proc')[1]   
	for y in range(len(source)):
		os.remove(source[y])
	output.to_csv(pathway + date_taken, index=False)
    #Once everything has been added to the local directory, the next phase is to copy the files and move them to a specified directory.
	copy_files(pathway, 'csv', '.csv', destination_folders['csv'])
	copy_files(pathway, 'plot', '.pdf', destination_folders['plot'])
	copy_files(pathway, 'plot', '.png', destination_folders['plot'])    
    
    
    
    
            