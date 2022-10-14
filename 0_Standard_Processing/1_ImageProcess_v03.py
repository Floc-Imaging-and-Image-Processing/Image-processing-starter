#............................................................................
# Organizes files and then runs an ImageJ macro to obtain particle size information
#............................................................................

#user inputs
#............................................................................
# where = 'local' # use 'local' if you are running the file from within the folder you want to process
where = '0_Paths.csv' # manage paths with the 0_Paths.csv file if you don't run the file locally

ts = 0 # enter 1 if you want PSD time series output. Enter zero for a single PSD. If "1" then a dir_size must be specified
dir_size = 10  # number of images to include per distribution (will assume # per minute), ignored if ts = 0
image_type = '.jpg' # enter the file extention for your images

IJcode = 'ImageJ_code_diff_v02.txt'
# IJcode = "ImageJ_code_diff_v02_suppressoutput.txt"  
# IJcode = "ImageJ_code_diff_watershed_v02.txt"

rmsrc = 0 # remove source images variable. enter 0 for "no" and "1" for yes
subdirp = 0 # process all sub directories (subdirp = 1) or just the local one (subdirp = 0)
#............................................................................

# imports

import numpy as np
import os
from os.path import exists
from glob import glob
from shutil import copyfile
import shutil
import pandas as pd
import time
time1 = time.time()
from datetime import datetime

if where == 'local':
    master_path = os.getcwd()
    CodePath = master_path
else:
    paths = pd.read_csv(where)
    master_path = paths.iloc[0,1]
    CodePath = os.getcwd()

# for processing all folders in directory ------------------------------------------

if(subdirp == 1):

    folders = sorted(glob(master_path+"/*/"))

    for i in range(0,len(folders)):

        startTime = datetime.now()

        os.chdir(folders[i])
        print('working on',folders[i])

        path_main = folders[i]

        # get rid of previous processing files... will delete all folders in working directory

        if(exists('filemod_list.csv')==True):
            profolders = glob(path_main+"/*/")
            os.remove('filemod_list.csv')
            for i in range(0,len(profolders)):
                shutil.rmtree(profolders[i])

        # time series or no?

        if(ts == 0):
            dir_size = len(sorted(glob(path_main+'/*'+image_type)))

        # find files and re-organize them into directories with each directory containing the images to aggregate into a single PSD

        files_1 = [f for f in os.listdir(path_main) if os.path.isfile(os.path.join(path_main, f)) if (f.lower().endswith(image_type))]
        sorted_files = sorted(files_1)

        modtime_arr = np.zeros(len(files_1))
        for ig in np.arange(len(files_1)):
            file_path = path_main+'/'+sorted_files[ig]
            stats = os.stat(file_path)
            mod_time = stats.st_mtime
            modtime_arr[ig] = mod_time
        modtime_list =modtime_arr.tolist()

        src_list = []
        dst_list = []

        pad = len(str(dir_size)) # padding length for renaming and sorting of images within each directory

        for f in np.arange(np.floor(len(sorted_files)/dir_size)):
            print('Moving files to directory', int(f+1), 'of', int(np.floor(len(sorted_files)/dir_size)))
            path = path_main+'/'+str(int(f+1)).zfill(3)
            os.mkdir(path)
            path1 = path_main
            path2 = path
            for g in np.arange(dir_size):
                src = os.path.join(path1, sorted_files[int(f)*dir_size+g])
                dst = os.path.join(path2, (str(int(g+1)).zfill(pad)+image_type))
                copyfile(src, dst)
                src_list.append(src)
                dst_list.append(dst)

                if(rmsrc==1):
                    os.remove(src)

        # save dataframe of source, destination and modified time

        files_mod_df = pd.DataFrame(list(zip(sorted_files, modtime_list, src_list, dst_list)), columns =['file_Name', 'T_mod', 'Src', 'Dst'])
        file_modlist_csv = path_main+'/filemod_list.csv'
        files_mod_df.to_csv(file_modlist_csv,index=False)

        # run ImageJ

        str1 = 'java -jar '+CodePath+'/ij-1.52p.jar -batch '+ CodePath+ '/'+IJcode
        os.system(str1)

        # go back to the main folder
        print('run time:', datetime.now() - startTime, '(',len(files_1),' images in folder)')

        os.chdir(master_path)

    os.chdir(CodePath)

# for only the files in the current directory ------------------------------------------

if(subdirp == 0):

    startTime = datetime.now()

    os.chdir(master_path)

    path_main=master_path

    # get rid of previous processing files... will delete all folders in working directory

    if(exists('filemod_list.csv')==True):
        profolders = glob(path_main+"/*/")
        os.remove('filemod_list.csv')
        for i in range(0,len(profolders)):
            shutil.rmtree(profolders[i])

    # time series or no?

    if(ts == 0):
        dir_size = len(sorted(glob(path_main+'/*'+image_type)))

    # find files and re-organize them into directories with each directory containing the images to aggregate into a single PSD

    files_1 = [f for f in os.listdir(path_main) if os.path.isfile(os.path.join(path_main, f)) if (f.lower().endswith(image_type))]
    sorted_files = sorted(files_1)

    modtime_arr = np.zeros(len(files_1))
    for ig in np.arange(len(files_1)):
        file_path = path_main+'/'+sorted_files[ig]
        stats = os.stat(file_path)
        mod_time = stats.st_mtime
        modtime_arr[ig] = mod_time
    modtime_list =modtime_arr.tolist()

    src_list = []
    dst_list = []

    pad = len(str(dir_size)) # padding length for renaming and sorting of images within each directory

    for f in np.arange(np.floor(len(sorted_files)/dir_size)):
        print('Moving files to directory', int(f+1), 'of', int(np.floor(len(sorted_files)/dir_size)))
        path = path_main+'/'+str(int(f+1)).zfill(3)
        os.mkdir(path)
        path1 = path_main
        path2 = path
        for g in np.arange(dir_size):
            src = os.path.join(path1, sorted_files[int(f)*dir_size+g])
            dst = os.path.join(path2, (str(int(g+1)).zfill(pad)+image_type))
            copyfile(src, dst)
            src_list.append(src)
            dst_list.append(dst)

            if(rmsrc==1):
                os.remove(src)

    # save dataframe of source, destination and modified time

    files_mod_df = pd.DataFrame(list(zip(sorted_files, modtime_list, src_list, dst_list)), columns =['file_Name', 'T_mod', 'Src', 'Dst'])
    file_modlist_csv = path_main+'/filemod_list.csv'
    files_mod_df.to_csv(file_modlist_csv,index=False)

    # run ImageJ

    str1 = 'java -jar '+CodePath+'/ij-1.52p.jar -batch '+ CodePath+ '/'+IJcode
    os.system(str1)

    print('run time:', datetime.now() - startTime, '(',len(files_1),' images in folder)')
    os.chdir(CodePath)
