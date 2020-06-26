#!/usr/bin/env python.
# Script to create metadata for export of SOCAT files.
# Adds < filename, hashsum, platformcode, filepath,level,export_date,expocode,startDate,endDate,nRows >  for each L2 file and linked L0 file in current directory

# linked L0s fetched from database generated by L0_links.py and Nuka_L0.py

# Maren K Karlsen 20200518

import os
import sys
import datetime
import sqlite3
import logging
import hashlib
import toml
import json
from scipy.io import netcdf

DB = 'SOCATexport.db'
logging.basicConfig(filename='console.log',format='%(asctime)s %(message)s', level=logging.DEBUG)

def main():
  now = datetime.datetime.now().strftime('%Y%m%d')
  c = create_connection() 
  L2_filenames=[] 
  
  for root,dirs,files in os.walk('.'):
    if 'L2' in root:
      for file in files:
        if file.endswith('.csv'):
          L2_filenames += [root + '/' + file]
        

  for L2_filepath in L2_filenames:
    L2_filename = L2_filepath.split('/')[-1]
    expocode = L2_filename.split('_')[0]
    print(expocode)
    platform_code = expocode[0:4]
    export_needed = True

    previous_entry = False
    try:
      c.execute("SELECT uploaded,hashsum FROM export WHERE filepath=?",[L2_filepath])
      previous_entry = c.fetchone()
    except Exception as e:
      raise Exception(f'Adding/Updating database failed: {L2_file}', exc_info=True)

    if not previous_entry:
      # Fetch information from L2 file
      L2_hashsum = get_hashsum(L2_filepath)
      with open(L2_filepath,'r',encoding="utf8",errors='ignore') as file:
        content = file.readlines()

      nRowsL2 = len(content)-1
      if nRowsL2 == 0: continue

      c.execute("INSERT INTO export \
        (filename, hashsum, platform, filepath,level,export_date,expocode,startDate,endDate,nRows)\
        VALUES (?,?,?,?,?,?,?,?,?,?)",
        (L2_filename, L2_hashsum, platform_code, L2_filepath,'L2',now, expocode,'startDate','endDate',nRowsL2))


    # Fetch information from linked L0s from L0-tbable
    c.execute("SELECT L0 FROM L0Links WHERE expocode=? ",[expocode])
    L0_links = c.fetchone()
    if L0_links:
      L0_files = L0_links[0].split(';')
    else:
      L0_files = []
    for L0_file in L0_files:
      print(L0_file)
      for root,dirs,files in os.walk('.'):
        for file in files:
          if file == L0_file:
            L0_filepath = root + '/' + file

      # check db
      previous_entry = False
      try:
        c.execute("SELECT uploaded,hashsum FROM export WHERE filepath=?",[L0_filepath])
        previous_entry = c.fetchone()
      except Exception as e:
        raise Exception(f'Adding/Updating database failed: {L0_file}', exc_info=True)


      if not previous_entry:
        # create metadata
        L0_hashsum = get_hashsum(L0_filepath)
        with open(L0_filepath,'r',encoding="utf8",errors='ignore') as file:
          content = file.readlines()

        nRowsL0 = len(content)-1
        if nRowsL0 == 0: continue

        ## L0 specific metadata 

        startDate = None
        endDate = None
        try:
          if 'SBE37' in L0_file: # SBE37SMP-ODO-RS232_03710886_2018_0108-0511_2.asc PALOMA
            t0 = '-'.join(content[1].split('\t')[1:3])
            t1 = '-'.join(content[-1].split('\t')[1:3])
            startDate = datetime.datetime.strptime(t0,'%m/%d/%Y-%H:%M:%S')
            endDate = datetime.datetime.strptime(t1,'%m/%d/%Y-%H:%M:%S')
          elif 'SD_datafile' in L0_file: #SD_datafile_20170309_110742CO2-0212-001.txt PALOMA
            t0 = '-'.join(content[5].split(';')[0:2]  )
            t1 = '-'.join(content[-1].split(';')[0:2])
            startDate = datetime.datetime.strptime(t0,'%Y-%m-%d-%H:%M:%S')
            endDate = datetime.datetime.strptime(t1,'%Y-%m-%d-%H:%M:%S')
          elif 'dat.txt' in L0_file:
            t0 = '-'.join(content[1].split('\t')[2:4])
            t1 = '-'.join(content[-1].split('\t')[2:4])
            startDate = datetime.datetime.strptime(t0,'%d/%m/%y-%H:%M:%S')
            endDate = datetime.datetime.strptime(t1,'%d/%m/%y-%H:%M:%S')
          elif '.cnv' in L0_file:
            t0 = content[42].split(' ')[3:7]
            duration = round(float(content[-1][34:45].strip())*3600)
            startDate = datetime.datetime.strptime('-'.join(t0),'%b-%d-%Y-%H:%M:%S')
            endDate = startDate + datetime.timedelta(0,duration)
          elif '.tsgqc' in L0_file:
            t0 = '-'.join(content[7].split(' ')[0:6])
            t1 = '-'.join(content[-1].split(' ')[0:6])
            startDate = datetime.datetime.strptime(t0,'%Y-%m-%d-%H-%M-%S')
            endDate = datetime.datetime.strptime(t1,'%Y-%m-%d-%H-%M-%S')
          elif '.nc' in L0_file:
            nc = netcdf.netcdf_file(L0_filepath,'r')
            t0 = nc.variables['DATE'][0,:].tostring().decode('utf-8')
            t1 = nc.variables['DATE'][-1,:].tostring().decode('utf-8')
            startDate = datetime.datetime.strptime(t0,'%Y%m%d%H%M%S')
            endDate = datetime.datetime.strptime(t1,'%Y%m%d%H%M%S')
          else:
            startDate = '-9'
            endDate = '-9'
        except:
          print(L0_file + '*** Missing start and/or end date ***')

        print(L0_file + ' ' + L0_hashsum)  
        c.execute("INSERT INTO export \
          (filename, hashsum, platform, filepath,level,export_date,expocode,startDate,endDate,nRows)\
          VALUES (?,?,?,?,?,?,?,?,?,?)",
          (L0_file, L0_hashsum, platform_code, L0_filepath,'L0',now, expocode,startDate,endDate,nRowsL0))


def create_connection():
  ''' creates connection and database if not already created '''
  conn = sqlite3.connect(DB, isolation_level=None)
  c = conn.cursor()
  c.execute(''' CREATE TABLE IF NOT EXISTS export (
              filename TEXT PRIMARY KEY,
              expocode TEXT, 
              hashsum TEXT NOT NULL UNIQUE,
              platform TEXT NOT NULL,
              filepath TEXT,
              level TEXT, 
              export_date TEXT,
              uploaded TEXT DEFAULT 'False',
              startDate TEXT,
              endDate TEXT,
              nRows TEXT
              )''')
  return c


def sql_investigate(filename, hashsum,level,platform):
  '''  Checks the sql database for identical filenames and hashsums
  returns 'exists', 'new', old_hashsum if 'update' and 'error' if failure. 
  '''
  c = create_connection(DB)
  status = {}
  try:
    c.execute("SELECT hashsum FROM export WHERE filename=? ",[filename])
    
    filename_exists = c.fetchone()
    if filename_exists: 
      if filename_exists[0] == hashsum:
        logging.info(f'{filename}: PREEXISTING entry')
        status = {'status':'EXISTS', 'info':'No action required'}
      else:
        logging.info(f'{filename}: UPDATE')
        status = ({'status':'UPDATE',
          'info':'Previously exported, updating entry',
          'old_hashsum':filename_exists[0]})
    else:
      logging.info(f'{filename}: NEW entry.')
      status = {'status':'NEW','info':'Adding new entry'}
  except Exception as e:
    logging.error(f'Checking database failed:  {filename} ', exc_info=True)
    status = {'status':'ERROR','info':e}

def get_hashsum(filename):
  ''' returns a 256 hashsum corresponding to input file. '''
  logging.debug(f'Generating hashsum for datafile {filename}')
  if filename.endswith('.nc'):
    with open(filename,'rb') as f: content = f.read()
    hashsum =  hashlib.sha256(content).hexdigest()
  else:
    with open(filename,'r',encoding="utf8",errors='ignore') as f: content = f.read()
    hashsum = hashlib.sha256(content.encode('utf-8')).hexdigest()

  return hashsum

if __name__ == '__main__':
  main()