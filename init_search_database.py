#!/usr/bin/env python
import os
import sqlite3 as lite
try:
    DB_FILE = os.environ['SEARCH_DB']
except:
    print("environmental variable CMD_BS_DB_DEF_FILE must be defined")
    
con = lite.connect(DB_FILE)
with con:
    #PulsarSearch Table: stores all major values for the mwa search pipeline for each pointing
    #Rownum: ID to reference
    #Comment: What the pipeline's purpose (eg testing or main search centered at x y)
    #SampleNum: Number of time samples in fits files used
    #*Proc: total prossing wall time in seconds
    #*Errors: total errors for each section
    #*DS: Number of jobs that didn't start
    #*DE: number of jobs that didn't end
    #*JobComp: Number of jobs that completed succesfully
    #*JobExp: Expected number of jobs to be completed
    #CandTotal: total cands
    #CandOverNoise: total cands over noise (normaly 3 sigma but I may change it)
    #CandDect: cands that were confirmed as pulsars (known or hopefully new)
    cur = con.cursor()
    cur.execute("""CREATE TABLE PulsarSearch(
                Rownum INTEGER PRIMARY KEY AUTOINCREMENT, 
                Started date, 
                Ended date, 
                Obsid INT, 
                Pointing TEXT,
                UserID TEXT,
                Version TEXT,
                Comment TEXT, 
                SampleNum INT,
                TotalProc FLOAT, 
                TotalErrors INT, 
                TotalDS INT,
                TotalDE INT,
                TotalJobComp INT,
                TotalJobExp INT,
                BeamformProc FLOAT, 
                BeamformErrors INT,  
                BeamformDS INT,
                BeamformDE INT,
                BeamformJobComp INT,
                BeamformJobExp INT,
                PrepdataProc FLOAT, 
                PrepdataErrors INT,  
                PrepdataDS INT,
                PrepdataDE INT,
                PrepdataJobComp INT,
                PrepdataJobExp INT,
                FFTProc FLOAT, 
                FFTErrors INT,  
                FFTDS INT,
                FFTDE INT,
                FFTJobComp INT,
                FFTJobExp INT,
                AccelProc FLOAT, 
                AccelErrors INT,  
                AccelDS INT,
                AccelDE INT,
                AccelJobComp INT,
                AccelJobExp INT,
                FoldProc FLOAT, 
                FoldErrors INT,  
                FoldDS INT,
                FoldDE INT,
                FoldJobComp INT,
                FoldJobExp INT,
                CandTotal INT, 
                CandOverNoise INT, 
                CandDect INT)""")
    
    #A table for each job type
    
    #Rownum: ID to reference
    #BSID: PulsarSearch id Rownum reference
    #Command: the code being run
    #Arguments: arguments put into the code
    #Exit: error code
    #CPUs: Number of cpus being used to be multiplied by processing time (end-start)
    #DMFileInt: Dm file reference eg 1: DM_002_004
    
    cur.execute("""CREATE TABLE Beamform(
                Rownum INT NOT NULL, 
                AttemptNum INT NOT NULL,
                BSID INT NOT NULL, 
                Command TEXT, 
                Arguments TEXT, 
                Started date, 
                Ended date, 
                Proc FLOAT,
                ExpProc FLOAT,
                Exit INT, 
                CPUs INT, 
                FOREIGN KEY (BSID) REFERENCES PulsarSearch(Rownum),
                primary key (Rownum, AttemptNum, BSID))""")
    
    cur.execute("""CREATE TABLE Prepdata( 
                Rownum INT NOT NULL, 
                AttemptNum INT NOT NULL,
                BSID INT NOT NULL, 
                Command TEXT, 
                Arguments TEXT, 
                Started date, 
                Ended date,  
                Proc FLOAT,
                ExpProc FLOAT,
                Exit INT, 
                CPUs INT, 
                FOREIGN KEY(BSID) REFERENCES PulsarSearch(Rownum),
                primary key (Rownum, AttemptNum, BSID))""")
    
    cur.execute("""CREATE TABLE FFT( 
                Rownum INT NOT NULL, 
                AttemptNum INT NOT NULL,
                BSID INT NOT NULL, 
                Command TEXT, 
                Arguments TEXT, 
                Started date, 
                Ended date,  
                Proc FLOAT,
                ExpProc FLOAT,
                Exit INT, 
                CPUs INT,
                FOREIGN KEY(BSID) REFERENCES PulsarSearch(Rownum),
                primary key (Rownum, AttemptNum, BSID))""")
    
    cur.execute("""CREATE TABLE Accel( 
                Rownum INT NOT NULL, 
                AttemptNum INT NOT NULL,
                BSID INT NOT NULL, 
                Command TEXT, 
                Arguments TEXT, 
                Started date, 
                Ended date,  
                Proc FLOAT,
                ExpProc FLOAT,
                Exit INT, 
                CPUs INT,
                FOREIGN KEY(BSID) REFERENCES PulsarSearch(Rownum),
                primary key (Rownum, AttemptNum, BSID))""")
    
    cur.execute("""CREATE TABLE Fold( 
                Rownum INT NOT NULL, 
                AttemptNum INT NOT NULL,
                BSID INT NOT NULL, 
                Command TEXT, 
                Arguments TEXT, 
                Started date, 
                Ended date,  
                Proc FLOAT,
                ExpProc FLOAT,
                Exit INT, 
                CPUs INT,
                FOREIGN KEY(BSID) REFERENCES PulsarSearch(Rownum),
                primary key (Rownum, AttemptNum, BSID))""")
    
    #Result: if not yet checked by a human eye unchecked, otherwise under, RFI, Noise, Potential, Pulsar
    cur.execute("""CREATE TABLE Candidates(
                Rownum integer primary key autoincrement, 
                BSID INT, 
                FileLocation TEXT, 
                RESULT TEXT, 
                FOREIGN KEY(BSID) REFERENCES PulsarSearch(Rownum))""")

