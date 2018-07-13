import os
import sqlite3 as lite
try:
    DB_FILE = os.environ['CMD_BS_DB_FILE']
except:
    print "environmental variable JOB_DB_FILE must be defined"
    
con = lite.connect(DB_FILE)
with con:
    #Blindsearch Table: stores all major values for the blind search pipeline
    #Rownum: ID to reference
    #Comment: What the pipeline's purpose (eg testing or main search centered at x y)
    #*Proc: total crossing time in cpu hours
    #*Errors: total errors for each section
    #CandTotal: total cands
    #CandOverNoise: total cands over noise (normaly 3 sigma but I may change it)
    #CandDect: cands that were confirmed as pulsars (known or hopefully new)
    cur = con.cursor()
    cur.execute("CREATE TABLE Blindsearch(Rownum integer primary key autoincrement, Started date, Ended date, Obsid INT, Pointing TEXT, Comment TEXT, TotalProc FLOAT, TotalErrors INT, BeamformProc FLOAT, BeamformErrors INT, RFIProc FLOAT, RFIErrors INT, PrepdataProc FLOAT, PrepdataErrors INT, FFTProc FLOAT, FFTErrors INT, AccelProc FLOAT, AccelErrors INT, FoldProc FLOAT, FoldErrors INT, CandTotal INT, CandOverNoise INT, CandDect INT)")
    
    #A table for each job type
    
    #Rownum: ID to reference
    #BSID: Blindsearch id Rownum reference
    #Command: the code being run
    #Arguments: arguments put into the code
    #Exit: error code
    #CPUs: Number of cpus being used to be multiplied by processing time (end-start)
    
    cur.execute("CREATE TABLE Beamform(Rownum integer primary key autoincrement, BSID INT, Command TEXT, Arguments TEXT, Started date, Ended date, Exit INT, CPUs INT, FOREIGN KEY(BSID) REFERENCES Blindsearch(Rownum))")
    
    cur.execute("CREATE TABLE RFI(Rownum integer primary key autoincrement, BSID INT, Command TEXT, Arguments TEXT, Started date, Ended date, Exit INT, CPUs INT, FOREIGN KEY(BSID) REFERENCES Blindsearch(Rownum))")
    
    cur.execute("CREATE TABLE Prepdata(Rownum integer primary key autoincrement, BSID INT, Command TEXT, Arguments TEXT, Started date, Ended date, Exit INT, CPUs INT, FOREIGN KEY(BSID) REFERENCES Blindsearch(Rownum))")
    
    #DMFileInt: Dm file reference eg 1: DM_002_004
    cur.execute("CREATE TABLE FFT(Rownum integer primary key autoincrement, BSID INT, Command TEXT, Arguments TEXT, Started date, Ended date, Exit INT, CPUs INT, DMFileInt INT, FOREIGN KEY(BSID) REFERENCES Blindsearch(Rownum))")
    
    cur.execute("CREATE TABLE Accel(Rownum integer primary key autoincrement, BSID INT, Command TEXT, Arguments TEXT, Started date, Ended date, Exit INT, CPUs INT, DMFileInt INT, FOREIGN KEY(BSID) REFERENCES Blindsearch(Rownum))")
    
    cur.execute("CREATE TABLE Fold(Rownum integer primary key autoincrement, BSID INT, Command TEXT, Arguments TEXT, Started date, Ended date, Exit INT, CPUs INT, DMFileInt INT, FOREIGN KEY(BSID) REFERENCES Blindsearch(Rownum))")
    
    #Result: if not yet checked by a human eye unchecked, otherwise under, RFI, Noise, Potential, Pulsar
    cur.execute("CREATE TABLE Candidates(Rownum integer primary key autoincrement, BSID INT, DMFileInt INT, FileLocation TEXT, RESULT TEXT, FOREIGN KEY(BSID) REFERENCES Blindsearch(Rownum))")

