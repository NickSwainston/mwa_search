import os
import sqlite3 as lite
try:
    DB_FILE = os.environ['CMD_DB_FILE']
except:
    print "environmental variable JOB_DB_FILE must be defined"
    
con = lite.connect(DB_FILE)
with con:
    cur = con.cursor()
    # Datadir - path to main data store
    # Project - overarching identifier for this processing (e.g. ips/160303)
    # Obsid - Observation ID (gpsseconds)
    # JobID - ID of job which ran the command
    # TaskID - ID of Task which ran the command (NULL except for array jobs)
    # Command  - e.g. obsdownload, cotter, 
    # Arguments - all arguments to the task
    # Channels - coarse channel or channels that were operated on (e.g. 076-077) NULL means all for this obsid)
    # Username - user who ran the job
    # Started - start time
    # Ended - end time NULL implies that the job is running or was recently terminated 
    # Exit - Command exit code (NULL and (Ended != NULL) means the job was terminated.
    cur.execute("CREATE TABLE Commands(Rownum integer primary key autoincrement, Datadir TEXT, Project TEXT, Obsid INT, JobId INT, TaskID INT, Command TEXT, Channels TEXT, Arguments TEXT, UserId TEXT, Started date, Ended date, Exit INT)")
