#!/usr/bin/env python
import os, datetime, logging
import sqlite3 as lite
import glob
from optparse import OptionParser #NB zeus does not have argparse!

DB_FILE = os.environ['CMD_BS_DB_DEF_FILE']
#how many seconds the sqlite database conection takes until it times out
TIMEOUT=120

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d
    
def database_blindsearch_start(obsid, pointing, comment):
        con = lite.connect(DB_FILE, timeout = TIMEOUT)
        with con:
                cur = con.cursor()
                cur.execute("""INSERT INTO Blindsearch(Started, Obsid, Pointing, Comment, 
                        TotalProc, TotalErrors, TotalDS, TotalDE, TotalJobComp,
                        BeamformProc, BeamformErrors, BeamformDS, BeamformDE, BeamformJobComp,
                        PrepdataProc, PrepdataErrors, PrepdataDS, PrepdataDE, PrepdataJobComp,
                        FFTProc, FFTErrors, FFTDS, FFTDE, FFTJobComp,
                        AccelProc, AccelErrors, AccelDS, AccelDE, AccelJobComp,
                        FoldProc, FoldErrors, FoldDS, FoldDE, FoldJobComp,
                        CandTotal, CandOverNoise, CandDect) VALUES(?,?,?,?,
                               ?,?,?,?,?,
                               ?,?,?,?,?,
                               ?,?,?,?,?,
                               ?,?,?,?,?,
                               ?,?,?,?,?,
                               ?,?,?,?,?,
                               ?,?,?)""",
                          (datetime.datetime.now(), obsid, pointing, comment,
                          0.0, 0, 0, 0, 0,
                          0.0, 0, 0, 0, 0,
                          0.0, 0, 0, 0, 0,
                          0.0, 0, 0, 0, 0,
                          0.0, 0, 0, 0, 0,
                          0.0, 0, 0, 0, 0,
                          0, 0, 0))
                vcs_command_id = cur.lastrowid
        return vcs_command_id
    

def database_script_list(bs_id, command, arguments_list, threads, expe_proc_time,
                         attempt=1):
    """
    Will create all the rows in the database for each job
    """
    #works out the table from the command
    if command == 'make_beam':
        table = 'Beamform'
    if command == 'prepsubband':
        table = 'Prepdata'
    elif command == 'realfft':
        table = 'FFT'
    elif command == 'accelsearch':
        table = 'Accel'
    elif command == 'prepfold':
        table = 'Fold'
    con = lite.connect(DB_FILE, timeout = TIMEOUT)
    with con:
        cur = con.cursor()
        row_id_list = []
        for ai, arguments in enumerate(arguments_list):
            cur.execute("INSERT OR IGNORE INTO "+table+" (Rownum, AttemptNum, BSID, Command, Arguments, CPUs, ExpProc) VALUES(?, ?, ?, ?, ?, ?, ?)", (ai, attempt, bs_id, command, arguments, threads, expe_proc_time))
        #update expected jobs
        if attempt == 1:
            cur.execute("UPDATE Blindsearch SET "+table+"JobExp=? WHERE Rownum=?", (str(len(arguments_list)),bs_id))
        else:
            cur.execute("SELECT "+table+"JobExp FROM Blindsearch WHERE Rownum=?", (str(bs_id)))
            table_job_exp = cur.fetchone()[0]
            cur.execute("UPDATE Blindsearch SET "+table+"JobExp=? WHERE Rownum=?", (str(len(arguments_list) + table_job_exp),bs_id))
        cur.execute("SELECT TotalJobExp FROM Blindsearch WHERE Rownum={}".format(bs_id))
        blind_job_exp = cur.fetchone()[0]
        if blind_job_exp is None:
            blind_job_exp = 0
        cur.execute("UPDATE Blindsearch SET TotalJobExp=? WHERE Rownum=?", (str(len(arguments_list) + blind_job_exp),bs_id))

    return 


def database_script_start(table, bs_id, rownum, attempt_num, time=datetime.datetime.now()):
    
    con = lite.connect(DB_FILE, timeout = TIMEOUT)
    with con:
        cur = con.cursor()
        cur.execute("UPDATE "+table+" SET Started=? WHERE Rownum=? AND AttemptNum=? AND BSID=?",
                    (time, rownum, attempt_num, bs_id))
        row_id = cur.lastrowid
    return row_id

def database_script_stop(table, bs_id, rownum, attempt_num, errorcode,
                         end_time=datetime.datetime.now()):

    con = lite.connect(DB_FILE, timeout = TIMEOUT)
    con.row_factory = lite.Row
    with con:
        cur = con.cursor()
        #get script data
        cur.execute("SELECT * FROM "+table+" WHERE Rownum=? AND AttemptNum=? AND BSID=?",
                    (rownum, attempt_num, bs_id))
        columns = cur.fetchone()
        #get blindsearch data
        cur.execute("SELECT * FROM Blindsearch WHERE Rownum="+str(bs_id))
        bs_columns = cur.fetchone()

        if int(errorcode) == 0:
            #add processing times and job completion count
            end_s = date_to_sec(str(end_time))
            start_s = date_to_sec(columns['Started'])
            processing = (end_s - start_s)

            cur.execute("UPDATE "+table+\
                        " SET Proc=?, Ended=?, Exit=? WHERE Rownum=? AND AttemptNum=? AND BSID=?",
                                (processing, end_time, errorcode, rownum, attempt_num, bs_id))

            tot_proc = float(bs_columns['TotalProc']) + processing
            job_proc = float(bs_columns[table+'Proc']) + processing
            tot_jc = int(bs_columns['TotalJobComp']) + 1
            job_jc = int(bs_columns[table+'JobComp']) + 1

            cur.execute("UPDATE Blindsearch SET TotalProc=?, "+table+\
                        "Proc=?, TotalJobComp=?, "+table+\
                        "JobComp=? WHERE Rownum=?",
                        (str(tot_proc)[:9], str(job_proc)[:9], str(tot_jc)[:9],
                         str(job_jc)[:9], bs_id))
        else:    
            tot_er = int(bs_columns['TotalErrors']) + 1
            job_er = int(bs_columns[table+'Errors']) + 1

            cur.execute("UPDATE "+table+\
                        " SET Ended=?, Exit=? WHERE Rownum=? AND AttemptNum=? AND BSID=?",
                              (end_time, errorcode, rownum, attempt_num, bs_id))
                
            cur.execute("UPDATE Blindsearch SET TotalErrors=?, "+table+\
                        "Errors=? WHERE Rownum=?",
                        (tot_er,job_er, bs_id))
    return


def database_script_check(table, bs_id, attempt_num):
    """
    Searches for any jobs that didn't work and return the data needed to send 
    them off again
    """
    con = lite.connect(DB_FILE, timeout = TIMEOUT)
    con.row_factory = lite.Row
    with con:
        cur = con.cursor()
        #get script data
        cur.execute("SELECT * FROM "+table+" WHERE AttemptNum=? AND BSID=?",
                    (attempt_num, bs_id))
        rows = cur.fetchall()
        
        error_data = []
        for row in rows:
            if row['Started'] == None or row['Ended'] == None or row['Exit'] != 0:
                error_data.append([row['Command'], row['Arguments'], row['ExpProc']])
    return error_data            


def database_mass_update(table,file_location):
    """
    Takes a csv from short period jobs and makes a large amount of updates
    """
    with open(file_location,'r') as csv:
        con = lite.connect(DB_FILE, timeout = TIMEOUT)
        con.row_factory = lite.Row
        with con:
            cur = con.cursor()
            lines = csv.readlines()
            for i,l in enumerate(lines):
                l = l.split(',')
                if len(l) > 2:
                    started = l[0]
                    rownum = l[2]
                    attempt_num = l[3]
                    bs_id = l[1]
                
                else:
                    end_time = l[0]
                    errorcode = l[1]

                    cur.execute("UPDATE "+table+\
                            " SET Started=? WHERE Rownum=? AND AttemptNum=? AND BSID=?",
                            (started, rownum, attempt_num, bs_id))

                    cur.execute("SELECT * FROM "+table+" WHERE Rownum=? AND AttemptNum=? AND BSID=?",
                                (rownum, attempt_num, bs_id))
                    columns = cur.fetchone()
                    #get blindsearch data
                    cur.execute("SELECT * FROM Blindsearch WHERE Rownum="+str(bs_id))
                    bs_columns = cur.fetchone()

                    if int(errorcode) == 0:
                        #add processing times and job completion count
                        end_s = date_to_sec(str(end_time))
                        start_s = date_to_sec(columns['Started'])
                        processing = (end_s - start_s)

                        cur.execute("UPDATE "+table+\
                                    " SET Proc=?, Ended=?, Exit=? WHERE Rownum=? AND AttemptNum=? AND BSID=?",
                                            (processing, end_time, errorcode, rownum, attempt_num, bs_id))

                        tot_proc = float(bs_columns['TotalProc']) + processing
                        job_proc = float(bs_columns[table+'Proc']) + processing
                        tot_jc = int(bs_columns['TotalJobComp']) + 1
                        job_jc = int(bs_columns[table+'JobComp']) + 1

                        cur.execute("UPDATE Blindsearch SET TotalProc=?, "+table+\
                                    "Proc=?, TotalJobComp=?, "+table+\
                                    "JobComp=? WHERE Rownum=?",
                                    (str(tot_proc)[:9], str(job_proc)[:9], str(tot_jc)[:9],
                                     str(job_jc)[:9], bs_id))
                    else:    
                        tot_er = int(bs_columns['TotalErrors']) + 1
                        job_er = int(bs_columns[table+'Errors']) + 1

                        cur.execute("UPDATE "+table+\
                                    " SET Ended=?, Exit=? WHERE Rownum=? AND AttemptNum=? AND BSID=?",
                                          (end_time, errorcode, rownum, attempt_num, bs_id))
                            
                        cur.execute("UPDATE Blindsearch SET TotalErrors=?, "+table+\
                                    "Errors=? WHERE Rownum=?",
                                    (tot_er,job_er, bs_id))
    return


def database_wrap_up(rownum, cand_val, end_time=datetime.datetime.now()):
    """
    Updated the cand numbers
    """
    cand_total, cand_over_noise, cand_detect = opts.cand_val.split()
    con = lite.connect(DB_FILE, timeout = TIMEOUT)
    with con:
        cur = con.cursor()
        cur.execute("UPDATE Blindsearch SET Ended=?, CandTotal=?, CandOverNoise=?, CandDect=? WHERE Rownum=?",
                    (end_time, cand_total, cand_over_noise, cand_detect, rownum))
    return


def database_beamform_find(file_location, bs_id):
    time_now = datetime.datetime.now()
    #go through the batch file for info
    arguments_list = []
    batch_file_list = glob.glob('{0}*.batch'.format(file_location))
    for batch_file in batch_file_list:
        with open(batch_file,'r') as batch:
            lines = batch.readlines()
            for l in lines:
                if l.startswith("#SBATCH --time="):
                    time_str = l.split("=")[1]
                    hr, mi, se = time_str.split(":")
                    expe_proc_time = float(se) + float(mi)*60. + float(hr)*3600.
                if l.startswith("export OMP_NUM_THREADS"):
                    nodes = l.split("=")[1]
                if l.startswith("srun"):
                    command = "make_beam" #I don't think these needs to be more robust
                    arguments_list.append(l.split(command)[1])
    #set up the beamform database
    database_script_list(bs_id, 'make_beam', arguments_list, nodes, expe_proc_time)
    
    #go through the output files for start stop times
    out_file_list = glob.glob('{0}*.out'.format(file_location))
    for rownum, out_file in enumerate(out_file_list):
        with open(out_file,'r') as output:
            lines = output.readlines()
            find_check = False
            for l in lines:
                if "**FINISHED BEAMFORMING**" in l:
                    find_check = True
                    time_seconds = float(l[1:10])
                    #So this may be inaccurate because now isn't when the job 
                    #finished but should get the right delta
                    time_then = datetime.datetime.now() - datetime.timedelta(seconds=time_seconds)
                    #TODO add attemp number options
                    database_script_start('Beamform', bs_id, rownum, 1, time=time_then)
                    database_script_stop('Beamform', bs_id, rownum, 1, 0, end_time=time_now)
            if not find_check:
                #no finshed string so likely it failed:
                #TODO make this more robust to work out how long it ran before it died
                database_script_start('Beamform', bs_id, rownum, 1, time=time_now)
                database_script_stop('Beamform', bs_id, rownum, 1, 1, end_time=time_now)
    return 

def date_to_sec(string):
    #just an approximation (doesn't even use year and month
    date, time = string.split(' ')
    y,m,d = date.split('-')
    h,mi,s = time.split(':')
    s_out = ((float(d)*24. + float(h))*60. + float(mi))*60. + float(s)
    #print "hours " + str(float(d)*24. + float(h))
    #print "Minutes " + str((float(d)*24. + float(h))*60. + float(mi))
    #print "secounds " + str(s_out)
    return s_out
    

if __name__ == '__main__':
    from optparse import OptionParser, OptionGroup, SUPPRESS_HELP
    parser = OptionParser(usage = "usage: %prog <options>" +
    """
    Script used to manage the VCS database by recording the scripts process_vcs.py uses and prints the database.
    Common commands:
    blindsearch_database.py -m vc
    blindsearch_database.py -m vs -c <presto command>
    
    """)
    parser.add_option("-m", "--mode", dest="mode", metavar="mode", default='vc', type=str, help='This script has three modes: "vc" used to view the database commands, "vs" used to view the database scripts, "vp" view processing and error statistics, "s" used to start a record of a script on the database, "e" used to record the end time and error code of a script on the database, "p" counts errors and processing time for one id and "b" is a special mode for receiving the total time of beamforming jobs. Default mode is v')
    parser.add_option("-f", "--file_location", dest="file_location", metavar="file_location", type=str, help='mass update csv file location.')
    parser.add_option("--cand_val", dest="cand_val", metavar="cand_val", type=str, help='Cand values in the form of "<total> <overnoise> <detected>. Cand detected not yet implimented.')
    
    view_options = OptionGroup(parser, 'View Options')
    view_options.add_option("--recent", dest="recent", metavar="HOURS", default=None, type=float, help="print only jobs started in the last N hours")
    view_options.add_option("--number", dest="n", metavar="N", default=20, type=int, help="number of jobs to print [default=%default]")
    view_options.add_option("--all", dest="all", action="store_true", help="print all lines of the database")
    view_options.add_option("-s", "--startrow", dest="startrow", default=0, type=int, help="ignore any row earlier than this one")
    view_options.add_option("-e", "--endrow", dest="endrow", default=None, type=int, help="ignore any row later than this one")
    view_options.add_option("-o", "--obsid", dest="obsid", default=None, type=str, help="Only prints one obsid's jobs.")
    
    start_options = OptionGroup(parser, 'Script Start Options')
    start_options.add_option("-b", "--bs_id", dest="bs_id", default=None, type=str, help="The row number of the blindsearch command of the databse")
    start_options.add_option("-c", "--command", dest="command", default=None, type=str, help="The script name being run. eg volt_download.py.")
    start_options.add_option("-a", "--attempt_num", dest="attempt_num", default=None, type=str, help="The attempt number of a script.")
    start_options.add_option("-n", "--nodes", dest="nodes", default=None, type=int, help="The number of cpu nodes used.")
    start_options.add_option("-d", "--dm_file_int", dest="dm_file_int", default=None, type=int, help="The DM file reference eg 1 = DM_002_004.")
    
    end_options = OptionGroup(parser, 'Script End Options')
    end_options.add_option("--errorcode", dest="errorcode", default=None, type=int, help="Error code of scripts.")
    end_options.add_option("-r", "--rownum", dest="rownum", default=None, type=str, help="The row number of the script.")
    parser.add_option_group(view_options)
    parser.add_option_group(start_options)
    parser.add_option_group(end_options)
    (opts, args) = parser.parse_args()
    
    #work out table
    if opts.command == 'rfifind':
        table = 'RFI'
    elif opts.command == 'prepsubband':
        table = 'Prepdata'
    elif opts.command == 'realfft':
        table = 'FFT'
    elif opts.command == 'accelsearch':
        table = 'Accel'
    elif opts.command == 'prepfold':
        table = 'Fold'
    elif opts.mode == 'vc' or opts.mode == 'vp' or opts.mode == 'vprog':
        table = 'Blindsearch'
    elif opts.mode == 'b' or opts.command == 'make_beam':
        table = 'Beamform'
        
    
    if opts.mode == "s":
        vcs_row = database_script_start(table, opts.bs_id, opts.rownum, opts.attempt_num)
    elif opts.mode == "e":
        database_script_stop(table, opts.bs_id, opts.rownum, opts.attempt_num, opts.errorcode)
    elif opts.mode == 'w':
        database_wrap_up(opts.bs_id, opts.cand_val)
    elif opts.mode == 'm':
        if opts.file_location:
            file_loc = opts.file_location
        else:
            file_loc = opts.command + '_temp_database_file.csv'
        database_mass_update(table,file_loc)
    elif opts.mode == 'b':
        database_beamform_find(opts.file_location, opts.bs_id)
    elif opts.mode.startswith("v"):
        con = lite.connect(DB_FILE, timeout = TIMEOUT)
        con.row_factory = dict_factory
    
        query = "SELECT * FROM " + table

        if opts.obsid:
            query += " WHERE Arguments LIKE '%" + str(opts.obsid) + "%'"

        if opts.recent is not None:
            query += ''' WHERE Started > "%s"''' % str(datetime.datetime.now() - relativedelta(hours=opts.recent))
            logging.debug(query)
        if opts.bs_id and opts.mode == 'vs':
            query += " WHERE BSID=" + str(opts.bs_id)
        elif opts.bs_id and not opts.mode == 'vs':
            query += " WHERE Rownum=" + str(opts.bs_id) 
            
        if opts.dm_file_int:
            query += " WHERE DMFileInt='" + str(opts.dm_file_int) + "'"
        
        if opts.attempt_num:
            if "WHERE" in query:
                query += " AND AttemptNum='" + str(opts.attempt_num) + "'"
            else:
                query += " WHERE AttemptNum='" + str(opts.attempt_num) + "'"
        
        if opts.errorcode:
            if "WHERE" in query:
                query += " AND Exit='" + str(opts.errorcode) + "'"
            else:
                query += " WHERE Exit='" + str(opts.errorcode) + "'"

        print query
        with con:
            cur = con.cursor()
            cur.execute(query)
            rows = cur.fetchall()

        if opts.startrow and opts.endrow is None:
            rows = rows[opts.startrow:]
        elif opts.endrow is not None:
            rows = rows[opts.startrow:opts.endrow+1]
        elif not (opts.all or opts.recent):
            rows = rows[-opts.n:]
        
        
        if opts.mode == "vc": 
            print 'Row# ','Obsid       ','Pointing                      ','Started               ','Ended                 ','Comments'
            print '--------------------------------------------------------------------------------------------------'
            for row in rows:
                print '%-5s' % (str(row['Rownum']).rjust(4)),
                print '%-12s' % (row['Obsid']),
                print '%-30s' % (row['Pointing']),
                print '%-22s' % (row['Started'][:19]),
                if row['Ended'] is None:
                    print '%-22s' % (row['Ended']), 
                else:
                    print '%-22s' % (row['Ended'][:19]),
                print row['Comment']
                #print "\n"
                
                
        if opts.mode == "vs":
            print 'BDIS ','Row# ','Atm#','Started               ','Ended                 ','Exit_Code','ProcTime ','ExpecTime ','Arguments'
            print '--------------------------------------------------------------------------------------------------'
            for row in rows:
                #BSID INT, Command TEXT, Arguments TEXT, Started date, Ended date, Exit
                print '%-5s' % (row['BSID']),
                print '%-5s' % (str(row['Rownum']).rjust(4)),
                print '%-5s' % (row['AttemptNum']),
                if row['Started'] is None:
                    print '%-22s' % (row['Started']),
                else:
                    print '%-22s' % (row['Started'][:19]),
                if row['Ended'] is None:
                    print '%-22s' % (row['Ended']),
                else:
                    print '%-22s' % (row['Ended'][:19]),
                print '%-7s' % (row['Proc']),
                print '%-7s' % (row['ExpProc']),
                if str(row['Exit']).endswith('\n'):
                    print '%-5s' % str(row['Exit'])[:-1],
                else:
                    print '%-5s' % (row['Exit']),
                print '%-5s' % (row['CPUs']),
                print row['Arguments'],
                print "\n"
                
        if opts.mode == "vp":
            for ri, row in enumerate(rows):
                if ri%20 == 0:
                    print 'Row# | Total proc | err# | Beamform proc | err# | Prep proc | err# | FFT proc | err# | Accel proc | err# | Fold proc | err# |'
                    print '-----|------------|------|---------------|------|-----------|------|----------|------|------------|------|-----------|------|'

                #TotalProc FLOAT, TotalErrors INT, RFIProc FLOAT, RFIErrors INT, PrepdataProc FLOAT, PrepdataErrors INT, FFTProc FLOAT, FFTErrors INT, AccelProc FLOAT, AccelErrors INT, FoldProc FLOAT, FoldErrors INT,
                print '{:4s} |{:11.2f} |{:5d} | {:13.2f} |{:5d} | {:9.2f} |{:5d} | {:8.2f} |{:5d} | {:10.2f} |{:5d} | {:9.2f} |{:5d} |'.format(str(row['Rownum']).rjust(4),row['TotalProc']/3600.,row['TotalErrors'],row['BeamformProc']/3600.,row['BeamformErrors'],row['PrepdataProc']/3600.,row['PrepdataErrors'],row['FFTProc']/3600.,row['FFTErrors'],row['AccelProc']/3600.,row['AccelErrors'],row['FoldProc']/3600.,row['FoldErrors'])
        
        if opts.mode == "vprog":
            for ri, row in enumerate(rows):
                if ri%20 == 0:
                    print 'Row# |  Total Jobs |Beamform Jobs|Prepdata Jobs|'+\
                            '   FFT Jobs  |  Accel Jobs |  Fold Jobs  |'
                    print '-----|-------------|-------------|-------------|'+\
                            '-------------|-------------|-------------|'
                print '{:4d} | {:5d}/{:5s} | {:5d}/{:5s} | {:5d}/{:5s} | {:5d}/{:5s} | {:5d}/{:5s} | {:5d}/{:5s} |'.\
                        format(row['Rownum'], row['TotalJobComp'], str(row['TotalJobExp']),
                               row['BeamformJobComp'], str(row['BeamformJobExp']),
                               row['PrepdataJobComp'], str(row['PrepdataJobExp']),
                               row['FFTJobComp'], str(row['FFTJobExp']),
                               row['AccelJobComp'], str(row['AccelJobExp']),
                               row['FoldJobComp'], str(row['FoldJobExp']) )
       
