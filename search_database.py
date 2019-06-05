#!/usr/bin/env python
import os, datetime, logging
import sqlite3 as lite
import glob
import argparse
import textwrap

DB_FILE = os.environ['SEARCH_DB']
#how many seconds the sqlite database conection takes until it times out
TIMEOUT=120

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d
    
def database_search_start(obsid, pointing, comment):
    import version
    con = lite.connect(DB_FILE, timeout = TIMEOUT)
    with con:
        cur = con.cursor()
        cur.execute("""INSERT INTO PulsarSearch(Started, Obsid, Pointing,
                UserID, Version, Comment, 
                TotalProc, TotalErrors, TotalDS, TotalDE, TotalJobComp,
                BeamformProc, BeamformErrors, BeamformDS, BeamformDE, BeamformJobComp,
                PrepdataProc, PrepdataErrors, PrepdataDS, PrepdataDE, PrepdataJobComp,
                FFTProc, FFTErrors, FFTDS, FFTDE, FFTJobComp,
                AccelProc, AccelErrors, AccelDS, AccelDE, AccelJobComp,
                FoldProc, FoldErrors, FoldDS, FoldDE, FoldJobComp,
                CandTotal, CandOverNoise, CandDect) VALUES(?,?,?,?,?,?,
                       ?,?,?,?,?,
                       ?,?,?,?,?,
                       ?,?,?,?,?,
                       ?,?,?,?,?,
                       ?,?,?,?,?,
                       ?,?,?,?,?,
                       ?,?,?)""",
                  (datetime.datetime.now(), obsid, pointing, 
                  os.environ['USER'], version.__version__, comment,
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
        for ai, arguments in enumerate(arguments_list):
            cur.execute("INSERT OR IGNORE INTO %s (Rownum, AttemptNum, BSID, Command, Arguments, CPUs, ExpProc) VALUES(?, ?, ?, ?, ?, ?, ?)" % table, (ai, attempt, bs_id, command, arguments, threads, expe_proc_time))
        #update expected jobs
        if attempt == 1:
            cur.execute("UPDATE PulsarSearch SET %sJobExp=? WHERE Rownum=?" % table, (len(arguments_list),bs_id))
        else:
            cur.execute("SELECT %sJobExp FROM PulsarSearch WHERE Rownum=?" % table, (bs_id,))
            table_job_exp = cur.fetchone()[0]
            cur.execute("UPDATE PulsarSearch SET %sJobExp=? WHERE Rownum=?" % table, (len(arguments_list) + table_job_exp, bs_id))
        cur.execute("SELECT TotalJobExp FROM PulsarSearch WHERE Rownum=?", (bs_id,))
        search_job_exp = cur.fetchone()[0]
        if search_job_exp is None:
            search_job_exp = 0
        cur.execute("UPDATE PulsarSearch SET TotalJobExp=? WHERE Rownum=?", (len(arguments_list) + search_job_exp, bs_id))

    return 


def database_script_start(table, bs_id, rownum, attempt_num, time=datetime.datetime.now()):
    
    con = lite.connect(DB_FILE, timeout = TIMEOUT)
    with con:
        cur = con.cursor()
        cur.execute("UPDATE %s SET Started=? WHERE Rownum=? AND AttemptNum=? AND BSID=?" % table,
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
        cur.execute("SELECT * FROM %s WHERE Rownum=? AND AttemptNum=? AND BSID=?" % table,
                    (rownum, attempt_num, bs_id))
        columns = cur.fetchone()
        #get search data
        cur.execute("SELECT * FROM PulsarSearch WHERE Rownum=?",(bs_id,))
        bs_columns = cur.fetchone()

        if int(errorcode) == 0:
            #add processing times and job completion count
            end_s = date_to_sec(str(end_time))
            start_s = date_to_sec(columns['Started'])
            processing = (end_s - start_s)

            cur.execute("UPDATE %s SET Proc=?, Ended=?, Exit=? WHERE Rownum=? AND AttemptNum=? AND BSID=?" % table, (processing, end_time, errorcode, rownum, attempt_num, bs_id))

            tot_proc = float(bs_columns['TotalProc']) + processing
            job_proc = float(bs_columns[table+'Proc']) + processing
            tot_jc = int(bs_columns['TotalJobComp']) + 1
            job_jc = int(bs_columns[table+'JobComp']) + 1

            cur.execute("UPDATE PulsarSearch SET TotalProc=?, ?Proc=?, TotalJobComp=?, ?JobComp=? WHERE Rownum=?", (str(tot_proc)[:9], table, str(job_proc)[:9], str(tot_jc)[:9], table,
                         str(job_jc)[:9], bs_id))
        else:    
            tot_er = int(bs_columns['TotalErrors']) + 1
            job_er = int(bs_columns[table+'Errors']) + 1

            cur.execute("UPDATE %s SET Ended=?, Exit=? WHERE Rownum=? AND AttemptNum=? AND BSID=?" % table,
                              (end_time, errorcode, rownum, attempt_num, bs_id))
                
            cur.execute("UPDATE PulsarSearch SET TotalErrors=?, ?Errors=? WHERE Rownum=?",
                        (tot_er, table, job_er, bs_id))
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
        cur.execute("SELECT * FROM %s WHERE AttemptNum=? AND BSID=?" % table,
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
            for l in lines:
                l = l.split(',')
                if len(l) > 2:
                    started = l[0]
                    rownum = l[2]
                    attempt_num = l[3]
                    bs_id = l[1]
                
                else:
                    end_time = l[0]
                    errorcode = l[1]

                    cur.execute("UPDATE %s SET Started=? WHERE Rownum=? AND AttemptNum=? AND BSID=?" % table, (started, rownum, attempt_num, bs_id))

                    cur.execute("SELECT * FROM %s WHERE Rownum=? AND AttemptNum=? AND BSID=?" % table, (rownum, attempt_num, bs_id))
                    columns = cur.fetchone()
                    #get search data
                    cur.execute("SELECT * FROM PulsarSearch WHERE Rownum=?", (str(bs_id),))
                    bs_columns = cur.fetchone()

                    if int(errorcode) == 0:
                        #add processing times and job completion count
                        end_s = date_to_sec(str(end_time))
                        start_s = date_to_sec(columns['Started'])
                        processing = (end_s - start_s)

                        cur.execute("UPDATE %s SET Proc=?, Ended=?, Exit=? WHERE Rownum=? AND AttemptNum=? AND BSID=?" % table, (processing, end_time, errorcode, rownum, attempt_num, bs_id))

                        tot_proc = float(bs_columns['TotalProc']) + processing
                        job_proc = float(bs_columns[table+'Proc']) + processing
                        tot_jc = int(bs_columns['TotalJobComp']) + 1
                        job_jc = int(bs_columns[table+'JobComp']) + 1

                        cur.execute("UPDATE PulsarSearch SET TotalProc=?, %sProc=?, TotalJobComp=?, %sJobComp=? WHERE Rownum=?" % table,
                                    (str(tot_proc)[:9], str(job_proc)[:9], str(tot_jc)[:9],
                                     str(job_jc)[:9], bs_id))
                    else:
                        tot_er = int(bs_columns['TotalErrors']) + 1
                        job_er = int(bs_columns[table+'Errors']) + 1

                        cur.execute("UPDATE %s SET Ended=?, Exit=? WHERE Rownum=? AND AttemptNum=? AND BSID=?" % table, (end_time, errorcode, rownum, attempt_num, bs_id))
                            
                        cur.execute("UPDATE PulsarSearch SET TotalErrors=?, ?Errors=? WHERE Rownum=?",
                                    (tot_er,job_er, bs_id))
    return


def database_wrap_up(rownum, cand_val, end_time=datetime.datetime.now()):
    """
    Updated the cand numbers
    """
    cand_total, cand_over_noise, cand_detect = args.cand_val.split()
    con = lite.connect(DB_FILE, timeout = TIMEOUT)
    with con:
        cur = con.cursor()
        cur.execute("UPDATE PulsarSearch SET Ended=?, CandTotal=?, CandOverNoise=?, CandDect=? WHERE Rownum=?",
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
    _, _, d = date.split('-')
    h, mi, s = time.split(':')
    s_out = ((float(d)*24. + float(h))*60. + float(mi))*60. + float(s)
    #print "hours " + str(float(d)*24. + float(h))
    #print "Minutes " + str((float(d)*24. + float(h))*60. + float(mi))
    #print "secounds " + str(s_out)
    return s_out
    

if __name__ == '__main__':
    mode_options = ['vc', 'vs', 'vp', 'vprog', 's', 'e', 'm', 'b', 'w']
    parser = argparse.ArgumentParser(description="""A script used to keep track of the scripts run by the pulsar search database.""",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-m", "--mode", dest="mode", metavar="mode", default='vc', type=str, help=textwrap.dedent('''This script has the following modes {}. All modes starting with v are used for viewing parts of the database, other options will be used by the pipeline not the user.
"vc" used to view the commands that start the pipeline.
"vs" used to view the the individual scripts processing time and arguments.
"vp" used to view the total processing for each command type in walltime hours.
"vprog" used to view the progress of how many scripts have finished.
"s" used to record the start time of a script.
"e" used to record the end time, error code and processing time of a script.
"m" used to record the start and stop time using an csv file. Used for quick commands.
"b" is a special mode for receiving the total time of beamforming jobs. 
"w" is a wrap up for the pipeline that records candidate statistics and the end time.
Default mode is vc'''.format(mode_options)))
    
    view_options = parser.add_argument_group('View Options')
    view_options.add_argument("--recent", default=None, type=float, help="Print only jobs started in the last N hours")
    view_options.add_argument("--number", default=20, type=int, help="Number of jobs to print. Default is 20")
    view_options.add_argument("--all", action="store_true", help="Print all lines of the database")
    view_options.add_argument("-s", "--startrow", default=0, type=int, help="Ignore any row earlier than this one")
    view_options.add_argument("-e", "--endrow", default=None, type=int, help="Ignore any row later than this one")
    view_options.add_argument("-o", "--obsid", default=None, type=str, help="Only prints one obsid's jobs.")
    
    start_options = parser.add_argument_group('Script Start Options')
    start_options.add_argument("-b", "--bs_id", default=None, type=str, help="The row number of the pulsar search command of the databse")
    start_options.add_argument("-c", "--command", default=None, type=str, help="The script name being run. eg volt_download.py.")
    start_options.add_argument("-a", "--attempt_num", default=None, type=str, help="The attempt number of a script.")
    start_options.add_argument("-n", "--nodes", default=1, type=int, help="The number of cpu nodes used.")
    
    end_options = parser.add_argument_group('Script End Options')
    end_options.add_argument("--errorcode", dest="errorcode", default=None, type=int, help="Error code of scripts.")
    end_options.add_argument("-r", "--rownum", dest="rownum", default=None, type=str, help="The row number of the script.")
    end_options.add_argument("-f", "--file_location", type=str, help='Mode "m" csv file location.')
    end_options.add_argument("--cand_val", type=str, help='Cand values in the form of "<total> <overnoise> <detected>. Cand detected not yet implimented.')
    args=parser.parse_args()

    #argument parsing
    if args.mode not in mode_options:
        print("Unrecognised mode, please use one of the following {}. Exiting".format(mode_options))
        quit()

    #work out table
    if args.command == 'rfifind':
        table = 'RFI'
    elif args.command == 'prepsubband':
        table = 'Prepdata'
    elif args.command == 'realfft':
        table = 'FFT'
    elif args.command == 'accelsearch':
        table = 'Accel'
    elif args.command == 'prepfold':
        table = 'Fold'
    elif args.mode == 'vc' or args.mode == 'vp' or args.mode == 'vprog':
        table = 'PulsarSearch'
    elif args.mode == 'b' or args.command == 'make_beam':
        table = 'Beamform'
        
    
    if args.mode == "s":
        vcs_row = database_script_start(table, args.bs_id, args.rownum, args.attempt_num)
    elif args.mode == "e":
        database_script_stop(table, args.bs_id, args.rownum, args.attempt_num, args.errorcode)
    elif args.mode == 'w':
        database_wrap_up(args.bs_id, args.cand_val)
    elif args.mode == 'm':
        if args.file_location:
            file_loc = args.file_location
        else:
            file_loc = args.command + '_temp_database_file.csv'
        database_mass_update(table,file_loc)
    elif args.mode == 'b':
        database_beamform_find(args.file_location, args.bs_id)
    elif args.mode.startswith("v"):
        con = lite.connect(DB_FILE, timeout = TIMEOUT)
        con.row_factory = dict_factory
    
        query = "SELECT * FROM %s" % table

        if args.obsid:
            query += " WHERE Arguments LIKE '%" + str(args.obsid) + "%'"

        if args.recent is not None:
            query += ''' WHERE Started > "%s"''' % str(datetime.datetime.now() - relativedelta(hours=args.recent))
            logging.debug(query)
        if args.bs_id and args.mode == 'vs':
            query += " WHERE BSID=" + str(args.bs_id)
        elif args.bs_id and not args.mode == 'vs':
            query += " WHERE Rownum=" + str(args.bs_id)
            
        if args.attempt_num:
            if "WHERE" in query:
                query += " AND AttemptNum='" + str(args.attempt_num) + "'"
            else:
                query += " WHERE AttemptNum='" + str(args.attempt_num) + "'"
        
        if args.errorcode:
            if "WHERE" in query:
                query += " AND Exit='" + str(args.errorcode) + "'"
            else:
                query += " WHERE Exit='" + str(args.errorcode) + "'"

        print(query)
        with con:
            cur = con.cursor()
            cur.execute(query)
            rows = cur.fetchall()

        if args.startrow and args.endrow is None:
            rows = rows[args.startrow:]
        elif args.endrow is not None:
            rows = rows[args.startrow:args.endrow+1]
        elif not (args.all or args.recent):
            rows = rows[-args.number:]
        
        
        if args.mode == "vc":
            print('{:4s} | {:10s} | {:26s} | {:19s} | {:19s} | {}'.format('Row',
                  'Obsid', 'Pointing', 'Started', 'Ended', 'Comments'))
            print('-----------------------------------------------------------------------------------------------------')
            for row in rows:
                if row['Ended'] is None:
                    temp_ended = '{:19s}'.format('None')
                else:
                    temp_ended = '{:.19}'.format(row['Ended'])
                print('{:4d} | {:10d} | {:26s} | {:19s} | {:19s} | {}'.format(row['Rownum'], row['Obsid'], row['Pointing'], row['Started'][:19], temp_ended, row['Comment']))
                
                
        if args.mode == "vs":
            print('BDIS ','Row# ','Atm#','Started               ','Ended                 ','Exit_Code','ProcTime ','ExpecTime ','Arguments')
            print('--------------------------------------------------------------------------------------------------')
            for row in rows:
                #BSID INT, Command TEXT, Arguments TEXT, Started date, Ended date, Exit
                print('%-5s' % (row['BSID']),)
                print('%-5s' % (str(row['Rownum']).rjust(4)),)
                print('%-5s' % (row['AttemptNum']),)
                if row['Started'] is None:
                    print('%-22s' % (row['Started']),)
                else:
                    print('%-22s' % (row['Started'][:19]),)
                if row['Ended'] is None:
                    print('%-22s' % (row['Ended']),)
                else:
                    print('%-22s' % (row['Ended'][:19]),)
                print('%-7s' % (row['Proc']),)
                print('%-7s' % (row['ExpProc']),)
                if str(row['Exit']).endswith('\n'):
                    print('%-5s' % str(row['Exit'])[:-1],)
                else:
                    print('%-5s' % (row['Exit']),)
                print('%-5s' % (row['CPUs']),)
                print(row['Arguments'])
                
        if args.mode == "vp":
            for ri, row in enumerate(rows):
                if ri%20 == 0:
                    print('Row# | Total proc | err# | Beamform proc | err# | Prep proc | err# | FFT proc | err# | Accel proc | err# | Fold proc | err# |')
                    print('-----|------------|------|---------------|------|-----------|------|----------|------|------------|------|-----------|------|')

                #TotalProc FLOAT, TotalErrors INT, RFIProc FLOAT, RFIErrors INT, PrepdataProc FLOAT, PrepdataErrors INT, FFTProc FLOAT, FFTErrors INT, AccelProc FLOAT, AccelErrors INT, FoldProc FLOAT, FoldErrors INT,
                print('{:4s} |{:11.2f} |{:5d} | {:13.2f} |{:5d} | {:9.2f} |{:5d} | {:8.2f} |{:5d} | {:10.2f} |{:5d} | {:9.2f} |{:5d} |'.format(str(row['Rownum']).rjust(4),row['TotalProc']/3600.,row['TotalErrors'],row['BeamformProc']/3600.,row['BeamformErrors'],row['PrepdataProc']/3600.,row['PrepdataErrors'],row['FFTProc']/3600.,row['FFTErrors'],row['AccelProc']/3600.,row['AccelErrors'],row['FoldProc']/3600.,row['FoldErrors']))
        
        if args.mode == "vprog":
            for ri, row in enumerate(rows):
                if ri%20 == 0:
                    print('Row# |  Total Jobs |Beamform Jobs|Prepdata Jobs|'+\
                            '   FFT Jobs  |  Accel Jobs |  Fold Jobs  |')
                    print('-----|-------------|-------------|-------------|'+\
                            '-------------|-------------|-------------|')
                print('{:4d} | {:5d}/{:5s} | {:5d}/{:5s} | {:5d}/{:5s} | {:5d}/{:5s} | {:5d}/{:5s} | {:5d}/{:5s} |'.\
                        format(row['Rownum'], row['TotalJobComp'], str(row['TotalJobExp']),
                               row['BeamformJobComp'], str(row['BeamformJobExp']),
                               row['PrepdataJobComp'], str(row['PrepdataJobExp']),
                               row['FFTJobComp'], str(row['FFTJobExp']),
                               row['AccelJobComp'], str(row['AccelJobExp']),
                               row['FoldJobComp'], str(row['FoldJobExp']) ))
