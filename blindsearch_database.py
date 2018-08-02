#!/usr/bin/env python
import os, datetime, logging
import sqlite3 as lite
from optparse import OptionParser #NB zeus does not have argparse!

DB_FILE = os.environ['CMD_BS_DB_FILE']

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d
    
def database_blindsearch_start(obsid, pointing, comment):
        DB_FILE = os.environ['CMD_BS_DB_FILE']
                        
        con = lite.connect(DB_FILE)
        with con:
                cur = con.cursor()
                
                cur.execute("INSERT INTO Blindsearch(Started, Obsid, Pointing, Comment, TotalProc, TotalErrors,BeamformProc, BeamformErrors, RFIProc, RFIErrors, PrepdataProc, PrepdataErrors, FFTProc, FFTErrors, AccelProc, AccelErrors, FoldProc, FoldErrors, CandTotal, CandOverNoise, CandDect) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (datetime.datetime.now(), obsid, pointing, comment, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
                vcs_command_id = cur.lastrowid
        return vcs_command_id
    
    
def database_script_start(table, bs_id, command, arguments,nodes,dm_file_int,time=datetime.datetime.now()):
    
    con = lite.connect(DB_FILE)
    with con:
        cur = con.cursor()
        if dm_file_int == None:
            cur.execute("INSERT INTO "+table+" (BSID, Command, Arguments, Started, CPUs) VALUES(?, ?, ?, ?, ?)", (bs_id, command, arguments, time, nodes))
        else:
            cur.execute("INSERT INTO "+table+" (BSID, Command, Arguments, Started, CPUs, DMFileInt) VALUES(?, ?, ?, ?, ?, ?)", (bs_id, command, arguments, time, nodes, dm_file_int))
        row_id = cur.lastrowid
    return row_id

def database_script_stop(table, rownum, errorcode,end_time=datetime.datetime.now()):

    con = lite.connect(DB_FILE)
    with con:
        cur = con.cursor()
        cur.execute("SELECT Ended FROM "+table+" WHERE Rownum=?", (rownum,))
        ended = cur.fetchone()[0]
        if ended is not None:
            logging.warn("Overwriting existing completion time: %s" % ended)
        cur.execute("UPDATE "+table+" SET Ended=?, Exit=? WHERE Rownum=?", (end_time, errorcode, rownum))
    return


def database_mass_update(table,file_location):
    with open(file_location,'r') as csv:
        lines = csv.readlines()
        for i,l in enumerate(lines):
            l = l.split(',')
            if i % 2 == 0:
                if len(l) == 6:
                    row_num = database_script_start(table, l[3], l[1], l[2], l[4], l[5],time=l[0])
                else:
                    row_num = database_script_start(table, l[3], l[1], l[2], l[4], None,time=l[0])
            else:
                database_script_stop(table, row_num, l[1], end_time=l[0])
    return

def database_beamform_find(table,file_location, bs_id):
    time_now = datetime.datetime.now()
    #go through the batch file for info
    with open("{0}.batch".format(file_location),'r') as batch:
        lines = batch.readlines()
        for l in lines:
            if l.startswith("export OMP_NUM_THREADS"):
                nodes = l.split("=")[1]
            if l.startswith("srun"):
                command = "make_beam" #I don't think these needs to be more robust
                arguments = l.split(command)[1]
    with open("{0}.out".format(file_location),'r') as output:
        lines = output.readlines()
        find_check = False
        for l in lines:
            if "**FINISHED BEAMFORMING**" in l:
                find_check = True
                time_seconds = float(l[1:10])
                #So this may be inaccurate because now isn't when the job 
                #finished but should get the right delta
                time_then = datetime.datetime.now() - datetime.timedelta(seconds=time_seconds)
                row_num = database_script_start(table, bs_id, command, arguments, nodes, None,\
                                                time_then)
                database_script_stop(table, row_num, 0, end_time=time_now)
        if not find_check:
            #no finshed string so likely it failed:
            #TODO make this more robust to work out how long it ran before it died
            row_num = database_script_start(table, bs_id, command, arguments, nodes, None, time_now)
            database_script_stop(table, row_num, 1, end_time=time_now)
 

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
    parser.add_option("-m", "--mode", dest="mode", metavar="mode", default='v', type=str, help='This script has three modes: "vc" used to view the database commands, "vs" used to view the database scripts, "vp" view processing and error statistics, "s" used to start a record of a script on the database, "e" used to record the end time and error code of a script on the database, "p" counts errors and processing time for one id and "b" is a special mode for receiving the total time of beamforming jobs. Default mode is v')
    parser.add_option("-f", "--file_location", dest="file_location", metavar="file_location", type=str, help='mass update csv file location.')
    
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
    start_options.add_option("-a", "--argument", dest="argument", default=None, type=str, help="The arguments that script used.")
    start_options.add_option("-n", "--nodes", dest="nodes", default=None, type=int, help="The number of cpu nodes used.")
    start_options.add_option("-d", "--dm_file_int", dest="dm_file_int", default=None, type=int, help="The DM file reference eg 1 = DM_002_004.")
    
    end_options = OptionGroup(parser, 'Script End Options')
    end_options.add_option("--errorcode", dest="errorcode", default=None, type=str, help="Error code of scripts.")
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
    elif opts.mode == 'vc' or opts.mode == 'vp':
        table = 'Blindsearch'
    elif opts.mode == 'b' or opts.command == 'make_beam':
        table = 'Beamform'
        
    
    if opts.mode == "s":
        vcs_row = database_script_start(table,opts.bs_id, opts.command, opts.argument,opts.nodes,opts.dm_file_int)
        print vcs_row
    elif opts.mode == "e":
        database_script_stop(table,opts.rownum, opts.errorcode)
    elif opts.mode == 'm':
        if opts.file_location:
            file_loc = opts.file_location
        else:
            file_loc = opts.command + '_temp_database_file.csv'
        database_mass_update(table,file_loc)
    elif opts.mode == 'b':
        database_beamform_find(table,opts.file_location, opts.bs_id)
    elif opts.mode.startswith("v"):
        con = lite.connect(DB_FILE)
        con.row_factory = dict_factory
    
        query = "SELECT * FROM " + table
        

        if opts.obsid:
            query += " WHERE Arguments LIKE '%" + str(opts.obsid) + "%'"

        if opts.recent is not None:
            query += ''' WHERE Started > "%s"''' % str(datetime.datetime.now() - relativedelta(hours=opts.recent))
            logging.debug(query)
        if opts.bs_id and not opts.mode == 'vp':
            query += " WHERE BSID='" + str(opts.bs_id) + "'"
        elif opts.mode == 'vp' and opts.bs_id:
            query += " WHERE Rownum='" + str(opts.bs_id) + "'"
            
        if opts.dm_file_int:
            query += " WHERE DMFileInt='" + str(opts.dm_file_int) + "'"
        
        with con:
            cur = con.cursor()
            cur.execute(query)
            rows = cur.fetchall()

        if opts.startrow or opts.endrow:
            rows = rows[opts.startrow:]
            if opts.endrow is not None:
                rows = rows[:opts.endrow+1]
        elif not (opts.all or opts.recent):
            rows = rows[-opts.n:]
        
        
        if opts.mode == "vc": 
            print 'Row# ','Obsid       ','Pointing                      ','Started               ','Comments'
            print '--------------------------------------------------------------------------------------------------'
            for row in rows:
                print '%-5s' % (str(row['Rownum']).rjust(4)),
                print '%-12s' % (row['Obsid']),
                print '%-30s' % (row['Pointing']),
                print '%-22s' % (row['Started'][:19]),
                print row['Comment'],
                print "\n"
                
                
        if opts.mode == "vs":
            if (table =='RFI' or table == 'Prepdata'):
                print 'Row# ','BDIS ','Started               ','Ended                 ','Exit_Code','Arguments'
            else:
                print 'Row# ','BDIS ','DM_i ','Started               ','Ended                 ','Err_Code','CPUs','Arguments'
            print '--------------------------------------------------------------------------------------------------'
            for row in rows:
                #BSID INT, Command TEXT, Arguments TEXT, Started date, Ended date, Exit
                print '%-5s' % (str(row['Rownum']).rjust(4)),
                print '%-5s' % (row['BSID']),
                if not (table =='RFI' or table == 'Prepdata' or table == 'Beamform'):
                    if str(row['DMFileInt']).endswith('\n'):
                        print '%-5s' % str((row['DMFileInt']))[:-1],
                    else:
                        print '%-5s' % (row['DMFileInt']),
                print '%-22s' % (row['Started'][:19]),
                if row['Ended'] is None:
                    print '%-22s' % (row['Ended']),
                else:
                    print '%-22s' % (row['Ended'][:19]),
                if str(row['Exit']).endswith('\n'):
                    print '%-5s' % str(row['Exit'])[:-1],
                else:
                    print '%-5s' % (row['Exit']),
                print '%-5s' % (row['CPUs']),
                print row['Arguments'],
                print "\n"
                
        if opts.mode == "vp":
            print 'Row# ', 'Total proc error# ', 'Beamform proc error# ', 'RFI proc   error# ','Prep proc  error# ','FFT proc   error# ','Accel proc error# ','Fold proc  error#'
            print '--------------------------------------------------------------------------------------------------'
            for row in rows:
                #TotalProc FLOAT, TotalErrors INT, RFIProc FLOAT, RFIErrors INT, PrepdataProc FLOAT, PrepdataErrors INT, FFTProc FLOAT, FFTErrors INT, AccelProc FLOAT, AccelErrors INT, FoldProc FLOAT, FoldErrors INT,
                print '%-5s' % (str(row['Rownum']).rjust(4)),
                print '%-10s' % (row['TotalProc']),
                print '%-7s' % (row['TotalErrors']),
                print '%-10s' % (row['BeamformProc']),
                print '%-7s' % (row['BeamformErrors']),
                print '%-10s' % (row['RFIProc']),
                print '%-7s' % (row['RFIErrors']),
                print '%-10s' % (row['PrepdataProc']),
                print '%-7s' % (row['PrepdataErrors']),
                print '%-10s' % (row['FFTProc']),
                print '%-7s' % (row['FFTErrors']),
                print '%-10s' % (row['AccelProc']),
                print '%-7s' % (row['AccelErrors']),
                print '%-10s' % (row['FoldProc']),
                print '%-7s' % (row['FoldErrors'])
                #print "\n",
                
    elif opts.mode == 'p':
        #goes through jobs of that id and command type and counts errors and processing time
        query = "SELECT * FROM " + table + " WHERE BSID=" + str(opts.bs_id)
        if opts.dm_file_int:
            query += " AND DMFileInt=" + str(opts.dm_file_int)
        print query
        con = lite.connect(DB_FILE)
        con.row_factory = dict_factory
        cur = con.cursor()
        cur.execute(query)
        rows = cur.fetchall()
        
        processing = 0.
        errors = 0
        
        for row in rows:
            #print row['Ended'], row['Started']
            #processsing += 
            if not row['Ended'] == None and not row['Started'] == None:
                end_s = date_to_sec(row['Ended'])
                start_s = date_to_sec(row['Started'])
                if (end_s - start_s) >= 0.:
                    processing += (end_s - start_s)
                else:
                    print "error in processing calc"
                    print "row num: " + str(row)
                    print "End secs: " + str(end_s)
                    print "Strat secs: " + str(start_s)
                    exit()
                    
            if str(row['Exit']).endswith("\n"):
                if not str(row['Exit'])[:-1] == "0":
                        errors += 1
            else:
                if not row['Exit'] == 0:
                        errors += 1
        #TODO make a better fix
        try:
            nodes = float(rows[0]['CPUs'])
        except:
            nodes =1
        print "Processing time (s): " + str(processing)
        print "Errors number: " + str(errors)
        print "Nodes: " + str(nodes)
        print "Database Table name: " + str(table)
        
        query = "SELECT * FROM Blindsearch WHERE Rownum='" + str(opts.bs_id) + "'"
        cur.execute(query)
        row = cur.fetchall()
        
        
        tot_proc = row[0]['TotalProc']
        if tot_proc:
            new_total_proc = tot_proc + processing/3600.*nodes
        else:
            new_total_proc = 0. + processing/3600.*nodes
            
        tot_er = row[0]['TotalErrors']
        if tot_er:
            new_total_er = tot_er + errors
        else:
            new_total_er = 0 + errors
            
        job_proc = row[0][table+'Proc']
        if job_proc:
            new_job_proc = job_proc + processing/3600.*nodes
        else:
            new_job_proc = 0. + processing/3600.*nodes
        
        job_er = row[0][table+'Errors']
        if job_er:
            new_job_er = job_er + errors
        else:
            new_job_er = 0 + errors
            
        print new_total_proc, new_total_er, new_job_proc, new_job_er
        
        
        con = lite.connect(DB_FILE)
        with con:
            cur = con.cursor()
            cur.execute("UPDATE Blindsearch SET TotalProc=?, TotalErrors=?, "+table+"Proc=?, "+table+"Errors=? WHERE Rownum=?", (str(new_total_proc)[:9], new_total_er, str(new_job_proc)[:9], new_job_er, opts.bs_id))
        
