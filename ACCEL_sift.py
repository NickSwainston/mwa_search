#! /usr/bin/env python
import sifting, re, glob, argparse

def dm_i_to_file(dm_i):
    if 0 <= dm_i <= 3:
        dm_file = "DM_00" + str(dm_i*2) + "-00" + str((dm_i+1)*2)
    elif dm_i == 4:
        dm_file = "DM_008-010"
    elif 5 <= dm_i <= 48:
        dm_file = "DM_0" + str(dm_i*2) + "-0" + str((dm_i+1)*2)
    elif dm_i == 49:
        dm_file = "DM_098-100"
    elif 50 <= dm_i <= 149:
        dm_file = "DM_" + str(dm_i*2) + "-" + str((dm_i+1)*2)
    else:
        print dm_i
    return dm_file
    
parser = argparse.ArgumentParser(description="""
Finds candidates with DM paterns.
""")
parser.add_argument('dir',type=str,help="Work dir. Should be my blindsearch sub dirs, eg DM_000-002.")
args=parser.parse_args()

# Note:  You will almost certainly want to adjust
#        the following variables for your particular search
 
d=args.dir
print d
if d.endswith("/"):
    d = d[:-1]
if d[-10:-7] == "DM_":
    i = int(d[-7:-4])
    j = int(d[-3:])
    # glob for ACCEL files
    globaccel = d+"/*ACCEL_*0"
    # glob for .inf files
    globinf = d+"/*DM*.inf"
    inffiles = glob.glob(globinf)
    candfiles = glob.glob(globaccel)
    if i > 0:
        #search for append the previous few DM steps
        dm_i = i/2 -1
        dm_file = dm_i_to_file(dm_i)
        globaccel_append = d[:-10] + dm_file + "/*DM"+str(i-1)+".9*ACCEL_*0"
        globinf_append = d[:-10] + dm_file + "/*DM"+str(i-1)+".9*.inf"
        #print d[:-10] + dm_file + "/*DM"+str(i-1)+".9*.inf"
        
        #append this extra files
        inftemp = glob.glob(globinf_append)
        for i in inftemp:
            inffiles.append(i)
        candtemp = glob.glob(globaccel_append)
        for i in candtemp:
            candfiles.append(i)
else:
    # glob for ACCEL files
    globaccel = d+"*ACCEL_*0"
    # glob for .inf files
    #globinf = "../*/*DM*.inf"
    globinf = d+"*DM*.inf"
    inffiles = glob.glob(globinf)
    candfiles = glob.glob(globaccel)
#print candfiles
#print inffiles

# In how many DMs must a candidate be detected to be considered "good"
min_num_DMs = 2
# Lowest DM to consider as a "real" pulsar
low_DM_cutoff = 2.0
# Ignore candidates with a sigma (from incoherent power summation) less than this
sifting.sigma_threshold = 4.0
# Ignore candidates with a coherent power less than this
sifting.c_pow_threshold = 100.0

# If the birds file works well, the following shouldn't
# be needed at all...  If they are, add tuples with the bad
# values and their errors.
#                (ms, err)
sifting.known_birds_p = []
#                (Hz, err)
sifting.known_birds_f = []

# The following are all defined in the sifting module.
# But if we want to override them, uncomment and do it here.
# You shouldn't need to adjust them for most searches, though.

# How close a candidate has to be to another candidate to                
# consider it the same candidate (in Fourier bins)
sifting.r_err = 1.1
# Shortest period candidates to consider (s)
sifting.short_period = 0.0005
# Longest period candidates to consider (s)
sifting.long_period = 15.0
# Ignore any candidates where at least one harmonic does exceed this power
sifting.harm_pow_cutoff = 8.0

#--------------------------------------------------------------

# Try to read the .inf files first, as _if_ they are present, all of
# them should be there.  (if no candidates are found by accelsearch
# we get no ACCEL files...

# Check to see if this is from a short search
print inffiles
if len(re.findall("_[0-9][0-9][0-9]M_" , inffiles[0])):
    dmstrs = [x.split("DM")[-1].split("_")[0] for x in candfiles]
else:
    dmstrs = [x.split("DM")[-1].split(".inf")[0] for x in inffiles]
dms = map(float, dmstrs)
dms.sort()
dmstrs = ["%.2f"%x for x in dms]

# Read in all the candidates
cands = sifting.read_candidates(candfiles)

# Remove candidates that are duplicated in other ACCEL files
if len(cands):
    cands = sifting.remove_duplicate_candidates(cands)

# Remove candidates with DM problems
if len(cands):
    cands = sifting.remove_DM_problems(cands, min_num_DMs, dmstrs, low_DM_cutoff)

# Remove candidates that are harmonically related to each other
# Note:  this includes only a small set of harmonics
if len(cands):
    cands = sifting.remove_harmonics(cands)

# Write candidates to STDOUT
if len(cands):
    cands.sort(sifting.cmp_sigma)
    if d[-10:-7] == "DM_":
        dm_file = dm_i_to_file(int(d[-7:-4])/2)
        print "Finished for foulder: " + dm_file
        sifting.write_candlist(cands,'cand_files/cands_'+ d.split("/")[0]+"_"+dm_file+'.txt')
        print 'cand_files/cands_'+ d.split("/")[0]+"_"+dm_file+'.txt'
    else:
        sifting.write_candlist(cands,'ACCEL_sift_cands.txt')

