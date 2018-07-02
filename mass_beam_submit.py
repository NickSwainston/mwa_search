import subprocess, os, time
all_files = os.listdir("/lustre/projects/p125_astro/DATA/1166459712/pointings/")
#print all_files


old = []
old_round_files = ["grid_positions_dec_limited_r1.txt","grid_positions_dec_limited_r2.txt",\
                   "grid_positions_dec_limited_r3.txt","grid_positions_dec_limited_r4.txt",\
                   "grid_positions_dec_limited_r5.txt","grid_positions_dec_limited_r6.txt"]
for o in old_round_files:
    with open('/lustre/projects/p125_astro/blindsearch/'+o,'r') as lines:
        for l in lines:
            if not l.startswith('#'):
                ra, dec, temp, temp, temp, temp = l.split()
                if not ra.startswith("0"):
                    ra ="0" + ra
                if len(ra) > 11:
                    ra = ra[:11]
                if len(dec) > 12:
                    dec = dec[:12]
                    
                if len(ra) == 8:
                    ra = ra + '.00'
                if len(dec) == 9:
                    dec = dec + '.00'
                pointing = ra + "_" + dec
                old.append(pointing)
    
        
new = []
with open('/lustre/projects/p125_astro/blindsearch/grid_positions_dec_limited_r7.txt','r') as lines:
    for l in lines:
        if not l.startswith('#'):
            ra, dec, temp, temp, temp, temp = l.split()
            if not ra.startswith("0"):
                ra ="0" + ra
            if len(ra) > 11:
                ra = ra[:11]
            if len(dec) > 12:
                dec = dec[:12]
                
            if len(ra) == 8:
                ra = ra + '.00'
            if len(dec) == 9:
                dec = dec + '.00'
            pointing = ra + "_" + dec
            new.append(pointing)
  
with open('/lustre/projects/p125_astro/blindsearch/grid_positions_round7.txt','w') as out:     
    for l in new:
        if l not in old:
            out.write(l+'\n')



sub_count = 0
with open('/lustre/projects/p125_astro/blindsearch/grid_positions_round7.txt','r') as lines:
    for i, line in enumerate(lines):
        if not line in all_files:
            ra, dec = line.split("_")
            print "ra dec: ",ra, dec
            print "i : ",i
            
            submit_line = 'qstat -u nswainst | wc -l' 
            submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
            q_num = ""
            for line in submit_cmd.stdout:
                    q_num += line
            print "q: " + str(int(q_num))
            while (int(q_num) > 100 ):
                print "waiting 100 s for queue to clear"
                time.sleep(100)
                submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
                q_num = ""
                for line in submit_cmd.stdout:
                        q_num += line
            
            submit_line = 'python /home/nswainst/vcstools/galaxy-scripts/scripts/process_vcs.py -m beamform -b 1166459719 -e 1166460019 -o 1166459712 --DI_dir=/lustre/projects/p125_astro/DATA/1166440392/ --flagged_tiles=/lustre/projects/p125_astro/DATA/1166440392/flagged_tiles.txt -w /lustre/projects/p125_astro/DATA -p ' + ra + ' ' + dec[:-1] + ' --bf_new --PBS'
            print submit_line
            submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
            for line in submit_cmd.stdout:
                    print line,
            """
            sub_count +=1
            if sub_count%10 == 0:
                print "Waiting 60 minutes before sending off more jobs"
                time.sleep(3600)
            """
submit_line = 'bash /lustre/projects/p125_astro/DATA/1166459712/mass_split.batch'
submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
#07:35:05.58_-28:26:43.88
