cd /lustre/projects/p125_astro/DATA/1166459712/pointings
j=1529
for i in $( cat /lustre/projects/p125_astro/blindsearch/grid_positions_round7.txt ); do 

q_num=`qstat -u nswainst | wc -l`
while [ $q_num -gt 100 ]; do
sleep 100
q_num=`qstat -u nswainst | wc -l`
done

echo "Mass submit r7 #${j}" | pbs_blindsearch_pipeline.py -o 1166459712 -p $i 
let j+=1

done
