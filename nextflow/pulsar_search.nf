params.obsid = null
params.dm_min = 1
params.dm_max = 250


process ddplan {
    output:
    file 'DDplan.txt' into ddplan
    
    """
    #!/usr/bin/env python

    from lfDDplan import dd_plan
    import csv
    
    output = dd_plan(150., 30.72, 3072, 0.1, $params.dm_min, $params.dm_max)
    with open("DDplan.txt", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        for o in output:
            spamwriter.writerow(o)
    """ 
}


ddlist_ch = ddplan.splitCsv()


process dedisperse {
    input:
    set val(lowdm), val(highdm), val(dmstep), val(ndms), val(timeres), val(downsamp) from ddlist_ch

    """
    echo $lowdm $highdm $dmstep $ndms $timeres $downsamp
    """
}
