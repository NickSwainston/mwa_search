.. _smart_processing:

Test Commands
=============

A bunch of command to run for testing using the first 600 seconds of 1301674968


beamform.nf
-----------

normal::

    beamform.nf --obsid 1301674968 --calid 1301739904 --begin 1301674969 --end 1301675568 --pointings 13:11:52.64_-12:28:01.63,14:18:50.28_-39:21:18.51 -w test_work --out_dir test_cands --vcstools_version devel --publish_fits

summed::

    beamform.nf --obsid 1301674968 --calid 1301739904 --begin 1301674969 --end 1301675568 --pointings 13:11:52.64_-12:28:01.63,14:18:50.28_-39:21:18.51 -w test_work --out_dir test_cands --vcstools_version devel --publish_fits --summed

ipfb::

    beamform.nf --obsid 1301674968 --calid 1301739904 --begin 1301674969 --end 1301675568 --pointings 13:11:52.64_-12:28:01.63,14:18:50.28_-39:21:18.51 -w test_work --out_dir test_cands --vcstools_version devel --publish_fits --ipfb


pulsar_search.nf
----------------
Run with outputs of beamform.nf

simple periodic::

    pulsar_search.nf --obsid 1301674968 --calid 1301739904 --fits_file /fred/oz125/vcs/1301674968/pointings/13:11:52.64_-12:28:01.63/*fits -w test_work --out_dir test_cands  --vcstools_version devel --dm_min 36 --dm_max 37

Just single pulse search::

    pulsar_search.nf --obsid 1301674968 --calid 1301739904 --fits_file /fred/oz125/vcs/1301674968/pointings/13:11:52.64_-12:28:01.63/*fits -w test_work --out_dir test_cands  --vcstools_version devel --dm_min 36 --dm_max 37 --sp