.. _candidate_follow_up:

SMART Candidate Follow Up
=========================

The following guide will teach you the different stages of following up a candidate


Following Up Multiple candidates
--------------------------------

Once most of an observation has been classified on the `SMART Pulsar Web Classifier<https://apps.datacentral.org.au/smart/>`_, you can claim/download the promising candidates for follow-up with the full 80 minutes of data.

Put all of the candidates' bestprof files into a directory on Garrawarla (make sure they are from the same obs ID), make sure the combined and calibration are downloaded and then run the command::

    mwa_search_pipeline.nf --obsid <obsid> --calid <calid> --all --bestprof_pointings <directory_with_bestprofs>

Once the processing is finished, the full-length obs candidates will be output to::

    /astro/mwavcs/pulsar_search/<obsid>_candidates

For each promising full-length obs candidate, perform the following steps.

Finding the Best SMART Positon
------------------------------

You can improve the full-length SMART detection by creating a grid of pointings around the candidate and try and find the highest signal-to-noise position.
There is a pipeline that will already do this for you, and you can run it with the following command::

    find_candidate_position.nf --obsid <obsid> --calid <calid> --all --pointing_gird <cand_pos HH:MM:SS_+DD:MM:SS> --fraction 0.6 --loops 4 --period <period in s> --dm <dm> --summed=false

This will make 4 loops/rings of pointings which is 61 pointings which may sound like a lot but will make sure the candidate was not a grating/side lobe detection.
The output of the position improvement will be put in::

    /astro/mwavcs/pulsar_search/<obsid>_candidate_follow_up

The predicted position is output in `predicted_pos.txt`, and you can see a heatmap of the detection SN in `residual.png`.
From this, you can assume the position accuracy is half a beam width (~9 arcmins).

Finding Archival Observations for Follow Up
-------------------------------------------

To improve the position accuracy further, you can use the MWA VCS archive to find Phase 1 and Phase 2 extended array observations, which have better resolution.
The command to all observations that contain your candidate you can use the following command::

    find_pulsar_in_obs.py --obs_for_source -c "<cand_pos HH:MM:SS_+DD:MM:SS>"

The output file will contain something like this::

    #Obs ID   |Dur |Enter|Exit |Power| OAP | Freq | Band
    1148652032 3600 0.000 0.509 0.808  P1  184.96  30.72
    1215080360 3600 0.213 1.000 0.524  P2E  215.68  30.72
    1301674968 4800 0.000 1.000 0.993  P2C  154.24  30.72
    1302712864 4800 0.000 0.275 0.418  P2C  154.24  30.72

From this information, you can decide which observations to follow up on first based on how likely you are to detect the pulsar (observation length and maximum zenith normalised power).
It is easiest to follow up on Phase 1 observations first as they have a larger resolution, so that it will take fewer beams.

Finding the Best Position in a Follow-Up Observation
----------------------------------------------------

Different observations will have different beam widths (FWHM), and you must take this into account when following up.
You can use the following command to estimate your follow-up observation's FWHM::

    mwa_metadb_utils.py <obsid> | grep FWHM

You can use this to work out how many loops of beams you need for the follow-up::

    loops = <position uncertainty radius> * 0.6 / <obs FWHM>

Round up this number to the nearest integer and use it in the following command::

    find_candidate_position.nf --obsid <obsid> --calid <calid> --all --pointing_gird <cand_pos HH:MM:SS_+DD:MM:SS> --fraction 0.6 --loops <loops> --period <period in s> --dm <dm> --summed=false
