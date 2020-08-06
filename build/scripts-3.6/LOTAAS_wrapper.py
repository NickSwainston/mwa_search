#!python

import logging
import argparse
import glob
import os
logger = logging.getLogger(__name__)

def categorize_classifier_files(out_dir):
    """
    To be run in a directory with LOTAAS classifier out files. Determines which pulsars are in at least 3 of the 5 models provided and lists them in a file

    Parameters:
    -----------
    otufile_name: string
    OPTIONAL - the name of the file to write the names of the detected pulsars to

    Returns:
    --------
    None
    """

    #sort all of the classifier files into a dictionary
    class_files = glob.glob("feature_extraction_m*")
    class_file_dict = {"positive":[], "negative":[]}
    class_cand_dict = {"m1":class_file_dict, "m2":class_file_dict, "m3":class_file_dict, "m4":class_file_dict, "m5":class_file_dict}
    
    for filename in class_files:
        split_name = filename.split("_")[-1].split(".")
        model_num = split_name[0]
        det = split_name[-1]
        class_cand_dict[model_num][det].append(filename)

    #get all of the pfd files into a list
    class_file_m1 = glob.glob("feature_extraction_m1*")
    pfd_files = []
    for afile in class_file_m1:
        f = open(afile, "r")
        for line in f.readlines():
            pfd_files.append(line)
    f.close()

    #fill a dictionary with pfds and a value for how many positive IDs each pfd has
    pulsar_pfds={}
    for key in pfd_files:
        pulsar_pfds[key]=0
    for model_num in class_cand_dict.keys():
        if class_cand_dict[model_num]["positive"]:
            print(class_cand_dict[model_num]["positive"])
            f = open(class_cand_dict[model_num]["positive"][0], "r")
            for line in f.readlines():
                pulsar_pfds[line]+=1
            f.close()

    #For each pfd with >=3 positive IDs, write that pfd to 'positive' file, else write to 'negative' file
    pos_f = open(os.path.join(out_dir, "LOTAAS_positive_detections.txt"), "w+")
    neg_f = open(os.path.join(out_dir, "LOTAAS_negative_detections.txt"), "w+")
    for pfd_key in pulsar_pfds.keys():
        if pulsar_pfds[pfd_key]>=3:
            print("detected pulsar: {}".format(pfd_key))
            pos_f.write(pfd_key.split("/")[-1])
        else:
            neg_f.write(pfd_key.split("/")[-1])
    pos_f.close()
    neg_f.close()

if __name__ == '__main__':

    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR = logging.ERROR)

    #Arguments
    parser = argparse.ArgumentParser(description="A script that handles pulsar folding operations")
    parser.add_argument("--out_dir", type=str, default="./", help="The name of the output path. Default: ./")
    parser.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level. Default: INFO", choices=loglevels.keys())
    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    categorize_classifier_files(args.out_dir)  