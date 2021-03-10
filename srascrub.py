#!/usr/bin/python3

"""
    srascrub utility - version: 0.01
    
    Author:         Matthew DePugh
    Organization:   Brown Lab, Dept. of Biology, Texas Tech University
    Contact:        matthew.depugh@ttu.edu
    
    Purpose:        The purpose of this software is to provide a command line utility that searches the NCBI sra databases,
                    with a user defined search term, and returns human readable results.
                    
    Goal:           Fit this into a larger automated utility to eventually download SRRs, decompress, and perform BLAST.
                    
    Dependencies:   Python Python 3.9.1
                    NCBI SRA Toolkit 2.10.9
                    
    Execution:
                    [prompt]$ python3 srascrub.py -e user@email.com -d /Use/complete/file/path/to/write/outfile -t "this is my search term"
"""

import getopt
import os
import re
import sys
import xml.etree.ElementTree as ET

from Bio import Entrez

"""
    sra_search - Queries NCBI SRA database with search term 'sterm'. Returns a list of SRA Study (SRX) IDs.
    
    Arguments:
        email: Required by NCBI, should be valid format *@*.*
        sterm: The user defined search term.
        resmax: Should be an integer that signifies the max number of results accepted.
    
    error checking
        Arguments: No. Should be verified as correct in caller.
        
        Entrez: Will attempt to connect to NCBI handle. Will complete function and return the "IdList" list if connection is established.
                Returns 71000 on connection failure.
"""
def sra_search(email, sterm, resmax):
    try:
        Entrez.email = email
        handle = Entrez.esearch(db="sra", term=sterm, retmax=resmax)
        record = Entrez.read(handle)
    except:
        print("Connectivity error: Could not establish connection with NCBI database")
        sys.exit(71000)
    
    handle.close()
    return record["IdList"]

"""
    sra_summaries - Uses the list of NCBI Study IDs (SRX values) to query the NCBI database and return corresponding XMLs that detail the
                    experiments within the study. This information is then parsed, added to a dictionary datastructure, and returned to
                    calling function.
    
    Arguments:
        email: Required by NCBI, should be valid format *@*.*
        idList: List of NCBI SRA study ID values (SRX).
        resmax: Should be an integer that signifies the max number of results accepted.

    error checking
        Entrez: Will attempt to connect to NCBI handle. Will complete function and return the "IdList" list if connection is established.
                Returns 71000 on connection failure.
"""
def sra_summaries(email, idList, resmax):
    study = {}
    experiment = {}
    # this is the return variable which will be populated with return
    retvar = {}

    for i in idList:
        # experiment dictionary
        #summary dictionary
        addDictS = {"study": "", "experiments": []}
        addDictE = {"name": "", "instrument": "", "organism": "", "submitter": {"subId": "", "name": "", "labname": ""}, "taxid": ""}
        
        try:
            Entrez.email = email
            handle = Entrez.esummary(db="sra", id=i)
            record = Entrez.read(handle)
        except:
            print("Connectivity error: Could not establish connection with NCBI database")
            sys.exit(71000)
        
        
        #variable text must be enclosed within 1 "root" tag
        tempstr0 = "<?xml version=\"1.0\"?><root>" + record[0]['ExpXml'] + "</root>"
        root = ET.fromstring(tempstr0)
        addKeyStudy = root[3].attrib['acc']
        addKeyExp = root[2].attrib['acc']
        #print(record[0]['ExpXml'])
        if addKeyStudy not in study.keys():
            addDictS["study"] = root[3].attrib["name"]
            addDictS["experiments"].append(root[2].attrib["acc"])
            study[addKeyStudy] = addDictS
        else:
            study[addKeyStudy]["experiments"].append(root[2].attrib["acc"])
        
        if addKeyExp not in experiment.keys(): # potentially unnecessary if condition ????
            addDictE["name"] = root[2].attrib["name"]
            addDictE["instrument"] = root[0][1].attrib["instrument_model"]
            # fix this
            if "ScientificName" in root[4].attrib:
                addDictE["organism"] = root[4].attrib["ScientificName"]
            addDictE["submitter"]["subId"] = root[1].attrib["acc"]
            addDictE["submitter"]["name"] = root[1].attrib["contact_name"]
            addDictE["submitter"]["labname"] = root[1].attrib["lab_name"]
            addDictE["taxid"] = root[4].attrib["taxid"]
            experiment[addKeyExp] = addDictE
        
        handle.close()
    
    retvar["studies"] = study
    retvar["experiments"] = experiment
    return retvar

"""
    human_read_file - This function accepts the experiment and study dictionaries and writes the contents to an output file in an
                      iterative manner.
                      
    Arguments:
        file_path: This is the name of the file (as a complete file path).
        experiments: This is a data structure holding "experiment" data.
        studies: This is a data structure holding "study" data.
        
    Error Checking: None
"""
def human_read_file(file_path, experiment, study):
    file = open(file_path, "w")
    for i in study.keys():
        file.write("STUDY: " + i + ": " + study[i]["study"] + "\n")
        file.write("Experiments: \n")
        for j in study[i]["experiments"]:
            file.write("\t\t" + j + "\n")
            file.write("\t\t\t\tName: " + experiment[j]["name"] + "\n")
            file.write("\t\t\t\tInstrument : " + experiment[j]["instrument"] + "\n")
            file.write("\t\t\t\tOrganism: " + experiment[j]["organism"] + "\n")
            file.write("\t\t\t\tTax ID: " + experiment[j]["taxid"] + "\n")
            file.write("\t\t\t\tSUBMITTER: " + experiment[j]["submitter"]["subId"] + ": " + experiment[j]["submitter"]["name"] + ", " + experiment[j]["submitter"]["labname"] + "\n")
        file.write("\n\n")
    file.close()

def csv_file_write(file_path, experiment, study):
    file = open(file_path, "w")
    for i in study.keys():
        for j in study[i]["experiments"]:
            file.write(i + ", " + study[i]["study"] + ", " + j + ", " + experiment[j]["name"] + ", " + experiment[j]["instrument"] + ", " + experiment[j]["organism"] + ", " + experiment[j]["taxid"] + ", " + experiment[j]["submitter"]["subId"] + ", " + experiment[j]["submitter"]["name"] + ", " + experiment[j]["submitter"]["labname"] + "\n")
    file.close()
    
"""
    "These aren't the droids we're looking for... you can go about your business... move along"
"""
def DEBUG__PRINT_STUDIES(study, experiment) :
    for i in study.keys():
        print("STUDY: " + i + ": " + study[i]["study"])
        print("Experiments: ")
        for j in study[i]["experiments"]:
            print("\t\t" + j)
            print("\t\t\t\tName: " + experiment[j]["name"])
            print("\t\t\t\tInstrument : " + experiment[j]["instrument"])
            print("\t\t\t\tOrganism: " + experiment[j]["organism"])
            print("\t\t\t\tTax ID: " + experiment[j]["taxid"])
            print("\t\t\t\tSUBMITTER: " + experiment[j]["submitter"]["subId"] + ": " + experiment[j]["submitter"]["name"] + ", " + experiment[j]["submitter"]["labname"])
        print("\n\n")

# main system script --- includes argument parser

"""
    main - The key take away from main is that parsed XML information obtained by the sra_summaries function are stored in a data
           structure named "summary_results". This variable is passed to human_read_file which writes the contents of the
           "summary_results" to the outfile.
    
    Arguments: Passed from command line.
        -help: Prints help menu.
        -dir: User specified directory. Expects full path.
        -email: User specified email.
        -term: Search term, if space character is used then use "" characters
        
    Error Checking
        Command Line Variables: These are checked thoroughly to ensure filepaths exist and the email is of a valid and expected form.
"""
def main(argv):
    
    csvflag = 0
    #ncbi user email
    email = ""
    emailflag = 0
    # search term
    sterm = ""
    stermflag = 0
    # output directory
    dir = ""
    dirflag = 0
    # max results
    resmax = 40000
    
    csv_file = ""
    human_file_path = dir + "human_readable.txt"
    
    # argument parsing start
    try:
        opts, args = getopt.getopt(argv, "hc:d:e:H:t:", ["help", "csvfile", "dir", "email", "hfile", "term"])
    except getopt.GetoptError:
        print("Argument Input Error 55: Option Declaration Error")
        print("srascrub.py example:")
        print("\t python3 srascrub.py -e myemail@email.com -d /full/path/to/directory -t \"this is my search term\"")
        sys.exit(55)
    for opt, arg in opts:
        if opt == '-h' or opt == '--help':
            print ("srascrub.py arguments:\n\t-h, --help : print help menu\n\t-d, --dir : directory for results\n\t-e, --email : ncbi user email\n\t-t, --term : search term")
            sys.exit(104)
        elif opt in ("-d", "--dir"):
            if re.search(r'^/', arg) and os.path.isdir(arg):
                    if (re.search("/\Z", arg)):
                        dir = arg
                    else:
                        dir = arg + "/"
                    dirflag = 1
            else:
                print("Argument Input Error 100: Invalid File Path")
                print("Must begin with / and contain the full file path")
                sys.exit(100)
        elif opt in ("-e", "--email"):
            if re.match(r'(.+)@(.+)\..+', arg):
                email = arg
                emailflag = 1
            else :
                print("Argument Input Error 101: Invalid Email")
                print("emails must@havecharacters.here")
                sys.exit(101)
        elif opt in ("-t", "--term"):
            sterm = arg
            stermflag = 1
 
        elif opt in ("-H", "--hfile"):
            human_file_path = arg
        
        elif opt in ("-c", "csvfile"):
            if re.match(r'(.+)\.csv', arg):
                csvflag = 1
                csv_file = arg
            else:
                print("Argument Input Error 99: Invalid CSV File Name")
                print("CSV File Names must end in .csv")
                sys.exit(99)
 
    if stermflag == 0 or emailflag == 0 or dirflag == 0:
        print("Argument Input Error 55555: Require Flag(s) Missing")
        print("Required Flags are: -d, -t, -e")
        sys.exit(55555)

    # outfile definitions assigned here
    # human readable file path
    #human_file_path = dir + "human_readable.txt"
    # studies out file
    study_out_file_path = dir + "study"
    # experiment out file
    exp_out_file_path = dir + "exp"
    
    
    idList = sra_search(email, sterm, resmax)
    summary_results = sra_summaries(email, idList, resmax)
    human_read_file(human_file_path, summary_results["experiments"], summary_results["studies"])
    if csvflag != 0:
        csv_file_write(csv_file, summary_results["experiments"], summary_results["studies"])
    sys.exit(0)
    
if __name__ == "__main__":
   main(sys.argv[1:])
