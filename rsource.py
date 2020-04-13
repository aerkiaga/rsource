#!/usr/bin/python3

import os, sys, shutil, configparser

GRCh_revision = "38"
chromosomes = ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'mt')

conf = {
    'delete sequence' : True,
    'delete annotations' : True,
    'delete gaps' : False
}

script_path = os.path.realpath(__file__)
path = os.path.dirname(script_path)

python3_path = None
pypy3_path = None
def get_python_paths():
    global python3_path, pypy3_path
    python3_path = sys.executable
    python3_path = shutil.which("python3")
    if not python3_path:
        python3_path = shutil.which("python")
        s = os.popen(python3_path + " -V", 'r').read()
        if s == "": # < 3.4
            python3_path = None
    if not python3_path:
        for v in range(8, -1, -1):
            python3_path = shutil.which("python3.{}".format(v))
            print(python3_path)
            if python3_path:
                break
    pypy3_path = shutil.which("pypy3")
    if not pypy3_path:
        if not python3_path:
            print("Python 3: not found")
        else:
            print("Python 3: " + python3_path)
            pypy3_path = python3_path
        print("PyPy3: not found")
    else:
        if not python3_path:
            python3_path = pypy3_path
        print("Python 3: " + python3_path)
        print("PyPy3: " + pypy3_path)
    if python3_path:
        python3_path = "\"" + python3_path + "\""
    if pypy3_path:
        pypy3_path = "\"" + pypy3_path + "\""

def get_sequence(sequence_gz_path):
    sequence_url = "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh" + GRCh_revision + "_latest/refseq_identifiers/GRCh" + GRCh_revision + "_latest_genomic.fna.gz"
    command = "wget -O \"" + sequence_gz_path + "\" \"" + sequence_url + "\""
    print("Command: " + command)
    os.system(command)

def get_annotations(annotations_gz_path):
    annotations_url = "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh" + GRCh_revision + "_latest/refseq_identifiers/GRCh" + GRCh_revision + "_latest_genomic.gff.gz"
    command = "wget -O \"" + annotations_gz_path + "\" \"" + annotations_url + "\""
    print("Command: " + command)
    os.system(command)

def condense(sequence_gz_path):
    condense_script_path = os.path.join(path, "condense.py")
    command = "gunzip -dc \"" + sequence_gz_path + "\" | " + pypy3_path + " \"" + condense_script_path + "\""
    print("Command: " + command)
    os.system(command)

def comment(annotations_gz_path):
    comment_script_path = os.path.join(path, "comment.py")
    command = "gunzip -dc \"" + annotations_gz_path + "\" | " + pypy3_path + " \"" + comment_script_path + "\""
    print("Command: " + command)
    os.system(command)

def check_exist(extension):
    for ch in chromosomes:
        file_path = os.path.join(path, ch + extension)
        if not os.path.isfile(file_path):
            return False
    return True

def make_chromosomes():
    sequence_gz_path = os.path.join(path, "sequence.fna.gz")
    if not os.path.isfile(sequence_gz_path):
        print("No sequence found; downloading...")
        get_sequence(sequence_gz_path)
    else:
        print("Sequence found!")
    print("Packing chromosomes...")
    condense(sequence_gz_path)
    if conf['delete sequence']:
        command = "rm \"" + sequence_gz_path + "\""
        print("Command: " + command)
        os.system(command)

def make_metadata():
    annotations_gz_path = os.path.join(path, "annotations.gff.gz")
    if not os.path.isfile(annotations_gz_path):
        print("No annotations found; downloading...")
        get_annotations(annotations_gz_path)
    else:
        print("Annotations found!")
    print("Commenting...")
    comment(annotations_gz_path)
    if conf['delete annotations']:
        command = "rm \"" + annotations_gz_path + "\""
        print("Command: " + command)
        os.system(command)
    if conf['delete gaps']:
        for ch in chromosomes:
            gap_path = os.path.join(path, ch + ".gap")
            if os.path.isfile(gap_path):
                command = "rm \"" + gap_path + "\""
                print("Command: " + command)
                os.system(command)

def get_config(section, config):
    conf[config] = section.getboolean(config, conf[config])

def parse_config():
    config = configparser.ConfigParser()
    config_path = os.path.join(path, "config.ini")
    config.read(config_path)
    if 'Setup' in config:
        section = config['Setup']
        get_config(section, 'delete sequence')
        get_config(section, 'delete annotations')
        get_config(section, 'delete gaps')

get_python_paths()
parse_config()
chromosomes_exist = check_exist(".bin")
if not chromosomes_exist:
    print("No chromosomes present")
    make_chromosomes()
metadata_exist = check_exist(".dat")
gaps_exist = check_exist(".gap")
if not metadata_exist:
    print("No metadata present")
    if not gaps_exist:
        print("No gap annotations found; generating from sequence...")
        make_chromosomes()
    make_metadata()

if chromosomes_exist and metadata_exist:
    viewer_script_path = os.path.join(path, "viewer.py")
    args = " ".join(sys.argv[1:])
    command = python3_path + " \"" + viewer_script_path + "\" " + args
    print("Command: " + command)
    os.system(command)
else:
    print("Exiting. Run again to view.")
