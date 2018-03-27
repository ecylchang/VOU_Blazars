#!/usr/bin/env python

# This script parse the logfile from VOU-Blazars 
# and output the parsed fields to a JSON file.
#
# The "logfile" is simply the screen output to a file:
# $ ./bin/vou-blazars <ra> <dec> <radius> | tee "logfile"
#
# The script will read *every* directory under
# 'VOU-Blazars/Results'. In each 'Results/*'
# directory --besides VOU-Blazars regular output--
# it is expected to find the so called logfile:
#
# logfile :: "RUN*.log"
# 
# To exemplify the expected structure and output
# by this script, consider the following example 
# depicting two runs:
# 
# Results/
# |-- 17_15_m0_47_12/
# |   |-- 1_error_map.eps
# |   |-- 1_output2.csv
# |   |-- 1_sed.eps
# |   |-- 1_Sed.txt
# |   |-- all.pdf
# |   |-- candidates.ps
# |   |-- output1.csv
# |   |-- RUN_642_17_0.log
# |   `-- RX_map.ps
# `-- 21_34_m0_07_12/
#     |-- 1_error_map.eps
#     |-- 1_output2.csv
#     |-- 1_sed.eps
#     |-- 1_Sed.txt
#     |-- 2_error_map.eps
#     |-- 2_output2.csv
#     |-- 2_sed.eps
#     |-- 2_Sed.txt
#     |-- 3_error_map.eps
#     |-- 3_output2.csv
#     |-- 3_sed.eps
#     |-- 3_Sed.txt
#     |-- all.pdf
#     |-- candidates.ps
#     |-- output1.csv
#     |-- RUN_1_21_0.log
#     `-- RX_map.ps
#
# As of version "1.0" of `read_logs`, the
# output is:
#
# Reading output from 21_34_m0_07_12
# {
#     "run_label": "21_34_m0_07_12",
#     "run_label_lbl": "1_21_0",
#     "sources": {
#         "1": {
#             "ra": 21.29775,
#             "dec": -0.02356,
#             "known_blazar": false,
#             "known_3hsp": false,
#             "known_fsrq": false,
#             "deepsky_source": true
#         },
#         "2": {
#             "ra": 21.3137,
#             "dec": -0.10908,
#             "known_blazar": false,
#             "known_3hsp": false,
#             "known_fsrq": false,
#             "deepsky_source": false
#         },
#         "3": {
#             "ra": 21.37021,
#             "dec": -0.09889,
#             "known_blazar": true,
#             "known_3hsp": false,
#             "known_fsrq": true,
#             "deepsky_source": true
#         }
#     }
# }
# 
# 
# Reading output from 17_15_m0_47_12
# {
#     "run_label": "17_15_m0_47_12",
#     "run_label_lbl": "642_17_0",
#     "sources": {
#         "1": {
#             "ra": 17.11182,
#             "dec": -0.62339,
#             "known_blazar": true,
#             "known_3hsp": false,
#             "known_fsrq": true,
#             "deepsky_source": true
#         }
#     }
# }
# 
# 
# Known blazars:  2
# NEW Blazars:  2
#
# ================================================

# Glob/pattern name for run logfiles:
#
LOGFILE_GLOB='RUN*.log'

# Output file with parsed content:
#
OUTPUT_JSON = 'results.json'

import os
import glob
import string
import json
import re

vou_result_run_dirs = glob.glob('[0-9]*')

cnt_known_blazars = 0
cnt_new_blazars = 0
results = {}

for vou_dir in vou_result_run_dirs:

    if not os.path.isfile(os.path.join(vou_dir, '1_sed.eps')):
        continue

    run_label = vou_dir
    print('Reading output from {}'.format(run_label))

    result_run = {}
    result_run['run_label'] = run_label

    # Read the Logfile
    #
    logfile = glob.glob(os.path.join(vou_dir, LOGFILE_GLOB))
    assert len(logfile) == 1
    logfile = logfile[0]
    
    with open(logfile, 'r') as fp:

        sources = {}
        src = 0
        known_blazar = False
        known_3hsp = False
        known_fsrq = False

        for line in fp: #.readlines():
            line = line.strip()

            # There is a summary block at the first part of the
            # logfile with information about blazar candidates
            # regarding Radio and X-Ray emition.
            # We use this first block to count the number of
            # candidates we'll deal with.
            #
            if src == 0:
                match_source = re.match('.*Match nr.(.*)', line)
                if match_source:
                    ll = match_source.group(1).replace(',', '')
                    ll = ll.split()
                    sources[int(ll[0])] = {} #list(map(float, ll[-2:]))
                    continue
    
                match_candidate = re.match('.*Candidate nr.(.*)', line)
                if match_candidate:
                    ll = match_candidate.group(1).replace(',', '')
                    ll = ll.split()
                    sources[int(ll[0])] = {}
                    continue

            # Each source's block begins with '=====...',
            # when the first or a new block begins we
            # have to flush any information from the last
            # source and increment the source-id
            #
            if re.match('.*=================.*', line):
                if src in sources:
                    sources[src].update({'known_blazar':known_blazar,
                                         'known_3hsp':known_3hsp,
                                         'known_fsrq':known_fsrq})
                    known_blazar = False
                    known_3hsp = False
                    known_fsrq = False
                src += 1

            if src == 0 or src not in sources:
                continue

            # At the beginning of each source's block, get the
            # respective RA and Dec 
            #
            if src > 0:
                match_position = re.match('.*R\.A.*Dec.*=(.*)', line)
                if match_position:
                    ll = match_position.group(1).replace('\s','').split(',')
                    ll = list(map(float, ll))
                    sources[src].update({'ra':ll[0], 'dec':ll[1]})
                    continue

            # Here we check the messages at the end of each source's
            # block which tell whether the sources is a known blazar
            #
            if src > 0:
                match_known = re.match('Known blazar', line)
                if match_known:
                    known_blazar = True
                match_3hsp = re.match('Already in 3HSP', line)
                if match_3hsp:
                    known_3hsp = True
                match_fsrq = re.match('Flat radio spectrum source', line)
                if match_fsrq:
                    known_fsrq = True
                continue


        if src in sources:
            sources[src].update({'known_blazar':known_blazar,
                                 'known_3hsp':known_3hsp,
                                 'known_fsrq':known_fsrq})
        known_blazars = sum(dct['known_blazar'] for src,dct in sources.items())
        cnt_known_blazars += known_blazars
        cnt_new_blazars += len(sources) - known_blazars

        result_run['sources'] = sources

    # Read each source output (where 'X' is the source id "srv":
    # * X_output2.csv
    # * X_Sed.txt
    #
    for src,dct in result_run['sources'].items():

        sed_file = '{:d}_Sed.txt'.format(int(src))
        sed_file = glob.glob(os.path.join(vou_dir, sed_file))[0]
        with open(sed_file, 'r') as fp:
            sed_content = fp.read()

        match_deepsky = re.search('XRTDEEP', sed_content)
        deepsky_source = bool(match_deepsky)
        dct.update({'deepsky_source':deepsky_source})

        dct.update({'sed_file':sed_file})

    results[run_label] = result_run['sources']
    print(json.dumps(result_run, indent=4))
    print('\n')

with open(OUTPUT_JSON,'w') as fp:
    json.dump(results, fp, indent=4)

print('Known blazars: ', cnt_known_blazars)
print('NEW Blazars: ', cnt_new_blazars)

