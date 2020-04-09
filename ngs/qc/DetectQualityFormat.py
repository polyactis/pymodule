#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Guess the encoding of a stream of qual lines.

main function return list of possible quality format or 
return ['chenlin', gmin, gmax] if no encodings for range

Detect quality format of files with information in database. 
list to be detected is from table "individual_alignment.ind_seq_id"
so, some id in table "individual_sequence" is ignored.

"""

from __future__ import with_statement, division, print_function

import fileinput
import operator
import optparse
import sys, os

from collections import Counter
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.expanduser('~/script'))
sys.path.insert(0, os.path.expanduser('~/src'))

#  Note that the theoretical maximum for all encodings is 126.
#  The upper limits below are for "typical" data only.
RANGES = {
    'Sanger': (33, 73),
    'Illumina-1.8': (33, 74),
    'Solexa': (59, 104),
    'Illumina-1.3': (64, 104),
    'Illumina-1.5': (66, 105)
}

# The threshold to decide between Illumina-1.3 and Illumina-1.5
# based upon how common "B" is. The threshold insists it is
# within the Nth most common quality scores.
# N.B. needs to be conservative, as this is applied per input line.
# N_MOST_COMMON_THRESH = 4


def get_qual_range(qual_str):
    """
    >>> get_qual_range("DLXYXXRXWYYTPMLUUQWTXTRSXSWMDMTRNDNSMJFJFFRMV")
    (68, 89...)
    """

    qual_val_counts = Counter(ord(qual_char) for qual_char in qual_str)

    min_base_qual = min(qual_val_counts.keys())
    max_base_qual = max(qual_val_counts.keys())

    return (min_base_qual, max_base_qual, qual_val_counts)


def get_encodings_in_range(rmin, rmax, ranges=RANGES):
    valid_encodings = []
    for encoding, (emin, emax) in ranges.items():
        if rmin >= emin and rmax <= emax:
            valid_encodings.append(encoding)
    return valid_encodings


# def heuristic_filter(valid, qual_val_counts):
#     """Apply heuristics to particular ASCII value scores
#        to try to narrow-down the encoding, beyond min/max.
#     """

#     if 'Illumina-1.5' in valid:
#         # 64–65: Phread+64 quality scores of 0–1 ('@'–'A')
#         #        unused in Illumina 1.5+
#         if qual_val_counts[64] > 0 or qual_val_counts[65] > 0:
#             valid.remove('Illumina-1.5')

#         # 66: Phread+64 quality score of 2 'B'
#         #     used by Illumina 1.5+ as QC indicator
#         elif 66 in map(operator.itemgetter(0),
#                        qual_val_counts.most_common(N_MOST_COMMON_THRESH)):
#             print("# A large number of 'B' quality scores (value 2, ASCII 66) "
#                   "were detected, which makes it likely that this encoding is "
#                   "Illumina-1.5, which has been returned as the only option.",
#                   file=sys.stderr)
#             valid = ['Illumina-1.5']

#     return valid


def main(fastq_filename, end_row):
    # p = optparse.OptionParser(__doc__)
    # p.add_option("-n", dest="n", help="number of qual lines to test default:-1"
    #              " means test until end of file or until it it possible to "
    #              " determine a single file-type",
    #              type='int', default=-1)

    # opts, args = p.parse_args()

    # if len(args) > 1:
    #     print("Only a single input file is supported.", file=sys.stderr)
    #     sys.exit(1)

    gmin = 99
    gmax = 0
    valid = []

    err_exit = False

    input_file = fileinput.input(fastq_filename, openhook=fileinput.hook_compressed)

    for i, line in enumerate(input_file):      
        if (i+1)%4 ==0:
            lmin, lmax, qual_val_counts = get_qual_range(line.rstrip())
            if lmin < gmin or lmax > gmax:
                gmin, gmax = min(lmin, gmin), max(lmax, gmax)
                valid = get_encodings_in_range(gmin, gmax)
    
                # valid = heuristic_filter(valid, qual_val_counts)
            if len(valid) == 0:
                # print("no encodings for range: "
                #       "{}".format((gmin, gmax)), file=sys.stderr)
                err_exit = True
                break

            if end_row > 0 and i > end_row:
                # parsed up to specified portion; return current guess(es)
                break

    input_file.close()

    if err_exit:
        return ['chenlin', gmin, gmax]
        #sys.exit(1)
    else:
        return valid
        # print("{}\t{}\t{}".format(",".join(valid), gmin, gmax))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Detect quality format of fastq files')
    parser.add_argument('-n', '--end_row', type=int,default=-1,help='the number of row each fastq file read, default:-1')
    parser.add_argument('-o','--outputfilename',required=True, help='the output file name')
    parser.add_argument('-p','--password',  required=True, help='the password used in database connecting')
    parser.add_argument('-d','--dbname', default='pmdb', help='database name')
    parser.add_argument('-ho','--hostname', default='172.22.99.9', help='hostname of the db server')
    parser.add_argument('-s','--schema', default='sunset', help='the password used in database connecting')
    parser.add_argument('-u','--db_user', default='cl', help='database username')
    args = parser.parse_args()
    # from sqlalchemy import create_engine
    # engine = create_engine('postgresql://cl:%s@172.22.99.9/pmdb' %password)
    # from sqlalchemy.orm import sessionmaker
    # Session = sessionmaker(bind = engine)
    # session = Session()

    from pymodule.Sunset.db import SunsetDB
    db_main = SunsetDB.SunsetDB(drivername='postgresql', db_user=args.db_user,db_passwd=args.password, 
                            hostname=args.hostname, dbname=args.dbname, schema=args.schema)
    db_main.setup(create_tables=False)
    session = db_main.session

    from pymodule.db.SunsetDB import IndividualSequence, IndividualAlignment, IndividualSequenceFile
    f_out = open(args.outputfilename, 'w') 
    id_list = session.query(IndividualAlignment.ind_seq_id).all()
    for the_id in id_list:
        the_id = the_id[0]
        individual_id = session.query(IndividualSequence.individual_id)\
                           .filter(IndividualSequence.id == the_id).all()[0][0]
        raw_path = session.query(IndividualSequenceFile.path).filter(
            IndividualSequenceFile.individual_sequence_id == the_id).all()[1][0]

        ripe_path=os.path.join('/y/Sunset/db',raw_path)
        if not os.path.exists(ripe_path):
            ripe_path=os.path.join('/y/backup/db',raw_path)
            if not os.path.exists(ripe_path):
                sys.stderr.write("The path %s is not exist! "%ripe_path)
                exit(1)

        result = main(ripe_path, args.end_row)
        if 'Illumina-1.8' not in result and 'Sanger' not in result:
            f_out.write("%s\t" % the_id) 
            if 'chenlin' not in result:
                f_out.write("%s\n" % individual_id) 
            else:
                f_out.write("%s\t" % individual_id)
                f_out.write("no encodings for range: (%s,%s)\n" %(result[1],result[2]))

    f_out.close()  
    
