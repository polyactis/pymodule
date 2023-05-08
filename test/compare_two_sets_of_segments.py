#!/usr/bin/env python3
"""
2023/5/8 Yu S. Huang polyactis@gmail.com
A program that compares the query efficiency between two different databases of \
    identical content.
    1) One database is constructed using a list.
    2) The other database is a red-black tree/dictionary.

Prerequisites:
    1. pip3 install palos tqdm

Example:
    ./compare_two_sets_of_segments.py -q HCC1187_10X_0_3_chr1.segments.M10.T5.tsv -d HCC1187_10X_0_3_chr1.segments.M20.T10.tsv

"""

import os, sys
import time
from palos import RBDict
from palos.polymorphism.CNV import CNVCompare, CNVSegmentBinarySearchTreeKey, \
    is_reciprocal_overlap, leftWithinRightAlsoEqualCmp
from tqdm import tqdm



class compare_two_sets_of_segments(object):
    def __init__(self, query_segment_filepath=None, db_segment_filepath=None, \
                 min_reciprocal_overlap=0.01, debug=False) -> None:
        self.query_segment_filepath = query_segment_filepath
        self.db_segment_filepath = db_segment_filepath
        self.min_reciprocal_overlap = min_reciprocal_overlap
        self.debug = debug

    def construct_segment_list_from_file(self, filepath) -> list:
        """
        """
        print(f"Constructing a segment list from {filepath} ... ", \
              file=sys.stderr)
        start_time = time.time()
        segment_list = []
        input_file = open(filepath)
        for line in input_file:
            if line[0] == '#':
                continue
            else:
                row = line.strip().split()
                chr = row[0]
                start = int(row[1])
                stop = int(row[2])
                segment_list.append((chr, start, stop))

        print(f"\t {len(segment_list)} segments in "
            f"{time.time()-start_time:.3f} seconds.\n", file=sys.stderr)
        return segment_list

    def construct_segment_rbdict_from_segment_list(self, segment_list) -> RBDict:
        """
        """
        start_time = time.time()
        print(f'Constructing segment rb dict from a list of {len(segment_list)} '
              f'segments ...', file=sys.stderr)

        #pbar = tqdm(total=len(segment_list))

        # leftWithinRightAlsoEqualCmp: During node insertion, if a new segment (left)
        #   is found tpbarpbaro be wholly contained by an old segment (right, alreaydy
        #   inserted into the RB tree), an error will be raised.
        rb_dict = RBDict(cmpfn=leftWithinRightAlsoEqualCmp)
        for segment in segment_list:
            chromosome, start, stop = segment[:3]
            segment_key = CNVSegmentBinarySearchTreeKey(chromosome=chromosome,
                span_ls=[start, stop], \
                min_reciprocal_overlap=self.min_reciprocal_overlap)
            rb_dict[segment_key] = segment
            #pbar.update(1)
        print(f"\t Done in {time.time()-start_time:.3f} seconds.\n", file=sys.stderr)
        return rb_dict

    def construct_segment_rbdict_from_file(self, filepath) -> RBDict:
        """
        """
        print(f"Constructing segment rb dict from {filepath} ... ", file=sys.stderr)
        start_time = time.time()
        input_file = open(filepath)
        # leftWithinRightAlsoEqualCmp: During node insertion, if a new segment (left)
        #   is found to be wholly contained by an old segment (right, alreaydy
        #   inserted into the RB tree), an error will be raised.
        rb_dict = RBDict(cmpfn=leftWithinRightAlsoEqualCmp)
        for line in input_file:
            if line[0] == '#':
                continue
            else:
                row = line.strip().split()
                chr = row[0]
                start = int(row[1])
                stop = int(row[2])
                segment_key = CNVSegmentBinarySearchTreeKey(chromosome=chr,
                    span_ls=[start, stop], \
                    min_reciprocal_overlap=self.min_reciprocal_overlap)
                rb_dict[segment_key] = 1
        print(f"\t {len(rb_dict)} segments in the red-black tree"\
                f" in {time.time()-start_time:.3f} seconds.\n", file=sys.stderr)
        return rb_dict
    
    def query_against_rbtree(self, query_segment_ls, segment_rbtree_db) -> int:
        """
        """
        no_of_total_query_segments = len(query_segment_ls)
        print(f"Querying {no_of_total_query_segments} against a rbtree db of "
              f"{len(segment_rbtree_db)} entries ...", file=sys.stderr)
        start_time = time.time()
        compare_function = CNVCompare(min_reciprocal_overlap=self.min_reciprocal_overlap)
        no_of_query_with_hits = 0
        for query_segment in query_segment_ls:
            query_key = CNVSegmentBinarySearchTreeKey(chromosome=query_segment[0], \
                            span_ls=[query_segment[1], query_segment[2]])
            hit_node_ls = []
            segment_rbtree_db.findNodes(query_key, node_ls=hit_node_ls, \
                                        compareIns=compare_function)
            if len(hit_node_ls)>0:
                no_of_query_with_hits += 1
        print(f"\t {no_of_query_with_hits}/{no_of_total_query_segments} query "
            f"segments have hits in the db. "\
            f"Cost {time.time()-start_time:.3f} seconds.\n", file=sys.stderr)
        return no_of_query_with_hits
    
    def query_against_list(self, query_segment_ls, db_segment_ls):
        no_of_total_query_segments = len(query_segment_ls)
        print(f"Querying {no_of_total_query_segments} against a list of "
              f"{len(db_segment_ls)} entries ...", file=sys.stderr)
        start_time = time.time()
        pbar = tqdm(total=len(query_segment_ls)*len(db_segment_ls))
        no_of_query_with_hits = 0
        for query_segment in query_segment_ls:
            query_chr = query_segment[0]
            for db_segment in db_segment_ls:
                db_chr = db_segment[0]
                if query_chr == db_chr and is_reciprocal_overlap(span1_ls=query_segment[1:3], \
                        span2_ls=db_segment[1:3], 
                        min_reciprocal_overlap=self.min_reciprocal_overlap):
                    no_of_query_with_hits += 1
                    #break
            pbar.update(len(db_segment_ls))
        print(f"\t {no_of_query_with_hits}/{no_of_total_query_segments} query "
            f"segments have hits in the db. "\
            f"Cost {time.time()-start_time:.3f} seconds.\n", file=sys.stderr)
        return no_of_query_with_hits

    def run(self):
        db_segment_ls = self.construct_segment_list_from_file(self.db_segment_filepath)
        db_segment_rbdict = self.construct_segment_rbdict_from_segment_list(db_segment_ls)
        query_segment_ls = self.construct_segment_list_from_file(self.query_segment_filepath)
        # query_segment_rbdict = self.construct_segment_rbdict_from_file(self.query_segment_filepath)
        # query_segment_rbdict_2 = self.construct_segment_rbdict_from_segment_list(query_segment_ls)
        self.query_against_rbtree(query_segment_ls, db_segment_rbdict)
        self.query_against_list(query_segment_ls, db_segment_ls)


if __name__ == '__main__':
    from argparse import ArgumentParser
    ap = ArgumentParser()
    ap.add_argument("-q", "--query_segment_filepath", type=str, required=True,
        help="the path to the query file containing genomic segments.")
    ap.add_argument("-d", "--db_segment_filepath", type=str, required=True,
        help="the path to the database file containing genomic segments.")
    
    ap.add_argument("-m", "--min_reciprocal_overlap", type=float, default=0.0001,
        help='Default: %(default)s. '
        'If a job virtual memory (usually 1.2X of JVM resident memory) exceeds request, '
        "it will be killed on some clusters. "
        "This will make sure your job requests enough memory.")
    
    
    ap.add_argument("--debug", action='store_true',
        help='Toggle debug mode.')
    args = ap.parse_args()
    instance = compare_two_sets_of_segments(
        query_segment_filepath=args.query_segment_filepath, \
        db_segment_filepath=args.db_segment_filepath, \
        min_reciprocal_overlap=args.min_reciprocal_overlap,
        debug=args.debug)
    instance.run()
