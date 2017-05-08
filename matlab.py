#! /usr/bin/env python

import os, sys
import argparse


def get_parser():
    parser = argparse.ArgumentParser(description='Matlab_Rewrite')
    subparsers = parser.add_subparsers(dest='subcommand', help='Select one of the following sub-commands')

    parser_a = subparsers.add_parser('dictionary', help='estimate purity using adaptive windows and dictionaries')

    parser_b = subparsers.add_parser('gene', help='analysis of all single genes in the sample')

    parser_c = subparsers.add_parser('window', help='based on adaptive window approach')
    parser_c.add_argument('-size', default=10, help='the size of window')

    print parser.parse_args()

get_parser()
