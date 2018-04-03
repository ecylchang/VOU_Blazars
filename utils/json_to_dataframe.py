#!/usr/bin/env python

import pandas
import json
import sys
import os


def read_json(filename):

    with open(filename,'r') as fp:
        results = json.load(fp)

    return results


def dict_2_df(results):

    dict_rows = {(i,j): results[i][j]
                        for i in results.keys()
                        for j in results[i].keys()}

    df = pandas.DataFrame.from_dict(dict_rows, orient='index')
    df.index.names = ['RUN_LABEL','SOURCE']
    print(df)
    
    return df


def json_2_df(json_file):

    d = read_json(json_file)
    
    df = dict_2_df(d)

    return df


if __name__ == '__main__':
    import sys

    try:
        json_file = sys.argv[1]
        df = json_2_df(json_file)
    except:
        print('\n'
              'Usage: {} <results.json> [results.csv]'
              '\n'
              .format(os.path.basename(sys.argv[0]))
              )
        sys.exit(1)
    try:
        csv_file = sys.argv[2]
    except:
        csv_file = 'results.csv'

    df.to_csv(csv_file)

