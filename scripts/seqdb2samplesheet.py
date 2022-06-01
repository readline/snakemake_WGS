#!/usr/bin/env python
import os,sys
import pandas as pd

def main():
    try:
        dbpath = sys.argv[1]
        batch = sys.argv[2]
    except:
        sys.exit(sys.argv[0] + ' [Seq db path] [Batch to include, sep by commas]')
    batch = [i.strip() for i in batch.split(',')]
    db = pd.read_csv(dbpath, sep='\t')
    print('Sample\tLib\tRun\tChip\tLane\tPlatform\tRead1\tRead2')
    for i in db.index:
        tmp = []
        if db.loc[i,'Batch'] not in batch:
            continue
        tmp.append(db.loc[i,'SM'])
        tmp.append('%s.%s'%(db.loc[i,'SM'], db.loc[i,'LB']))
        tmp.append('%s.%s'%(db.loc[i,'PU'], db.loc[i,'ID']))
        tmp.append(db.loc[i,'PU'].split('_')[0])
        tmp.append(db.loc[i,'PU'].split('_')[1])
        tmp.append(db.loc[i,'PL'])
        tmp.append(db.loc[i,'Read1'])
        tmp.append(db.loc[i,'Read2'])
        print('\t'.join(tmp))
if __name__ == '__main__':
    main()
