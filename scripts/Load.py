#!/usr/bin/env python
import sys
import pandas as pd

def samplesheet(sspath):
    df = pd.read_csv(sspath, sep='\t', usecols=range(12))
    idx = []
    for i in df.index:
        idx.append('%s.%s.%s'%(df.loc[i,'SM'], df.loc[i,'LB'], df.loc[i,'PU']))
    df.index=idx
    sample = {}
    lib = {}
    run = {}
    for i in df.index:
        ss = str(df.loc[i,'SM'])
        li = str(df.loc[i,'LB'])
        if ss not in sample:
            sample[ss] = []
        if li not in lib:
            lib[li] = []
        if li not in sample[ss]:
            sample[ss].append(li)
        if i not in lib[li]:
            lib[li].append(i)
        run[i] = {n:df.loc[i,n] for n in df.columns}
    return sample,lib,run

if __name__ == '__main__':
    result = samplesheet(sys.argv[1])
    import json
    for i in result:
        print(json.dumps(i, indent=4, sort_keys=True))
        print()
