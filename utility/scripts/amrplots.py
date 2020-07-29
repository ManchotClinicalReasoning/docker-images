import os
import json
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description='Create categorical AMR plots')
parser.add_argument('-d', type=str,help='card json')
parser.add_argument('-c', type=str,help='categories')
parser.add_argument('-s', type=str,nargs='+',help='samples')
parser.add_argument('-f', type=str,help='figure')
parser.add_argument('-t', type=str,help='table')
args = parser.parse_args()


#Load CARD


CARD = json.load(open(args.d,'r'))
#Restructure, key by AROTerm

CARDByName = {CARD[k]['ARO_name'] : CARD[k] for k in CARD if 'ARO_name' in CARD[k]}
CARDByID = {CARD[k]['ARO_accession'] : CARD[k] for k in CARD if 'ARO_id' in CARD[k]}

drugs = json.load(open(args.c,'r'))

data = []

for f in args.s:
    time = int(f.split('/')[-3].split('_')[1])

    AROTermHits = {}

    with open(f,'r') as infile:
        for l in infile.read().splitlines():
            AROTerm = l.split('\t')[2]
            if not AROTerm in AROTermHits:
                AROTermHits[AROTerm] = 0
            AROTermHits[AROTerm] += 1
            
    print(AROTermHits)
    print('Total Hits: {}'.format(sum(AROTermHits.values())))

    drugHits = { d : 0 for d in drugs}

    for term in AROTermHits:
        print('Processing term: {}'.format(term))
        if not term in CARDByName:
            print('{} not contained in CARD ontology, skipping ... '.format(term))
            continue
        CARDEntry = CARDByName[term]
        for annotation in CARDEntry['ARO_category'].values():
            if annotation['category_aro_class_name'] in ['Drug Class','Antibiotic']:
                hitProcessed = {d : False for d in drugs}
                #Fetch corresponding drugs
                drugClassID = annotation['category_aro_accession']
                for drug in drugs:
                    if drugClassID in drugs[drug] and hitProcessed[drug] == False:
                        hitProcessed[drug] = True
                        drugHits[drug] += AROTermHits[term]
    for drug in drugHits:
        data.append((time,drug,drugHits[drug]))

df = pd.DataFrame(data, columns =['time', 'drug', 'count'])
df.to_csv(args.t)
plt.figure()
sns.lineplot(x='time',y='count',color='drug',data=df)
plt.savefig(args.f)
