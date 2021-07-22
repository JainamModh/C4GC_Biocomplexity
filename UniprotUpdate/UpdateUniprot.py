# Jainam Modh
# Andrew Warren, UVA Biocomplexity Insititute
# jhm3ab@virginia.edu

import requests
import pandas as pd
import json
import sys

# Read in feature information from PATRIC
patric = pd.read_table(sys.argv[1])
patric = patric.rename(columns = {'genome.genome_id':'genome_id', 
    'feature.patric_id':'patric_id', 
    'feature.aa_sequence_md5':'md5',
    'feature.pgfam_id':'pgfam_id',
    'feature.plfam_id':'plfam_id'
})
patric.dropna(subset = ['md5'], inplace = True)
patric['taxon_id'] = (patric['genome_id'].astype(str).apply(lambda x: x.split('.')[0])).astype(int)

#Extract md5 and taxon ID information
patricmd5 = patric['md5'].values
patrictaxon = patric['taxon_id'].drop_duplicates().values


# Create md5 and taxon ID queries
def md5toString(md5array):
    s = '('
    for md5 in md5array:
        s = s + '"' + md5 + '",'
    s = s[:-1] + ')'
    return s
def taxontoString(taxon_id):
    s = '('
    for org in taxon_id:
        s = s + 'taxon:' + str(org) + ','
    s = s[:-1] + ')'
    return s


# query to uniprot SPARQL endpoint for uniprot_id, aa_sequence, goterms, and golabels. 
search = f"""
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
        SELECT DISTINCT
            (CONCAT(SUBSTR(STR(?protein), 33)) AS ?uniprot)
            ?md5 
            ?aa_sequence
            (CONCAT(SUBSTR(STR(?taxon), 34)) AS ?taxon_id)
            (CONCAT(SUBSTR(STR(?goTerm), 32)) AS ?GO)
            ?goLabel
        WHERE {{
    VALUES (?go) {{("GO_0008150") ("GO_0005575") ("GO_0003674")}}
    BIND (IRI(CONCAT("http://purl.obolibrary.org/obo/", ?go)) AS ?aspect)
            ?protein a up:Protein ;
                up:organism ?taxon ;
                up:classifiedWith ?goTerm ;
                up:sequence ?sequence .
            ?sequence up:md5Checksum ?md5 ;
                      rdf:value ?aa_sequence .
            ?goTerm rdfs:subClassOf ?aspect ;
                    rdfs:label ?goLabel .
        FILTER (?md5 in {md5toString(patricmd5)})
        FILTER (?taxon in {taxontoString(patrictaxon)})
        }}
"""
r1 = requests.post("https://sparql.uniprot.org", headers={'accept': 'application/sparql-results+json'}, data={'query': search})
results = r1.json()


# Store information from query into dataframe
uniprot = pd.DataFrame(columns=('md5', 'uniprot_id', 'taxon_id', 'aa_sequence', 'goTerm', 'goLabel'))
for row in results['results']['bindings']:
    uniprot = uniprot.append({
        'md5': row['md5']['value'], 
        'uniprot_id': row['uniprot']['value'],
        'taxon_id': int(row['taxon_id']['value']),
        'aa_sequence': row['aa_sequence']['value'],
        'goTerm': row['GO']['value'],
        'goLabel': row['goLabel']['value']},
        ignore_index=True)


# Merge dataframes and create linkout crosslinks
patricGO = pd.merge(uniprot, patric)
crosslinks = patricGO[['uniprot_id', 'patric_id']].drop_duplicates()
crosslinks.to_csv(sys.argv[2], index=False, sep = '\t')


# Create json fasta file for PATRIC annotation
fasta = []
for md5 in patricGO['md5'].unique():
    frame = patricGO.loc[patricGO['md5'] == md5]
    data = {
        'md5': md5,
        'sequence': frame['aa_sequence'].iloc[0],
        'sequence_type': 'protein'
    }
    data['patric'] = []
    for id in frame['patric_id'].unique():
        frame2 = frame.loc[frame['patric_id'] == id]
        p_ids = {
            'patric_id': id,
            'pgfam_id': frame2['pgfam_id'].iloc[0],
            'plfam_id': frame2['plfam_id'].iloc[0]}
        p_ids['GO'] = []
        for i in range(len(frame2)):
            p_ids['GO'].append({
                'GO_id': frame2['goTerm'].iloc[i],
                'GO_label': frame2['goLabel'].iloc[i]
            })
        data['patric'].append(dict(p_ids))
    fasta.append(dict(data))
with open(sys.argv[3], 'w') as f:
    json.dump(fasta, f)
