#!/usr/bin/env python
# coding: utf-8

# In[1]:


from urllib.request import urlopen
from bs4 import BeautifulSoup
from Bio.KEGG.KGML.KGML_parser import read
import pandas as pd
import os
from collections import defaultdict
from tqdm import tqdm
import requests
import networkx as nx
import numpy as np


# In[2]:


# url = 'http://rest.kegg.jp/link/hsa/pathway'
# html = urlopen(url).read().decode("utf-8")


# In[3]:


# soup = BeautifulSoup(html, 'html.parser')


# In[4]:


# finding KEGG ID's of every gene (?) here
kegg_dict = {
    'Zap70' : 7535,
    'SHP2' : 5781,
    'Stat1' : 6772,
    'Akt' : 207,
    'Stat5a' : 6776,
    'Stat5b' : 6777,
    'p38' : 1432,
    'NFkB1' : 4790,
    'NFKB2' : 4791,
    'NFkBIA' : 4792,
    'S6' : 6194,
    'Lat' : 27040,
    'Erk' : 5594,
    'Plcg2' : 5336,
    'Btk' : 695,
    'Slp76' : 3937,
    'Stat3' : 6774
}


# In[ ]:





# In[6]:


# pathways = []
# for line in html.split('\n'):
#     for kegg_id in kegg_dict.values():
#         if ('hsa:' + str(kegg_id)) in line:
#             pathway = line.split('\t')[0]
#             pathways.append(pathway)
# pathways = list(set(pathways))


# In[7]:


# len(pathways)


# In[8]:


# pathways[:5]


# In[9]:


# # files persisted
# for pathway in tqdm(pathways):
#     url = 'http://rest.kegg.jp/get/' + pathway.split(':')[1] + '/kgml'
#     response = requests.get(url)
#     with open('pathways/' + pathway.split(':')[1] + '.xml', 'wb') as file:
#         file.write(response.content)


# In[10]:


test_pathway = read(open('pathways/hsa04270.xml', 'r'))


# In[11]:


test_pathway.relations[0].entry1.name


# In[12]:


gene_hsa = ['hsa:' + str(val) for val in kegg_dict.values()]
gene_hsa_dict = {('hsa:' + str(v)) : k for k,v in kegg_dict.items()}


# In[13]:


link_list = []
for kgml in tqdm(os.listdir('pathways/')):
    # check for xml
    if kgml[-4:] != '.xml':
        continue
    path_name = kgml[:-4]
    
    # read file
    path = read(open('pathways/' + kgml, 'r'))
    
    # iterate through edges
    for rel in path.relations:
        link_list.append({
            'pathway' : path_name,
            'entry1' : rel.entry1.name,
            'entry2' : rel.entry2.name,
            'type' : rel.type,
            'subtype' : rel.subtypes
        })


# In[14]:


links = pd.DataFrame(link_list)
links = links[(links['entry1'] != 'undefined') & (links['entry2'] != 'undefined')]
links


# In[ ]:





# In[15]:


# create a graph for each pathway from links
all_paths = {}
for kegg_path in tqdm(np.unique(links['pathway'])):
    G = nx.DiGraph()
    for i in list(links[links['pathway'] == kegg_path].index):
        row = links.loc[i]
        G.add_edge(row['entry1'], row['entry2'])
    all_paths[kegg_path] = (nx.algorithms.shortest_path(G))


# In[ ]:





# In[17]:


counter = 0
for pathway in tqdm(all_paths):
    pathway_paths = all_paths[pathway]
    keys = list(pathway_paths.keys())
    for key in keys:
        # check if source node is relevant protein
        relevant = False
        for gene in gene_hsa:
            if gene in key:
                relevant = True
                # now check if target relevant
                for key2 in list(pathway_paths[key].keys()):
                    relevant2 = False
                    for gene in gene_hsa:
                        if gene in key2:
                            relevant2 = True
                    if not relevant2:
                        # if target not relevant, remove it from dict
                        pathway_paths[key].pop(key2)
        if not relevant:
            # if source not relevant, remove it from dict
            pathway_paths.pop(key)


# In[18]:


def find_node(node_name):
    '''
    Given a string like 'hsa:2023 hsa:2026 hsa:2027 hsa:387712', return which proteins match
    '''
    proteins = []
    for protein in gene_hsa_dict.keys():
        if protein in node_name:
            proteins.append(gene_hsa_dict[protein])
    return proteins


# In[19]:


relevant_links = []
for pw in all_paths.keys():
    pathway = all_paths[pw]
    for node1 in pathway.keys():
        for node2 in pathway[node1].keys():
            node1_proteins = find_node(node1)
            node2_proteins = find_node(node2)
            for protein in node1_proteins:
                for protein2 in node2_proteins:
                    relevant_links.append({
                        'node1' : protein,
                        'node2' : protein2,
                        'pathway' : pw,
                        'path': pathway[node1][node2]
                    })
relevant_links = pd.DataFrame(relevant_links)
relevant_links = relevant_links[relevant_links['node1'] != relevant_links['node2']].reset_index(drop = True)
relevant_links.to_csv('kegg_links.csv', index = False)


# In[20]:


relevant_links


# In[ ]:





# In[ ]:




