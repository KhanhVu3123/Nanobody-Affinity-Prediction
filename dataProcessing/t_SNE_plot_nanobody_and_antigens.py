#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 12:27:02 2023

@author: khanhvu
"""

import ProcessPdbFile
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt


# We are now retreving the nanobody embeddings

df = pd.read_csv("/home/khanhvu/Desktop/antigen_emb.csv")

df['Embeddings'] = df['Embeddings'].apply(lambda x: np.fromstring(x.strip('[]'), sep = ','))
embeddings = np.stack(df["Embeddings"].values)

tsne = TSNE(n_components= 2, perplexity= 30, n_iter= 300)
tsne_results = tsne.fit_transform(embeddings)

plt.figure(figsize=(10,6))
plt.scatter(tsne_results[:,0], tsne_results[:, 1])

plt.xlabel('t-SNE component 1')
plt.ylabel('t-SNE component 2')
plt.title('t-SNE Plot of Embeddings')

plt.savefig('tSNEPlot_of_antigen.pdf', format = 'pdf', bbox_inches = 'tight')

plt.show()



"""
antigen_seq_dict = ProcessPdbFile.retrieve_all_antigenSeq_fromzippedFile("/home/khanhvu/Desktop/UpdatedPDBFile")
with open("/home/khanhvu/Desktop/Sequence/Antigen_seq.txt", "w") as file:
    for name in antigen_seq_dict.values():
        file.write(">" + name + "\n")
        
"""



