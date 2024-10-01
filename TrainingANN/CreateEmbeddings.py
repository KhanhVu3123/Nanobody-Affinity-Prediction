# -*- coding: utf-8 -*-
# This file is used to create embeddings value from a sequence.


import torch 
import csv
import gc
import numpy as np
model, alphabet = torch.hub.load("facebookresearch/esm", "esm2_t36_3B_UR50D")

file_list = ["/home/khanhvu/Desktop/Sequence/New_Wrong_concat_seq.txt", "/home/khanhvu/Desktop/Sequence/New_concat_seq.txt"]

for file_name in file_list:
    embeddings = []
    
    print(torch.cuda.is_available())
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    model = model.to(device)
    sequences = []
    seq = ""
    
    with open(file_name, 'r') as file:
      for line in file:
        if(line.startswith(">")):
          if(seq != ""):
             if(len(seq) <= 1024):
                 sequences.append(seq)
           #  else:
                 #sequences.pop()
          seq = line[1:].strip()
        else:
          seq = seq + line.strip()
      if(seq != ""):
        sequences.append(seq)
    
    
    
    for index in range(len(sequences)):
        
        child_sequences = sequences[index]
    
        batch_tokens = (alphabet.encode(child_sequences))
        batch_tokens = torch.tensor([batch_tokens], dtype= torch.long)
        
        
        batch_tokens = batch_tokens.to(device)
        
        with torch.no_grad():
          results = model(batch_tokens, repr_layers = [36])
          token_representations =results["representations"][36].squeeze(0)
        
        # for i, token_seq in enumerate(batch_tokens):
            
        """seq_len = len(child_sequences) +1
        seq_embs = token_representations[0, 1:seq_len].mean(0)
        emb_list = seq_embs.tolist()
        writer.writerow([index,emb_list])"""
        
        print(token_representations.cpu().numpy().shape)
        
        embeddings.append(token_representations.cpu().numpy())
        
        print(index/ len(sequences) * 100, end = " ")
        print("%")
        
        del batch_tokens
        del results
        del token_representations
        """del seq_embs
        del emb_list"""
        gc.collect()
        
    if(file_name == "/home/khanhvu/Desktop/Sequence/New_Wrong_concat_seq.txt"):
        np.save("/home/khanhvu/Desktop/Embeddings/New_Wrong_concat_embeddings.npy",embeddings)
    else:
        np.save("/home/khanhvu/Desktop/Embeddings/New_concat_embeddings.npy",embeddings)

print("end of program")
