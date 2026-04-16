#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 12:25:29 2024

@author: khanhvu
"""

import torch
import tensorflow as tf
from tensorflow.keras.models import load_model
import numpy as np

esmModel, alphabet = torch.hub.load("facebookresearch/esm","esm1b_t33_650M_UR50S")

torch.cuda.empty_cache()
myModel = load_model("NbAgBindingModel.h5")

print(torch.cuda.is_available())
device = torch.device("cpu")

esmModel = esmModel.to(device)


def predict(seq):
    if(len(seq) >1024):
        return None
    
    # Start of encoding sequence in 1280 dimensions vector
    batch_tokens = alphabet.encode(seq)
    batch_tokens = torch.tensor([batch_tokens], dtype= torch.long)
    batch_tokens = batch_tokens.to(device)
    
    with torch.no_grad():
        result = esmModel(batch_tokens, repr_layers = [33], return_contacts = True)
        token_representations = result["representations"][33]
    
    seq_len = len(seq) + 1
    seq_embs = token_representations[0, 1: seq_len].mean(0)
    emb_list = seq_embs.tolist()
    emb_array = np.array(emb_list)
    
    emb_array = np.array(emb_array).reshape(1, -1)
    prediction = myModel.predict(emb_array)[0][0]
    
    
    print(prediction)
    if(prediction >=0.5):
        prediction = "Yes"
    else:
        prediction = "No"
    
    return prediction
    

print("Enter the nanobody sequence: (please ensure that the concatenated seq of antigen and nanobodies is less than 1024 amino acids)")
nanobody_seq = input()

print("Enter the antigen sequence: (please ensure that the concatenated seq of antigen and nanobodies is less than 1024 amino acids)")
antigen_seq = input()

concat_seq = nanobody_seq + antigen_seq

print(predict(concat_seq))

    