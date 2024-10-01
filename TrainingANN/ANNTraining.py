#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 11:34:42 2023

@author: khanhvu
"""

import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import train_test_split

from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import keras_tuner as kt
from sklearn.model_selection import KFold
import optuna

# Start loading the dataset
truepair_file_path= '/home/khanhvu/Desktop/Embeddings/True_Embed15Jan.csv'
wrongpair_file_path = '/home/khanhvu/Desktop/Embeddings/False_Embed15Jan.csv'

true_embedding = []
truepair_df = pd.read_csv(truepair_file_path)
truepair_array = list(truepair_df["Embeddings"])

false_embedding = []
falsepair_df = pd.read_csv(wrongpair_file_path)
falsepair_array = list(falsepair_df["Embeddings"])


for line in truepair_array:
    values = line.strip()[1:-1].split(", ")
    true_embedding.append([float(value) for value in values])

for line in falsepair_array:
    values = line.strip()[1:-1].split(", ")
    false_embedding.append([float(value) for value in values])

        
true_embedding = np.array(true_embedding)
false_embedding = np.array(false_embedding)

# Finish loading the dataset


true_labels = np.ones(len(true_embedding))
false_labels = np.zeros(len(false_embedding))

# Combine the data and the labels

X = np.concatenate((true_embedding, false_embedding))
y = np.concatenate([true_labels,false_labels])

print(X.shape)


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)

X_tensor = torch.tensor(X, dtype = torch.float32)
y_tensor = torch.tensor(y, dtype = torch.float32)

# Split the data
X_train, X_test, y_train, y_test = train_test_split(X_tensor, y_tensor, test_size= 0.2, random_state= 42)

train_dataset = TensorDataset(X_train, y_train)
val_dataset = TensorDataset(X_test, y_test)

train_loader = DataLoader(train_dataset, batch_size= 64, shuffle = True)
val_loader = DataLoader(val_dataset, batch_size= 64, shuffle= False)

# Create a binary classification artificial neural network with 3 layers 
class Net(nn.Module):
    def __init__(self, dropout_rate = 0.5):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(1280, 256)
        self.dropout1 = nn.Dropout(dropout_rate)
        self.fc2 = nn.Linear(256, 128)
        self.dropout2 = nn.Dropout(dropout_rate)
        self.fc3 = nn.Linear(128,1)
        
    def forward(self,x):
        x = torch.relu(self.fc1(x))
        x = self.dropout1(x)
        x = torch.relu(self.fc2(x))
        x = self.dropout2(x)
        x = torch.sigmoid(self.fc3(x))
        return x
    
# Test different learning rate and dropout rate of each layer.
criterion = nn.BCELoss()

def objective(trial):
    lr = trial.suggest_float("lr", 1e-5, 1e-1, log = True)
    dropout_rate = trial.suggest_float("dropout_rate", 0.1, 0.5)
    
    model = Net(dropout_rate = dropout_rate).to(device)
    optimizer = optim.Adam(model.parameters(), lr = lr)


    
    
    num_epochs = 10
    for epoch in range(num_epochs):
        model.train()
        for inputs, labels in train_loader:
            inputs, labels = inputs.to(device), labels.to(device)
            optimizer.zero_grad()
            outputs = model(inputs)
            loss = criterion(outputs.squeeze(), labels)
            loss.backward()
            optimizer.step()
        
        model.eval()
        with torch.no_grad():
            val_loss = 0
            for inputs, labels in val_loader:
                inputs, labels = inputs.to(device), labels.to(device)
                outputs = model(inputs)
                loss = criterion(outputs.squeeze(), labels)
                val_loss += loss.item()
                
        print("Epoch " + str(epoch +1) )
        
    y_true = []
    y_scores = []
    threshold = 0.5
    with torch.no_grad():
        correct =0 
        total = 0
        for inputs, labels in val_loader:
            inputs, labels = inputs.to(device), labels.to(device)
            outputs = model(inputs)
            
            outputs = outputs.squeeze().cpu()
            labels = labels.cpu()
            
            y_scores.extend(outputs.numpy().flatten())
            y_true.extend(labels.numpy().flatten())
            predicted = (outputs > threshold).float()
            
            total += labels.size(0)
            
            correct += (predicted == labels).sum().item()
            accuracy = correct/total * 100
            
    return accuracy
    
study = optuna.create_study(direction = "maximize")
study.optimize(objective, n_trials = 100) # Try 100 times 

best_params = study.best_params
model = Net(dropout_rate = best_params['dropout_rate']).to(device)

optimizer = optim.Adam(model.parameters(), lr = best_params['lr'])
num_epochs = 10
for epoch in range(num_epochs):
    model.train()
    for inputs, labels in train_loader:
        inputs, labels = inputs.to(device), labels.to(device)
        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs.squeeze(), labels)
        loss.backward()
        optimizer.step()
    
    model.eval()
    with torch.no_grad():
        val_loss = 0
        for inputs, labels in val_loader:
            inputs, labels = inputs.to(device), labels.to(device)
            outputs = model(inputs)
            loss = criterion(outputs.squeeze(), labels)
            val_loss += loss.item()
            
    print("Epoch " + str(epoch +1) )
    
y_true = []
y_scores = []
threshold = 0.5
with torch.no_grad():
    correct =0 
    total = 0
    for inputs, labels in val_loader:
        inputs, labels = inputs.to(device), labels.to(device)
        outputs = model(inputs)
        
        outputs = outputs.squeeze().cpu()
        labels = labels.cpu()
        
        y_scores.extend(outputs.numpy().flatten())
        y_true.extend(labels.numpy().flatten())
        predicted = (outputs > threshold).float()
        
        total += labels.size(0)
        
        correct += (predicted == labels).sum().item()
        accuracy = correct/total * 100

fpr, tpr, _ = roc_curve(y_true, y_scores)
roc_auc = auc(fpr, tpr)

plt.figure()
plt.plot(fpr, tpr, lw =2, color = "darkorange", label = "ROC curve")
plt.plot([0,1], [0,1], color = 'navy', lw =2, linestyle = '--') 
plt.xlim([0.0,1.0])
plt.ylim([0.0, 1.0])
plt.xlabel("False positive rate ")
plt.ylabel("True positive rate")
plt.title("ROC curve")
plt.legend(loc = "lower right")
plt.show()



torch.save(model.state_dict(),"NbAgBindingModel.pth")

print("Accuracy: " + str(round(accuracy, 2)))

