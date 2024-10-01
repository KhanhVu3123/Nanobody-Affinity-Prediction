#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 12:17:57 2023

@author: khanhvu
"""
import ProcessPdbFile
import FindingCDRRegion
import pandas as pd
import numpy as np
from collections import Counter
import os

def add_percentage_for_DataFrame(dataframe):
    """
    

    Parameters
    ----------
    dataframe : The numerical dataframe
    Returns
    -------
    normalized_df : The processed DataFrame which element is a percentage 
    value with 3 decimal places.

    """
    total_sum = dataframe.to_numpy().sum()
    normalized_df = dataframe/total_sum *100
    normalized_df = normalized_df.round(3)
    normalized_df = normalized_df.applymap(lambda x: str(x) + "%")
    return normalized_df


def find_bonding_pair_list(pdb_file, chain, distance):
    """

    Parameters
    ----------
    pdb_file : File you want to investigate
    chain : The chain letter of the nanobody
    distance : Threshold that will be used to determine if this is a bond or not.

    Returns
    -------
    List of bonding pair of the nanobody sequence in the file

    """
    bonding_list = []
    big_dict = ProcessPdbFile.get_sequence(pdb_file) 
    
    if(chain not in big_dict.keys()):
        return "Check your data"
    aminolst = big_dict[chain]
    sequence = ""
    for d in aminolst:
        sequence = sequence + list(d.keys())[0]

    lst = FindingCDRRegion.find_CDR(sequence)
    aminolst = aminolst[lst[0] : lst[1] +1] + aminolst[lst[2] : lst[3] +1] + aminolst[lst[4] : lst[5] +1]


    # Aminolst is a list of dictionaries of coordiantes of every amino acids
    if(len(big_dict) <=1):
        return "Your dictionary is either wrong or you just have 1 chain"
    for chain_name in big_dict.keys():
        if chain_name == chain:
            continue
        aminolst2 = big_dict[chain_name]
        # AminoDict is a list, not a dictionary

        # Loop through every dictionary of the list.
        for amino_acid_coordinate in aminolst:
            # This is a dictionary
            amino_acid1 = next(iter(amino_acid_coordinate))
            lst1 = amino_acid_coordinate[next(iter(amino_acid_coordinate))]

            lst1 = ProcessPdbFile.to_float(lst1)
            # Change all of its elements to integers.

            for amino_acid_coordinate2 in aminolst2:
                amino_acid2 = next(iter(amino_acid_coordinate2))
                lst2 = amino_acid_coordinate2[next(iter(amino_acid_coordinate2))]
                lst2 = ProcessPdbFile.to_float(lst2)
                distance_squared = 0
                for i in range(3):
                    distance_squared = distance_squared+ ((lst1[i] - lst2[i]) ** 2)
                if(distance_squared < distance ** 2):
                    bonding_list.append(amino_acid1+ amino_acid2)
    """
    Function output examples: [AC, TK, VW, AC]
    Note that one bond can be repeated many times.
    """
    return bonding_list

def find_self_pairing_list(pdb_file, threshold_distance):
    
    """
    
    Parameters
    ----------
    pdb_file : The path to the PDB file
    threshold_distance: The value which any pair of amino acid
    that are less than this will be included in the output.

    Returns
    -------
    The list of bonding appears INTRAMOLECULAR

    """
    
    big_lst = []
    
    big_dict = ProcessPdbFile.get_sequence(pdb_file)
    for chain_name in big_dict.keys():
        all_coordinates_list = big_dict[chain_name]
        for amino_acid_dict1 in all_coordinates_list:
            amino_acid1 = next(iter(amino_acid_dict1))
            coordinate1 = amino_acid_dict1[amino_acid1]
            for amino_acid_dict2 in all_coordinates_list:
                amino_acid2 = next(iter(amino_acid_dict2))
                coordinate2 = amino_acid_dict2[amino_acid2]
                
                distance = ProcessPdbFile.coordinate_distance(coordinate1, coordinate2)
                
                if(distance <= threshold_distance and distance > 0):
                    string= amino_acid1 + amino_acid2
                    big_lst.append(string)
    return big_lst
    
def constructing_self_bonding_stats_table(base_dir, nanobody_file):
    """
    

    Parameters
    ----------
    base_dir : Full path to all the PDB files
    nanobody_file : The file that has PDB file with its nanobody chain

    Returns
    -------
    Pandas dataframe of the percentage of appearances of each of the pairing
    of the amino acids.

    """
    name_list = []
    with open(nanobody_file, 'r') as file:
        for line in file:
            name = line[:4]
            name_list.append(name)
    
    counting_variable = 0
    bonding_list = []
    for name in name_list:
        print(counting_variable/ len(name_list) * 100 , end = " ")
        print("%")
        
        this_directory = base_dir + "/" + name + ".pdb"
        bonding_list += find_self_pairing_list(this_directory, 8)
        counting_variable +=1
        
    print(bonding_list)
    bonding_list = clean_self_bonding_list(bonding_list)
    count_dict = Counter(bonding_list)
    
    rows = set(s[0] for s in bonding_list)
    cols = set(s[1] for s in bonding_list)
    
    df = pd.DataFrame(index = sorted(rows), columns = sorted(cols))
    
    for s, count in count_dict.items():
        row,col = s[0], s[1]
        df.at[row,col] = count
    df = df.fillna(0)
    df = add_percentage_for_DataFrame(df)
    df.to_csv("/home/khanhvu/Desktop/self_bonding.csv")
    
    return df

def constructing_bonding_stats_table(base_dir, nanobody_file):
    """

    Parameters
    ----------
    base_dir : Full path to all the PDB files
    nanobody_file : The file that has PDB file with its nanobody chain

    Returns
    -------
    Pandas dataframe of the percentage of appearances of each of the pairing
    of the amino acids.

    """
    lst= []
    with open(nanobody_file, "r") as file:
        for line in file:
            name = line[:4]+ line[22]
            lst.append(name)
    
    bonding_list = []
    i= 0
    for name in lst:
        print(i/len(lst) * 100, end =" ")
        print("%")
        dir_now = base_dir+ "/" + name[:4] + ".pdb"
        bonding_list += find_bonding_pair_list(dir_now, name[-1], 8)
        i = i+1
    count_dict = Counter(bonding_list)
    
    rows = set(s[0] for s in bonding_list)
    cols = set(s[1] for s in bonding_list)
    
    df = pd.DataFrame(index = sorted(rows), columns= sorted(cols))
    
    for s, count in count_dict.items():
        row,col = s[0],s[1]
        df.at[row, col] = count
        
    df = df.fillna(0)
    
    df = add_percentage_for_DataFrame(df)
    df.to_csv("/home/khanhvu/Desktop/out.csv")
    return df
    

def clean_self_bonding_list(bonding_list):
    new_list = []
    for member in bonding_list:
        if (member[::-1] not in new_list and len(member) == 2):
            new_list.append(member)
    
    return new_list
def main():
    #constructing_self_bonding_stats_table("/home/khanhvu/Desktop/PDBfile", "/home/khanhvu/Desktop/NearestChain.txt")
    #constructing_self_bonding_stats_table("/home/khanhvu/Desktop/PDBfile", "/home/khanhvu/Desktop/NearestChain.txt")
    
    df = pd.read_csv("/home/khanhvu/Desktop/Report/Analysis/intermolecular_bonding.csv")
    new_df = df.iloc[:,1:]
    print(new_df)
    new_df = add_percentage_for_DataFrame(new_df)
    new_df.to_csv("/home/khanhvu/Desktop/intermolecular_bonding.csv")
if __name__ == "__main__":
    main()

