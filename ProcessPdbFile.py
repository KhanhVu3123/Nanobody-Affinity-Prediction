import FindingCDRRegion
import gzip
import os

amino_acid_dict = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLN': 'Q',
    'GLU': 'E',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
    'UNK': ''
}

def clean(string):
    new_string = ""
    for char in string:
        if(char != " "):
            new_string = new_string + char
    return new_string

def retrieve_sequence_from_zippedFolder(gzipFolder):
    """
    

    Parameters
    ----------
    gzipFolder : Your folder to all the zipped PDB file you want to analyze
                with correct file name format
    Returns
    -------
    The dictionary which the key is the PDB ID and the value being the list of
    2 sequences (1 antigen and 1 nanobody)

    """
    
    pdb_gz_files = [file for file in os.listdir(gzipFolder) if file.endswith(".pdb.gz")]
    
    final_dictionary = {}
    for file in pdb_gz_files:
        name = file[:4]
        sequences = retrieve_sequence_from_zippedFile(file)
        final_dictionary[name] = sequences
        
    return final_dictionary

def retrieve_nanobodySeq_fromzippedFile(gzipFile):
    pdb_gz_file_name = gzipFile[-14:]
    nanobody_chain = pdb_gz_file_name[5]
    seq = ""
    with gzip.open(gzipFile, 'rt') as file:
        for line in file:
            this_type = line[:6]
            this_type = clean(this_type)
            
            atom_name = line[13:17]
            atom_name = clean(atom_name)
            
            chain_name = line[21]
            
            if(this_type == "HETATOM"):
                break
            if(this_type == "ATOM"): 
                if(chain_name == nanobody_chain and atom_name == "CA"):
                    amino_acid = line[17:20]
                    seq = seq+ amino_acid_dict[amino_acid]
    return seq

def retrieve_antigenSeq_fromzippedFile(gzipFile):
    pdb_gz_file_name = gzipFile[-14:]
    antigen_chain = pdb_gz_file_name[6]
    seq = ""
    with gzip.open(gzipFile, 'rt') as file:
        for line in file:
            this_type = line[:6]
            this_type = clean(this_type)
            
            atom_name = line[13:17]
            atom_name = clean(atom_name)
            
            chain_name = line[21]
            
            if(this_type == "HETATOM"):
                break
            
            if(this_type == "ATOM" and atom_name == "CA"):
                if(chain_name == antigen_chain):
                    amino_acid = line[17:20]
                    seq = seq+ amino_acid_dict[amino_acid]
    return seq


def retrieve_all_nanobodySeq_fromzippedFile(gzipFolder):
    """
    

    Parameters
    ----------
    gzipFolder : The base folder to all your zipped PDB file

    Returns
    -------
     A dictionary that has PDB ID as the key as the nanobody sequence
     as the value

    """
    pdb_gz_files = [file for file in os.listdir(gzipFolder) if file.endswith(".pdb.gz")]
    
    final_dictionary = {}
    for file in pdb_gz_files:
        name = file[:4]
        nanobody_seq = retrieve_nanobodySeq_fromzippedFile(os.path.join(gzipFolder, file))
        final_dictionary[name] = nanobody_seq
        
    return final_dictionary


def retrieve_all_antigenSeq_fromzippedFile(gzipFolder):
    """
    

    Parameters
    ----------
    gzipFolder : The base folder to all your zipped PDB file

    Returns
    -------
     A dictionary that has PDB ID as the key as the nanobody sequence
     as the value

    """
    pdb_gz_files = [file for file in os.listdir(gzipFolder) if file.endswith(".pdb.gz")]
    
    final_dictionary = {}
    for file in pdb_gz_files:
        name = file[:4]
        antigen_seq = retrieve_antigenSeq_fromzippedFile(os.path.join(gzipFolder, file))
        final_dictionary[name] = antigen_seq
        
    return final_dictionary



def retrieve_sequence_from_zippedFile(gzipfile):
    sequences = []
    seq = ""
    prev_chain = ""
    with gzip.open(gzipfile, "rt") as file:
        for line in file:
            this_type = line[:6]
            this_type = clean(this_type)
            
            atom_name = line[13:17]
            atom_name = clean(atom_name)
            
            chain_name = line[21]
            
            
            if(chain_name != prev_chain):
                if(seq != ""):
                    sequences.append(seq)
                    seq = ""
            if(this_type== "HETATOM"):
                break
            if (this_type == "TER"):
                sequences.append(seq)
                seq = ""
            if(atom_name != "CA"):
                continue
            
            if(this_type != "ATOM"):
                continue
            
        
            amino_acid = line[17:20]
            seq = seq + amino_acid_dict[amino_acid]
            prev_chain = chain_name 
        if(seq != ""):
            sequences.append(seq)

    return sequences

def retrieve_sequence_of_chain(chain, filename):
    seq = ""
    with open(filename, "r") as file:
        for line in file:
            this_type = line[:6]
            this_type = clean(this_type)
            
            atom_name = line[13:17]
            atom_name = clean(atom_name)

            # Line 21 indicates the chain name of the sequence
            if(line[21] != chain):
                continue
            if(this_type== "HETATOM"):
                break
            if (this_type == "TER"):
                break
                seq = ""
            if(atom_name != "CA"):
                continue
            if(this_type != "ATOM"):
                continue
            amino_acid = line[17:20]
            seq = seq + amino_acid_dict[amino_acid]
    return seq

def retrieve_sequence(file1):
    sequences = []
    seq = ""
    with open(file1, 'r') as file:
        for line in file:
            this_type = line[:6]
            this_type = clean(this_type)
            
            atom_name = line[13:17]
            atom_name = clean(atom_name)
            
            if(this_type== "HETATOM"):
                break
            if (this_type == "TER"):
                sequences.append(seq)
                seq = ""
            if(atom_name != "CA"):
                continue
            if(this_type != "ATOM"):
                continue
            amino_acid = line[17:20]
            seq = seq + amino_acid_dict[amino_acid]
    return sequences
            
def coordinate_distance(lst1, lst2):
    """
    

    Parameters
    ----------
    
    lst1 : first list of 3 element of the coordinates
    lst2 : second list of 3 element of the coordinates, length equal to 
    lst1

    Returns
    -------
    Distance of the one amino acid that has coordinates of lst1 to the
    amino acid that has the coordinates of lst2.

    """
    
    for i in range(len(lst1)):
        lst1[i] = float(lst1[i])
        lst2[i] = float(lst2[i])
    distance = 0
    
    if(len(lst1) != len(lst2)):
        return "The two lists must have same number of elements"
    
    for i in range(len(lst1)):
        distance = distance + (lst1[i] - lst2[i]) **2
    return distance




def get_sequence(file1):
    """

    File1 is the PDB file you want to investigate. This funcion will
    return a dictionary, with each key being the chain name in that file
    The key value of each key will be a list of dictionaries, each has the key of 
    the amino acid and the value being the coordinates of the c-alpha of that 
    amino acid.
    
    """
    
    big_dict = {}
    with open(file1, "r") as file:

        all_pos_lst =[]

        for line in file:
            aminoDict = {}
            pos_lst = []

            this_type = line[:6]
            this_type = clean(this_type)

            atom_name = line[13:17]
            atom_name = clean(atom_name)

            if(this_type == "HETATOM"):
                break
            if(this_type == "TER"):
                big_dict[line[21]] = all_pos_lst
                all_pos_lst = []
                continue
            if(this_type != "ATOM" ):
                continue
            
            if(atom_name != "CA"):
                continue

            x_coordinate = line[31:38]
            x_coordinate = clean(x_coordinate)
            pos_lst.append(x_coordinate)

            y_coordinate = line[39:46]
            y_coordinate = clean(y_coordinate)
            pos_lst.append(y_coordinate)

            z_coordinate = line[47:54]
            z_coordinate = clean(z_coordinate)
            pos_lst.append(z_coordinate)

            amino_acid = line[17:21]
            amino_acid = clean(amino_acid)
            if(len(amino_acid) == 4):
                amino_acid = amino_acid[1:]
            if(amino_acid) not in amino_acid_dict.keys():
                continue
            amino_acid = amino_acid_dict[amino_acid]
            aminoDict[amino_acid] = pos_lst
            all_pos_lst.append(aminoDict)
        
    
    """
    Example output: {"A":[{"A":[1,2,3]}, {"B":[4,5,6]}], "B": [{"A":[1,2,3]}, {"B":[4,5,6]}]}
    """
    return big_dict

def to_float(lst):
    """
        The function takes a list of numerical string and return the list of 
        float variables.
    """
    new_lst = []
    for num in lst:
        new_num = float(num)
        new_lst.append(new_num)
    return new_lst

def determine_chain(chain, file1, distance, CDRlst = None):
    big_dict = get_sequence(file1)
    
    if(chain not in big_dict.keys()):
        return "Check your data"
    aminolst = big_dict[chain]
    sequence = ""
    for d in aminolst:
        sequence = sequence + list(d.keys())[0]

    
    list_of_dict = big_dict[chain]
    text = "".join([next(iter(d)) for d in list_of_dict]) 
    if CDRlst == None:
        lst = FindingCDRRegion.find_CDR(sequence)
    else:
        lst = []
        for cdr in CDRlst:
            index = text.find(cdr)
            lst.append(index)
            lst.append(index + len(cdr) - 1)
    aminolst = aminolst[lst[0] : lst[1] +1] + aminolst[lst[2] : lst[3] +1] + aminolst[lst[4] : lst[5] +1]


    distance_dict = {}
    # Aminolst is a list of dictionaries of coordiantes of every amino acids
    if(len(big_dict) <=1):
        return "Your dictionary is either wrong or you just have 1 chain"
    for chain_name in big_dict.keys():
        if chain_name == chain:
            continue
        distance_dict[chain_name] = 0
        aminolst2 = big_dict[chain_name]
        # AminoDict is a list, not a dictionary

        # Loop through every dictionary of the list.
        for amino_acid_coordinate in aminolst:
            # This is a dictionary
            lst1 = amino_acid_coordinate[next(iter(amino_acid_coordinate))]

            lst1 = to_float(lst1)
            # Change all of its elements to integers.

            for amino_acid_coordinate2 in aminolst2:
                lst2 = amino_acid_coordinate2[next(iter(amino_acid_coordinate2))]
                lst2 = to_float(lst2)
                distance_squared = 0
                for i in range(3):
                    distance_squared = distance_squared+ ((lst1[i] - lst2[i]) ** 2)
                if(distance_squared < distance ** 2):
                    distance_dict[chain_name] +=1
    return distance_dict

def main():

    sequences = retrieve_sequence("/home/khanhvu/Desktop/PDBfile/4LSP.pdb")
    print(sequences)
if __name__ == "__main__":
    main()
    
