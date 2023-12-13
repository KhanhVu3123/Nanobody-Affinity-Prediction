from Bio import pairwise2
from Bio.Seq import Seq

def edit_distance(x, y):
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    for i in range(len(x)+1):
        D[i][0] = i
    for j in range(len(y)+1):
        D[0][j] = j
    for i in range(1,len(x)+1):
        for j in range(1,len(y)+1):
            Dia = D[i-1][j-1]
            Ver = D[i][j-1]
            Hori = D[i-1][j]
            delt = 0
            if(x[i-1] != y[j-1]):
                delt = 1
            D[i][j] = min(Dia+delt,Ver+1,Hori+1)
    return D[-1][-1]


def delete_dash(string):
    new_string = ""
    for char in string:
        if(char != "-"):
            new_string = new_string + char
    return new_string

def split(string):
    lst = []
    this = string.split("-")
    for char in this:
        if(len(char)>0):
            lst.append(char)
    return lst

def compare(seq,frameworks):

    alignment = pairwise2.align.globalxx(seq,frameworks)
    if not alignment:
        return None
    min_dash = 100
    for align in alignment:
        string = align[1]
        num = len(split(string))
        if(num<min_dash):
            min_dash = num

    for align in alignment:
        if(len(split(align[1])) == min_dash):
            return align



# The function to return the index of the starting and ending index of 3 CDRs in the nanobody.
def find_CDR(string):
    framework1 = "QVQLVESGGGLVQAGGSLRLSCAASG"
    framework2 = "WYRQAPGKQRELVA"
    framework3 = "DSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYC"
    framework4 = "WGQGTQVTVSS"
    lst = []
    lst.append(framework1)
    lst.append(framework2)
    lst.append(framework3)
    lst.append(framework4)

    frameworks = framework1 +framework2 + framework3 + framework4
    string = Seq(string)
    align = compare(string,frameworks)
    if align == None:
        return "Check your data"
    seqB = align[1] 
    seqA = align[0]
    # lst is the list of the frameworks.
    FR_framework = []
    index_lst = []

    # Find the matching k-mer in the sequence.
    for element in lst:
        length_of_k_mer = len(element) + 5
        min_point = 100
        index = 0
        for i in range(len(seqB) - length_of_k_mer +1):
            k_mer = seqB[i:i+length_of_k_mer]
            k_mer = delete_dash(k_mer)
            point = edit_distance(k_mer,element)
            if(point < min_point):
                min_point = point
                index = i

    # Check 6 of each k-mer
        min_point = 100
        j_value = 0
        j_limit = min(6,len(seqB) - index-len(element))
        for j in range(j_limit):
            j_mer = seqB[index:index+len(element)+j]
            modified_j_mer = delete_dash(j_mer)
            point = edit_distance(modified_j_mer,element)
            if(point < min_point):
                min_point = point
                j_value = j
        

        forward_index = index
        while(seqB[forward_index] == "-"):
            forward_index = forward_index+1
        
        backward_index = index+len(element) + j_value
        while(seqB[backward_index] == "-"):
            backward_index = backward_index -1

        index_lst.append(forward_index)
        index_lst.append(backward_index)

        FR_framework.append(seqB[forward_index:backward_index+1])

    CDR_lst = []
    CDR_lst.append(seqA[index_lst[1]+1:index_lst[2]])
    CDR_lst.append(seqA[index_lst[3]+1:index_lst[4]])
    CDR_lst.append(seqA[index_lst[5]+1:index_lst[6]])

    index1 = string.index(delete_dash(CDR_lst[0]))
    index2 = string.index(delete_dash(CDR_lst[1]))
    index3 = string.index(delete_dash(CDR_lst[2]))

    new_index_lst = []
    new_index_lst.append(index1)
    new_index_lst.append(index1 + len(CDR_lst[0]))

    new_index_lst.append(index2)
    new_index_lst.append(index2 + len(CDR_lst[1]))

    new_index_lst.append(index3)
    new_index_lst.append(index3 + len(CDR_lst[2]))

    return new_index_lst

def main():

    #seq = "QVQLVESGGGLVQPGGSLTLSCTASGFTLDHYDIGWFRQAPGKEREGVSCINNSDDDTYYADSVKGRFTIFMNNAKDTVYLQMNSLKPEDTAIYYCAEARGCKRGRYEYDFWGQGTQVTVSSKKKHHHHHH"
    seq = input()
    print(find_CDR(seq))

if __name__ == "__main__":
    main()


