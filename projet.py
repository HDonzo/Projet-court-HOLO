#projet court
import Bio
from Bio import AlignIO
import numpy as np
import scipy.cluster.hierarchy as sch
from Bio import SeqIO

#etape 1 : initialisation de la matrice blosum

def  read_blosum(matrix_file) : 
        blosum_dict = {}
        with open(matrix_file, 'r') as filin : 
            lines = filin.readlines()
        AA = lines[0].split()[0:]
        for line in lines[1:]:
            content = line.split()
            AA1 = content[0]
            
            for i, AA2 in enumerate(AA) : 
                key = (AA1,AA2)
                value = content[i+1]
                blosum_dict[key] = value
        return blosum_dict
        
        
#etape 2 : lire le fichier fasta 

def Read_fasta(file) : 
    numero_seq = 0
    dico = {'seq':[]}
    prot_dict = {}
    with open (file, 'r') as fasta_file :
        prot_id = ""

            
        for line in fasta_file : 
            if line.startswith(">"):
                numero_seq += 1
                prot_id = numero_seq
                prot_dict[prot_id] = ""
            else : 
                prot_dict[prot_id] += line.strip()
        for aa1 in prot_dict : 
            dico['seq'].append(prot_dict[aa1])
        return prot_dict      
 

             
#etape 3 : needleman et vunsch, matrice de score  
 
 
#creation d'une matrice de score entre 2 seq
 
def calculate_score_matrix(seq1, seq2, blosum_matrix): 
    
    score_matrix = []
    
    for res1 in seq1:
        row = []
        
        for res2 in seq2:
            pair = (res1, res2)
            
            if pair in blosum_matrix:
                row.append(int(blosum_matrix[pair]))
            else:
                row.append(0) #si on ne trouve pas la paire dans le dico on ajoute la valeur 0
                
        score_matrix.append(row)
        
    return np.array(score_matrix) 
 
# calcul de la matrice de distance à partir de la matrice de score

def calculate_distance_matrix(score_matrix): 
    
    max_score = max(max(row) for row in score_matrix)
    distance_matrix = []
    
    for row in score_matrix:
        distance_row = [1 - (score / max_score) for score in row]
        distance_matrix.append(distance_row)
     
    return np.array(distance_matrix)

#calcul de score et de distance optimisé pour un alignement par paire

def calcul_score_et_distance(sequence_dict, blosum_matrix): 
    
    liste_DM = []
    sequence_ids = list(sequence_dict.keys())
    
    for i in range(0,len(sequence_ids),2):
        seq1_ids = sequence_ids[i]
        seq1 = sequence_dict[seq1_ids]
        
        for j in range(i+1,len(sequence_ids)) : 
            seq2_ids = sequence_ids[j]
            seq2 = sequence_dict[seq2_ids] 
            
            score_matrix = calculate_score_matrix(seq1, seq2, blosum_matrix) 
            distance_matrix = calculate_distance_matrix(score_matrix) 
            
            liste_DM.append((seq1_ids,seq2_ids, distance_matrix))

            
            print(f"Matrice de score pour la sequence {seq1_ids} et {seq2_ids} : \n{score_matrix}")
            print(f"Matrice de distance pour la sequence {seq1_ids} et {seq2_ids} : \n{distance_matrix}")

            
            break
    with open('distance_liste.txt', 'w') as filout:
        for seq1_ids, seq2_ids, distance_matrix in liste_DM:
            filout.write(f"\n\n")
            for row in distance_matrix:
                filout.write(" ".join(map(str, row)) + "\n")
            filout.write("\n")
        print("sauvegarde okay")

#embranchement sequentiel

def embranchement_seq(file_distance):   
    
    with open(file_distance, 'r') as file:
        distance_matrices = []
        current_matrix = []
        
        for line in file:
            line = line.strip()
            
            if not line:

                if current_matrix:
                    distance_matrices.append(np.array(current_matrix))
                    current_matrix = []
            else:
                
                values = [float(value) for value in line.split()]
                current_matrix.append(values)

    
    if current_matrix:
        distance_matrices.append(np.array(current_matrix))

    sequence_ids = list(sequence_dict.keys())

    if distance_matrices : 
        linkage_matrix = sch.linkage(distance_matrices[0], method='single')
        
        dendrogram = sch.dendrogram(linkage_matrix, no_plot=True)
        
        nombre_de_clusters = 5
        
        clusters = sch.fcluster(linkage_matrix, t=nombre_de_clusters, criterion='maxclust')
        
        print("Clusters résultants :")
        
        for i in range(1, nombre_de_clusters + 1):
            cluster_i = [sequence_ids[j] for j in range(len(sequence_ids)) if clusters[j] == i]
            print(f"Cluster {i}: {cluster_i}")
        #sch.dendrogram(dendrogram) #n'arrive pas a s'afficher
    else:
        print(f"Aucune matrice de distance n'a été lue depuis le fichier {file_distance}.")
   

       
#main 

#lecture de la matrice blosum
matrix_file = "blosum62.txt"
print(read_blosum(matrix_file))

#Lecture du fichier FASTA
file = "test.fasta"
print(Read_fasta(file))

#Calcul de la matrice de score de deux séquences
sequence_dict = Read_fasta(file)
seq1 = sequence_dict[1]  # sequence 1 
seq2 = sequence_dict[2]  #sequence 2 
blosum_matrix = read_blosum(matrix_file)
print(calculate_score_matrix(seq1, seq2, blosum_matrix))

#Calcul de la matrice de distance de deux séquences 
score_matrix = calculate_score_matrix(seq1, seq2, blosum_matrix)
print(calculate_distance_matrix(score_matrix))

#calcul de la matrice de score et de distance pour un fichier contenant X sequences
print(calcul_score_et_distance(sequence_dict,blosum_matrix))

#embranchement sequentiel
file_distance = "distance_liste.txt"
print(embranchement_seq(file_distance))

#alignement multiple 
file = "test.fasta"

#le script d'alignement multiple n'a pas pu etre réalisé