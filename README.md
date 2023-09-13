# Projet-court-HOLO

# README - Projet Court

## Objectif

Ce projet court vise à réaliser différentes tâches liées à l'analyse de séquences protéique en utilisant les étapes de l'algorithme de Needleman et Wunsch. Les principales étapes de cet algorithme comprennent l'initialisation de la matrice BLOSUM, la lecture d'un fichier FASTA et le calcul de la matrice de score. Nous chercherons par la suite à représenter ces alignements sous forme d'arbre grace au calcul de la matrice de distance et l'embranchement séquentiel puis nous essaieront de réaliser un alignement multiple.

## Étapes du projet

### Étape 1 : Initialisation de la Matrice BLOSUM

La matrice BLOSUM (Block Substitution Matrix) est initialisée à partir d'un fichier de matrice spécifié. La fonction `read_blosum(matrix_file)` lit le fichier, extrait les valeurs de la matrice et les stocke dans un dictionnaire.

Commandes : 
matrix_file = 'blosum62.txt'

Lancer la fonction : 
print(read_blosum(matrix_file))

### Étape 2 : Lecture du Fichier FASTA

Le programme lit un fichier FASTA contenant les séquences à aligner par la suite. La fonction `Read_fasta(file)` parcourt le fichier, extrait les séquences et les stocke dans un dictionnaire qui renverra en sortie l'ensemble des séquences toutes précédés d'un numéro rendant la lecture plus simple.

Commandes : 
file = 'test.fasta'

Lancer la fonction : 
print(Read_fasta(file))


### Étape 3 : Calcul de la Matrice de Score

Le programme calcule la matrice de score entre deux séquences en utilisant la matrice BLOSUM précédemment initialisée et le dictionnaire contenant les séquences du fichier FASTA. La fonction `calculate_score_matrix(seq1, seq2, blosum_matrix)` prend deux séquences et la matrice BLOSUM comme entrées et renvoie la matrice de score.

Commandes : 
sequence_dict = Read_fasta(file)
seq1 = sequence_dict[1]  # sequence 1 
seq2 = sequence_dict[2]  #sequence 2
blosum_matrix = read_blosum(matrix_file)

Lancer la fonction : 
print(calculate_score_matrix(seq1, seq2, blosum_matrix))

### Étape 4 : Calcul de la Matrice de Distance

À partir de la matrice de score, le script calcule la matrice de distance. La fonction `calculate_distance_matrix(score_matrix)` va normaliser les valeurs de la matrice de score pour obtenir une matrice de distance.

Commandes:
score_matrix = calculate_score_matrix(seq1, seq2, blosum_matrix)

Lancer la fonction : 
print(calculate_distance_matrix(score_matrix))

### Étape 5 : Calcul de Score et de Distance Optimisé pour un Alignement par Paire

Ici, le programme effectue un alignement par paire entre plusieurs séquences. Les séquences sont alignées les unes après les autres, et pour chaque paire, une matrice de score et une matrice de distance sont calculées. Les résultats sont ensuite stockés dans un fichier de sortie.

Commandes : 
sequence_dict = Read_fasta(file)
blosum_matrix = blosum_matrix = read_blosum(matrix_file)

Lancer le programme : 
print(calcul_score_et_distance(sequence_dict,blosum_matrix))

### Étape 6 : Embranchement Séquentiel

Le programme lit les matrices de distance à partir du fichier de sortie généré à l'étape précédente. En utilisant ces matrices, il effectue un regroupement séquentiel en utilisant la méthode "single linkage" de la bibliothèque SciPy. Les clusters résultants sont affichés avec leur composition.

Commandes : 
file_distance = "distance_liste.txt"

Lancer le programme : 
print(embranchement_seq(file_distance))

### Étape 7 : Alignement Multiple

L'objectif était de réaliser un alignement multiple à l'aide des résultats du regroupement séquentiel, mais cette étape n'a malheureusement pas été achevée dans le script fourni.

## Conseil d'utilisation

1. Assurez-vous d'avoir bien installé les bibliothèques Biopython, NumPy et SciPy.
2. Préparez un fichier de matrice BLOSUM (par exemple ici "blosum62.txt") et un fichier FASTA contenant vos séquences (par exemple "test.fasta").
3. Avant d'exécuter le script faites bien attention de bien spécifier les fichiers de matrice BLOSUM et FASTA.
4. Veuillez bien suivre les étapes du projet dans l'ordre annoncé.

## Remarque

- La modélisation de l'arbre après l'embranchement séquentiel n'a pas pu etre affiché le terminal ubuntu. Nous vous conseillons d'utiliser le package matplotlib si l'arbre ne s'affiche pas non plus sur votre terminal.
- L'alignement multiple prévu dans le script n'a pas pu etre réalisé.
- Assurez-vous que les dépendances requises sont installées pour exécuter le script avec succès.

