# 23 mars 2024
# Aqel, Hamza 

# Le but de ce programme est de partir d’un brin d’ADN pour arriver aux 
# protéines codées par les gènes contenus dans ce brin d’ADN

from turtle import *

adn = "TCGACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGACG\
GCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGAGTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGAATGCCAGCCAGC\
CAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGAACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGA\
ACTCGACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCATCCCAGCGATACCC\
AGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAG\
CCAGCGAACTCGTCTGCGTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGC\
GATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGTATGCCAGCC\
AGCATCCCAGCGA"

codons_aa = {
    "UUU": "Phénylalanine",
    "UUC": "Phénylalanine",
    "UUA": "Leucine",
    "UUG": "Leucine",
    "CUU": "Leucine",
    "CUC": "Leucine",
    "CUA": "Leucine",
    "CUG": "Leucine",
    "AUU": "Isoleucine",
    "AUC": "Isoleucine",
    "AUA": "Isoleucine",
    "AUG": "Méthionine (Start)",
    "GUU": "Valine",
    "GUC": "Valine",
    "GUA": "Valine",
    "GUG": "Valine",
    "UCU": "Sérine",
    "UCC": "Sérine",
    "UCA": "Sérine",
    "UCG": "Sérine",
    "CCU": "Proline",
    "CCC": "Proline",
    "CCA": "Proline",
    "CCG": "Proline",
    "ACU": "Thrénine",
    "ACC": "Thrénine",
    "ACA": "Thrénine",
    "ACG": "Thrénine",
    "GCU": "Alanine",
    "GCC": "Alanine",
    "GCA": "Alanine",
    "GCG": "Alanine",
    "UAU": "Tyrosine",
    "UAC": "Tyrosine",
    "UAA": "Stop",
    "UAG": "Stop",
    "CAU": "Histidine",
    "CAC": "Histidine",
    "CAA": "Glutamine",
    "CAG": "Glutamine",
    "AAU": "Asparagine",
    "AAC": "Asparagine",
    "AAA": "Lysine",
    "AAG": "Lysine",
    "GAU": "Aspartate",
    "GAC": "Aspartate",
    "GAA": "Glutamate",
    "GAG": "Glutamate",
    "UGU": "Cystéine",
    "UGC": "Cystéine",
    "UGA": "Stop",
    "UGG": "Tryptophane",
    "CGU": "Arginine",
    "CGC": "Arginine",
    "CGA": "Arginine",
    "CGG": "Arginine",
    "AGU": "Sérine",
    "AGC": "Sérine",
    "AGA": "Arginine",
    "AGG": "Arginine",
    "GGU": "Glycine",
    "GGC": "Glycine",
    "GGA": "Glycine",
    "GGG": "Glycine"}

lettreAa = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "UGU": "C",
    "UGC": "C",
    "UGA": "*",
    "UGG": "W",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G"
}

# Utilisé pour la construction du brin d'ADN complementaire dans la fonction "antisens"
adnAppariement = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

# Utilisé pour faire la conversion d'ADN en ARN dans la fonction "transcrire"
adnTranscriptionDictionnaire = {
    "A" : "U",
    "T" : "A",
    "C" : "G",
    "G" : "C"
}

def antisens(brinAdn):
    # Part de brin ADN fourni et renvoie le brin d’ADN complémentaire.
    adnAntisens = ''
    for i in brinAdn:
        adnAntisens += adnAppariement[i]
    return adnAntisens

def testAntisens():
    # cas géneral:
    assert antisens("AAATTTCCCGGG") == "TTTAAAGGGCCC"



def reverseComplement(brinAdn):
    # Inversion du brin complémentaire d'ADN pour qu'il soit affiché dans le sense 5' -> 3'
    brinAdn = antisens(brinAdn) # Obtention du brin complémentaire d'ADN

    adnBrinReverseComplement = '' # Variable pour garder l'output

    for i in str(brinAdn[len(brinAdn)::-1]):
        adnBrinReverseComplement += i

    return adnBrinReverseComplement

def testReverseComplement():
    # cas géneral:
    assert reverseComplement("AAATTTCCCGGG") == "CCCGGGAAATTT"



def trouveDebut(brinAdn):
    # Recherche tous les codons de départ "TAC" sur un brin d’ADN et renvoie un tableau contenant
    # les positions du premier nucléotide de chacun des codons.

    resultat = [] # Garde les positions du premier nucléotide d'un start codon
    i = 0
    resultatSpecial = [] # Garde le tableau spécial avec TAC en position 3, 67, 89

    while i <= len(brinAdn)-len("TAC"):
        if "TAC" == brinAdn[i:i+len("TAC")]:
            resultat.append(i)
        i+=1

    return resultat



def testTrouveDebut():

    assert trouveDebut("-TAC--TAC---") == [1, 6]
    assert trouveDebut("---TAC------TAC----------------------------------------------------"
                       "TAC-------------------TAC") == [3, 12, 67, 89]



def trouveFin(brinAdn):
    # Même chose que la fonction précédente mais renvoie un tableau avec les positions de
    # tous les codons de terminaison "ATT", "ATC", "ACT" (attention, il y a trois possibilités de codons de
    # terminaison).

    resultat = []  # Garde les positions du premier nucléotide d'un stop codon
    i = 0
    

    while i <= len(brinAdn) - len("ATT"):
        if (   "ATT" == brinAdn[i:i + len("ATT")]
            or "ATC" == brinAdn[i:i + len("ATC")]
            or "ACT" == brinAdn[i:i + len("ACT")]   ):

            resultat.append(i)
        i += 1

    return resultat

def testTrouveFin():
    # 1 cas général pour le codon ATT:
    assert trouveFin("-ATT----") == [1]

    # 2 cas général pour le codon ATC:
    assert trouveFin("-ATT----ATC---") == [1, 8]

    # 3 cas général pour le codon ACT:
    assert trouveFin("-ATT----ATC---ACT---") == [1, 8, 14]



def trouveGene(debut, fin):
    
    # Prend en paramètre un tableau contenant les positions de tous les codons de départ et un
    # autre tableau contenant les positions de tous les codons de terminaison pour un brin
    # d’ADN et renvoie un tableau de tuples contenant la liste des gènes (début et fin) trouvés
    # sur un brin.
    
    # Ainsi, s’il y a trois gènes sur un brin, le tableau renvoyé ressemblera à :
    # [ (debutGene1, finGene1) , (debutGene2, finGene2) ,
    # (debutGene3, finGene3) ]
    
    # FinGene doit être supérieur à debutGene et finGene doit être situé à un multiple
    # de trois nucléotides de debutGene.
    
    # Hypothèse 1: Le nombre de génes ne peut pas être supérieur au nombre de codon "TAC" trouvés avec la fonction
    # trouveDebut(brinAdn). Si on trouve 3 codons "TAC" sur le brinAdn par exemple, on peut avoir un maximum
    # de 3 gènes
    # Hypothèse 2: Le chevauchement est acceptable uniquement pour les codons de fin

    listeDeGenes = []  # Garde les tuples avec les positions de début et fin d'un gène
    genesTemp = []  # Garde temporairemente les positions de genes qui vont être convertis en tuple
    
    # Cas spécial: si le paramètre debut contient une seule valeur
    if len(debut) == 1:
        max_j = 0
        for j in fin:
            if j > debut[0] and (j - debut[0]) % 3 == 0:
                max_j = j
                break
        genesTemp.append(debut[0])
        genesTemp.append(max_j)
        listeDeGenes.append(tuple(genesTemp))
        
    #Cas général:
    else:
        listeDeGenes = []
        genesTemp = []
        for i in debut:
            for j in fin:
                if j > i and (j - i) % 3 == 0:
                    genesTemp.append(i)
                    genesTemp.append(j)
                    listeDeGenes.append(tuple(genesTemp))
                    genesTemp = []
                    break
    return listeDeGenes

def testTrouveGene():
    #1 un gène avec un start codon et un stop codon:
    assert trouveGene([1],[4]) == [(1,4)]

    #2 deux gènes avec le même stop codon:
    assert trouveGene([1,4], [13]) == [(1,13), (4,13)]

    #3 un gène entre deux stop codons:
    assert trouveGene([5], [2, 11]) == [(5,11)]



def transcrire(brinAdn):
    
    # Prend en paramètre la sous-chaine de caractère du brin d’ADN débutant au début du gène
    # et se terminant à la fin du gène et renvoie le brin d’ARN correspondant sous forme d’une
    # chaine de caractères.
    
    arnTranscrit = ''
    for i in brinAdn:
        arnTranscrit += adnTranscriptionDictionnaire[i]
    return arnTranscrit

def testTranscrire():
    # cas géneral:
    assert transcrire("AAATTTCCCGGG") == "UUUAAAGGGCCC"



    
largeurCarre = 0 # variable globale utilisée dans plusieurs fonctions, 0 est une valeur initiale

def traduire(brinArn):
    
    # Prend en paramètre un brin d’ARN (chaine de caractères) et affiche la protéine sous
    # forme d’une chaine de caractères et la dessine à l’aide de la tortue.

    proteinString = 'Méthionine (Start)' # une protéine commence toujours par Méthionine (Start)
    proteinCorrige = 'M' # représente la protéine sous forme de lettres sans "-"

    # La boucle for est pour créer une chaine de caractères représentant 
    # la protéine (chaine d’acides aminés)
    for i in range(3, len(brinArn), 3): # La lecture se fait directement après le codon de début
        
        # ne pas afficher le stop
        if codons_aa[brinArn[i:i + 3]] == "Stop":
            pass
            
        # ajouter un "-" entre les acides juste après (Start)
        else:  
            proteinString += '-' + codons_aa[brinArn[i:i+3]]
    
    # On attribue une lettre du dictionnaire lettreArn à chaque acide aminé
    for i in range(3, len(brinArn), 3):
        
        if lettreAa[brinArn[i:i + 3]] == "*":
            pass
        
        else:
            proteinCorrige +=  lettreAa[brinArn[i:i+3]]


    # Dessin tortue: on dessine un carré correspondant à l'indice
    # n avec la fonction carre() et on écrit l'acide aminé sous forme 
    # d’une lettre avec la fonction milieuCarre()
    
    for n in range(len(proteinCorrige)):
        carre(20, n) # taille de chaque coté est 20, valeur suggerée pour une représentation optimale
        milieuCarre(n, proteinCorrige[n])


    return proteinString

def testTraduire():
    #1 cas général pour les codons codants:
    traduire("AUGCUU") == "Méthionine (Start)-Leucine"

    #2 séquences avec les trois stop codons:
    
    # séquence avec codon de fin "UGA":
    traduire("AUGCUUAUAUGAUAGUAA") == "Méthionine (Start)-Leucine-Isoleucine"

    # séquence avec codon de fin "UAG":
    traduire("AUGCUUAUAUAGUAGUAA") == "Méthionine (Start)-Leucine-Isoleucine"
    
    # séquence avec codon de fin "UAA":
    traduire("AUGCUUAUAUAAUAGUAA") == "Méthionine (Start)-Leucine-Isoleucine"
    


def carre(largeur, nombre):
    # Dessine un carré à l'indice "nombre"
    
    coordonnee_x = coordonne(nombre)[0]
    global largeurCarre
    largeurCarre = largeur
    coordonnee_y = coordonne(nombre)[1]
    
    positionner(coordonnee_x*largeur, coordonnee_y*largeur)
    for _ in range(4):
        fd(largeur); lt(90)
    positionner(-coordonnee_x * largeur, -coordonnee_y * largeur)

def positionner(coordonnee_x,coordonnee_y):
    # Cette fonction sert à positionner la tortue en fonction des parametres coordonne_x et coordonne_y
    pu(); fd(coordonnee_x); lt(270); fd(coordonnee_y); rt(270); pd()

def coordonne(nombre):
    # Cette fonction prend en parametre l'indice du carre à dessiner et retourne
    # les coordonnées x et y, exemples: indice=0 correspond au carre à 
    # la position [1,1]; indice=14 correspond au carre à la position [15,1]
    # le nombre maximal de colonnes est 15 (x_max=15), y varie en fonction de 
    # la largeur de la protéine
   
    coordonnee_x = nombre % 15 + 1
    coordonnee_y = nombre // 15 + 1
    
    return [coordonnee_x,coordonnee_y]

def testCoordonne():
    
    assert coordonne(0) == [1,1]
    
    assert coordonne(14) == [15,1]
    
    assert coordonne(15) == [1,2]
    
def milieuCarre(nombre, caractere):
    # Cette fonction sert à positionner la tortue au milieu du carre en fonction
    # de l'indice du carre et ecrit le caractere au milieu
    coordonnee_x = coordonne(nombre)[0]
    coordonnee_y = coordonne(nombre)[1]
    
    # positionner la tortue au centre du carre
    positionner(coordonnee_x * largeurCarre + largeurCarre/2, coordonnee_y * largeurCarre - largeurCarre/2)
    
    write(caractere)
    
    # retour à la position initale afin d'effectuer d'autres dessins
    positionner(-coordonnee_x * largeurCarre - largeurCarre/2, -coordonnee_y * largeurCarre + largeurCarre/2)


# Appel de fonctions test:

testAntisens()    
testReverseComplement()    
testTrouveDebut()
testTrouveFin()
testTrouveGene()
testTranscrire()    
testTraduire()
testCoordonne()

print("Les protéines codées par les gènes contenus dans ce brin d’ADN sont : \n")

# ------------------------------------------------------------------------
# test pour le brin adn

variable_fct = adn 
liste_genes = trouveGene(trouveDebut(variable_fct), trouveFin(variable_fct))
clear(800,600)
positionner(-100,-300)

for i in range(len(liste_genes)):
    
    gene = list(liste_genes[i])
    
    largeurGene = int((gene[1] - gene[0] + 3)/3) # nombre acides aminees
    
    print(traduire(transcrire(variable_fct[gene[0]:gene[1] + 3])))
    
    print('\n')
    
    # trois lignes suivantes servent à faire retour à la ligne et 
    # saut de ligne entre chaque protéine
    x = 0
    y = coordonne(largeurGene)[1] + 1
    positionner(x * largeurCarre, y * largeurCarre)

# -------------------------------------------------------------------------
# test pour le brin d’ADN complémentaire

variable_fct = reverseComplement(adn) 
liste_genes = trouveGene(trouveDebut(variable_fct), trouveFin(variable_fct))

for i in range(len(liste_genes)):
    
    gene = list(liste_genes[i])
   
    largeurGene = int((gene[1] - gene[0] + 3)/3) # nombre acides aminees
    
    print(traduire(transcrire(variable_fct[gene[0]:gene[1] + 3])))
    
    # trois lignes suivantes servent à faire retour à la ligne et 
    # saut de ligne entre chaque protéine
    x = 0
    y = coordonne(largeurGene)[1] + 1
    positionner(x * largeurCarre, y * largeurCarre)
