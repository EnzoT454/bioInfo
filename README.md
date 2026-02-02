# BioInfo — ADN → ARN → Protéines (avec visualisation Turtle)

Petit projet Python qui part d’un **brin d’ADN**, détecte les **gènes** (start/stop codons), transcrit en **ARN**, puis traduit en **protéines** (acides aminés).  
Le programme affiche les protéines trouvées et les **dessine sous forme de grille** (carrés + lettres) à l’aide de la tortue.

## Fonctionnalités

- Construction du **brin complémentaire** (antisens) et **reverse complement** (5' → 3')
- Détection des **codons de départ** (start) : `TAC`
- Détection des **codons stop** : `ATT`, `ATC`, `ACT`
- Extraction des **gènes** (start/stop valides, longueur multiple de 3)
- **Transcription** ADN → ARN (A→U, T→A, C→G, G→C)
- **Traduction** ARN → protéine (nom complet + code 1 lettre)
- **Visualisation** de la protéine avec Turtle (carrés + lettre au centre)
- Tests simples via `assert` (fonctions `test...()`)

## Prérequis

- Python 3.x
- Module Turtle (généralement inclus avec Python)

> Le programme utilise des fonctions Turtle comme `fd`, `lt`, `rt`, `pu`, `pd`, `write`, `clear`.  


## Exécution

Depuis la racine du projet :

```bash
python3 bioInfo.py
```

## Résultat attendu

### Affichage console

Les protéines codées par les gènes contenus dans ce brin d’ADN sont :


Puis une liste de protéines, par exemple :

Méthionine (Start)-Leucine-Isoleucine

### Visualisation

Ouverture d’une fenêtre **Turtle** qui dessine les protéines trouvées sous forme de carrés, avec une lettre représentant chaque acide aminé.

---

## Détails biologiques (simplifiés)

- **Start codon (ADN)** : `TAC`  
  (correspond au codon start `AUG` en ARN)

- **Stop codons (ADN)** : `ATT`, `ATC`, `ACT`

- Un **gène est valide si** :
  - `stop > start`
  - `(stop - start)` est un multiple de `3` (lecture par codons)

---

## Tests

Le script exécute automatiquement des tests à l’aide de `assert` :

- `testAntisens()`
- `testReverseComplement()`
- `testTrouveDebut()`
- `testTrouveFin()`
- `testTrouveGene()`
- `testTranscrire()`
- `testTraduire()`
- `testCoordonne()`

Si un test échoue, Python lèvera une erreur `AssertionError`.
















