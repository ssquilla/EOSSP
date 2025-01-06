ssquilla, Samuel Squillaci, samuel.squillaci@mailfence.com

**REMARQUES GENERALES**

Ce projet contient les algorithmes développés à l'Onera dans le cadre de la thèse de Samuel Squillaci (https://theses.fr/2023ESAE0071), présentant les méthodes publiées notamment dans les articles suivants :

- Squillaci, Samuel et Roussel, Stéphanie et Pralet, Cédric : «Parallel Scheduling of Complex Requests  for a Constellation of Earth Observing Satellites». Dans : Passerini, A., Schiex, T. (eds.). Frontiers in Artificial Intelligence and Applications. IOS Press 2022.
- Squillaci, Samuel et Roussel, Stéphanie et Pralet, Cédric : «Scheduling Complex Observation Requests for a Constellation of Satellites: Large Neighborhood Search Approaches». Dans : CPAIOR 2023. https://link.springer.com/chapter/10.1007/978-3-031-33271-5_29
- Squillaci, Samuel et Roussel, Stéphanie et Pralet, Cédric : «Comparison of time-dependent and time-independent scheduling approaches for a constellation of Earth observing satellites». Dans : IWPSS 2023. https://drive.com/file/d/1D_5v45e4jTcrSGwc9cQkHLRqb8xF0Tz2/view (p.96)

L'intégralité des algorithmes sont regroupés au sein d'un projet Python présentant diverses options. Le code fourni peut présenter des résultats légèrement différents des articles pour les raisons suivantes :
- la capacité de calcul des machines de l'Onera présentent des caractéristiques particulières (20 coeurs de calcul)
- certaines options ne sont pas fournies à l'heure actuelle (y compris le solver externe ImaxLNS* qui permet d'accélérer notamment le LNS glouton et les algorithmes gloutons stochastiques itérés)
- si vous ne disposez pas d'une version complète de CPlex (académique ou commerciale), la taille des problèmes à résoudre est limité. L'algorithme exacte CPSolver ne renverra pas de solution dans ce cas, et le CPLNS ne trouvera pas de solution pour l'exploration des voisinages concernés
- un refactoring récent a été réalisé

*lien vers l'article introduisant ImaxLNS: 
Pralet, C.: Iterated maximum large neighborhood search for the traveling salesman problem with time windows and its time-dependent version. Computers & Operations Research 150 (2022)

Le solver LKH3 est exploité en tant que solver TSPTW externe. Le code provient du site suivant : http://webhotel4.ruc.dk/~keld/research/LKH-3/

**ENVIRONNEMENT**

Si conda n'est pas installé, installer conda. Sur linux, miniconda peut s'installer de la manière suivante :

    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm -rf ~/miniconda3/miniconda.sh
    ~/miniconda3/bin/conda init bash

Mettre en place l'environnement conda :

    conda env create -f envFile/envEOSSP.yml -n EOSSP
    conda activate EOSSP
    conda install python=3.6

Mettre en place le solver externe LKH :

    make download_lkh
    make all

L'environnement conda doit être activé à chaque redémarrage du pc, mais pas nécessairement recréé.

**NOTICE DES ALGO**

Pensez à :

    - Compiler les sources de LKH-3 : make all
    
    - Vider les répertoires tmp de LKH de temps en temps (contiennent des fichiers tmp si les algo se sont interrompus) : make clean
    
    - Aller dans le sous-repertoire EOSSP (contenant le main.py) : cd EOSSP
    
    - Utiliser la commande suivante pour exécuter le code : mpiexec -n <Ncoeurs> python3 main.py -v <indice de verbosité> -t <time> --solver=<le solver à utiliser>


/!\ Exécuter le code depuis un autre répertoire ne fonctionne pas. Python utilise le répertoire courant comme référence pour les chemins relatifs, et non le répertoire des sources.

**UTILISATION**

Pour connaître les options disponibles, tapez "python3 main.py -v 1 -h" ou "python3 main.py --solver=<algo> -v 1 --help".
La commande renverra un texte similaire au texte suivant, listant les options disponibles POUR L'ALGORITHME INDIQUE PAR --solver=... :

    --sample : Indiquer un chemin pour sauvegarder les résultats.
    -v || --verbose : Indiquer le niveau d'affichage.
    -o || --verif : Activer pour effectuer des vérifications en cours d'execution (plus lent, plus sûr).
    -l || --light_sample : Activer pour alléger la quantité d'information sauvegardée.
    -t || --time : Indiquer le temps d'execution (en minutes).
    -m || --modes : Indiquer le nombre de modes considérés par requêtes.
    -w || --threads : Indiquer le nombre de threads travailleurs.
    -h || --help : Afficher l'aide.
    -c || --comm : Ajouter un commentaire à l'échantillon de résultats.
    -f || --folder : Indiquer le dossier de l'instance.
    -n || --file : Indiquer le fichier de l'instance.
    -s || --seed : Graine aléatoire.
    --include_systematic : Inclure les requêtes systématiques.

/!\ Un niveau de verbosité non nul doit être indiqué afin d'afficher les informations dans la console au cours de l'éxécution ( -v <niveau de verbosité>). Par défaut, l'indice est nul. Par défaut, les requêtes systématiques ne sont pas prises en compte. Les inclure à l'aide de l'option "--include_systematic".

**SAUVEGARDE DES DONNEES**

La sauvegarde des données peut se faire à l'aide de l'option "--sample=<nom du dossier cible>" qui va créé un dossier si nécessaire et fabriquer un pickle Python dans celui-ci. Les résultats peuvent être analysés à partir du notebook "analyse_resultats", ou simplement à l'aide de la bibliothèque python "pickle" depuis un script python.
