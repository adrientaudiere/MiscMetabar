# Session Context

## User Prompts

### Prompt 1

Under windows I obtain the following error :

Barcodes : 
16S bactéries (341F-785R) avec SILVA 
18S archées (Arch349-arch806) avec SILVA + KSGP 
ITS2 champignons (ITS86-ITS4) avec EUKaryome + UNITE
18S microeucaryotes (TAReuk…) avec PR2 et EUKaryome
18S gloméromycètes (AMV4.5-AMDGR) avec PR2 et MAARJAM 

Méthode de regroupement : DADA2 et SWARM 
Algorithme d’assignation : RDP et Blast


I think it is a matter of "/" instead of "\". How to deal with this in all functions from MiscMetabar package

### Prompt 2

Erreur dans (function (physeq = NULL, ref_fasta = NULL, seq2search = NULL,  : 
  Vsearch sintax failed with status 1.


Fatal error: Unrecognized string on command line (-)
De plus : Message d'avis :
Dans system2(vsearchpath, args = cmd_sintax, stdout = TRUE, stderr = TRUE) :
  l'exécution de la commande '"C:\Users\2024cb004\AppData\Roaming/R/data/R/MiscMetabar/bin/vsearch.exe"  --sintax C:\Users\2024CB~1\AppData\Local\Temp\RtmpAdSmYa/temp.fasta --db C:/Users/2024cb004/Desktop/Post-doc/7 - RE...

### Prompt 3

roll the shQuote() fix across all the vsearch/blast/mmseqs2/dada2
  shell-out functions. Krona and cutadapt are Unix-only functions.

