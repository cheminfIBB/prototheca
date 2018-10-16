from __future__ import unicode_literals
from copy import deepcopy
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-f" ,"--fasta" , nargs = "+" , help = "FASTA files with organisms proteoms.")
parser.add_argument("-r" ,"--rep_names_faa" , nargs = '+' ,help = "Representative names_faa of examined organisms.")
parser.add_argument("-db", "--databases", nargs = '*', help = "If analyzis type is base than here goes FASTA files of databases.")
parser.add_argument("-db_names", "--databases_names", nargs = '*', help = "Representative names of databases - these will be shown ")
parser.add_argument("-pc", "--prot_count", nargs = '*', help = "Number of proteins in examined proteomes.")
parser.add_argument("-e", "--e_cutoff", nargs = '?',default = 10, help = "Number of proteins in examined proteomes.")
args = parser.parse_args()

#nazwy "reprezentatywne" organizmów i baz
rep_names = args.rep_names_faa
base_names = args.databases_names
#nazwy plików FASTA organizmów i baz
org_fasta = args.fasta
db_names = args.databases
cutoff = float(args.e_cutoff)
#Rozmiar proteomów organizmów 
whole_protein_count = args.prot_count
#słownik "kodów" białek organizmów - napotykając białko w pliku wynikowym z blasta program dzięki temu słownikowi
#potrafi przyporządkować białko do konkretnego organizmu

#wczytuje kody białek organizmów
results = open("codes","r")
codes_org = dict()
for line in results:
    if line.strip().strip(">") in rep_names:
        current_org = line.strip().strip(">")
        codes_org[current_org] = []
    else:
        codes_org[current_org].append(line.strip())
results.close()

#wczytuje kody białek z baz
results = open("codes_bases","r")
codes = dict()
for line in results:
    if line.strip().strip(">") in base_names:
        current_base = line.strip().strip(">")
        codes[current_base] = []
    else:
        codes[current_base].append(line.strip().strip(">"))
results.close()

names = []
for base_name in base_names:
    for org_name in org_fasta:
        file_name = org_name.split(".")[0]
        file_name = file_name + "." + str(base_name).lower()
        names.append(file_name)

#lista zawierająca obiekty klasy organizm - czyli obrobione wyniki z blasta
organizmy = []

#lista zawierająca zbior id białek unikalnych dla każdego organizmu
unikalne = []

#w późniejszym etapie zawiera listę wygenerowaną funkcją make_merge na podstawie listy organizmy
#zawiera wszystkie możliwe kombinacje części wspólnych zbiorów kodów z każdego organizmu
merged = []
       
class Protein:
    def __init__(self,ids):
        self.id = ids       #identyfikator białka (jego kod z pliku blast)
        self.hits = []      #lista kodów białek które blast zaliczył jako hit
        self.e_val = []     #lista watości e-value kolejnych hitów
        
        #słownik zawierający pary nazwa rganizmu : lista kodów białek z tego organizmu które były hitami 
        self.taxa = {}
        for base_name in base_names:
            self.taxa[base_name] = []
        for p2 in codes_org:
            for l in codes_org[p2]:
                if self.id[0:len(l)] == l:
                    self.organism = p2

    def __repr__(self):
        return self.id
        
    #dodawanie wartości do list białka
    def add_hit(self,hit):
        self.hits.append(hit)
        
    def add_eval(self,val):
        self.e_val.append(val)
        
    #funkcja która usuwa informacje o hitach z e-value gorszym niż wskazane przez nas
    def eval_check(self):
        if len(self.e_val[-1].split("-")) > 1 and int(self.e_val[-1].split("-")[1]) >= int(cutoff): #ta linijka wyznacza minimalną dopuszczalną wartosć ujemnej potęgi e-value jakie zaakceptujemy
            for p2 in codes:
                for l in codes[p2]:
                    if self.hits[-1][0:len(l)] == l:
                        self.taxa[p2].append(self.id)

        elif len(self.e_val[-1].split("."))>1 and int(self.e_val[-1].split(".")[1]) == 0 and int(self.e_val[-1].split(".")[0]) == 0:
            for p2 in codes:
                    for l in codes[p2]:
                        if self.hits[-1][0:len(l)] == l:
                            self.taxa[p2].append(self.id)
        else:   #jeżeli hit ma zbyt niskie e-value to jest usuwany
            del(self.hits[-1])
            del(self.e_val[-1])           
         
#uogólnienie klasy Protein - przechowuję słownik zawierający:
#sumę analogicznych słowników ze wszystkich białek aktualnie badanego organizmu
#pozwala sprawdzić co właściwie dał nam blast - ile łącznie hitów było w pozostałych organizmach i jak się one mają do siebie
#można wyliczyć części wspólne
#pracuje na liście obiektów klasy Protein
class Organism:
    def __init__(self,prot_list):   
        self.id = prot_list[0].organism
        self.all_hits = prot_list[0].taxa
        for prot in prot_list:
            for each in prot.taxa:
                self.all_hits[each] = self.all_hits[each] + prot.taxa[each]
    def __repr__(self):
        return ""

#kod czytający pliki blastowe i tworzący na ich podstawie obiekty klasy Protein          
for i in range(0,len(names)):
    blast_results = []
    name = names[i] #aktualnie czytany plik
    print(name)
    rep_name = rep_names[i%len(rep_names)] #nazwa organizmu dla którego czytamy wyniki blasta
    file_name = open(name, "r")
    for k in file_name:
        if k[0:6] == "Query=":  #rozpoznaje wystąpienie nowego zapytania
            row = k.strip().split()
            blast_results.append(Protein(row[1]))
        if k[0:2] == "> ":  #wykrywa pojawienie się nowego hitu
            row = k.strip().split()
            blast_results[-1].add_hit(row[1])
        if k[0:6] == " Score" and len(blast_results[-1].hits)>len(blast_results[-1].e_val): #spisuje wartość e-value ostatniego hitu
            row = k.strip().split()
            blast_results[-1].add_eval(row[7].strip(","))
            blast_results[-1].eval_check() #sprawdza czy ostatnio dodany hit spełnia próg e-value narzucony przez funkcję eval_check
            
        #if k[0:6] == " Ident"  and len(blast_results[-1].hits)>len(blast_results[-1].identities): #spisuje wartość identyczności pomiędzy zapytaniem a ostatnim hitem
            #row = k.strip().split()
            #blast_results[-1].add_ident(row[2])
    file_name.close()
    org = Organism(blast_results) #tworzy obiekt klasy Organism na podstawie otrzymanej wyżej listy białek
    organizmy.append(org)   #zapisuje obiekt Organism w liście globalnej
    
for i in organizmy:
    for i2 in i.all_hits:
        i.all_hits[i2] = set(i.all_hits[i2])

organizmy2 = deepcopy(organizmy)

for i in organizmy2:
    for k, v in list(i.all_hits.items()):
         if len(v) == 0:
             del i.all_hits[k]

plt.rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

for i in rep_names:
    i = "\textit{0}".format(i)
    print(i)
sns.set_style("whitegrid")
dicti = dict()
dicti['whole'] = whole_protein_count
dicti2 = dict()
dicti['proteomes'] = rep_names
dicti2['proteomes'] = rep_names

wyniki = []     
temp = []
s = 0
current = ''
for i in organizmy2:
    #ilosc hitow z danej bazy
    l = 0
    if len(i.all_hits) == 1:     
        temp.append(len(i.all_hits[list(i.all_hits.keys())[0]]))
        current = list(i.all_hits.keys())[0]
    else:
        temp.append(l)
    if len(temp) == len(rep_names):
        sns.set()
        plt.figure(s)
        dicti2['blasthits'] = temp
        df = pd.DataFrame(dicti)
        df2 = pd.DataFrame(dicti2)
        sns.set(font_scale = 2.4)
        ax = sns.barplot(x=rep_names,y='blasthits',color="blue",data = df2)
        ax.set_ylabel("Number of organism proteins with hits from " + current,fontsize=24, alpha=0.8)
        ss = 0
        for p in ax.patches[0:int((len(ax.patches)))]:
            height = p.get_height()
            ax.text(p.get_x()+p.get_width()/2.,height+3,str(temp[ss]),ha="center")
            ss+=1      
        s+=1
        wyniki = wyniki + temp
        temp = []
        plt.gcf().subplots_adjust(bottom=0.1)
        #plt.savefig(os.path.dirname(os.path.abspath(__file__)) + "\\" + str(current).lower() + '.pdf')
        plt.show()
        plt.close()     
            
results = open('wyniki.txt','w')
l=0
s=0
for k in wyniki:
    if l == 0:
        results.write("> " + base_names[s] + " \n")
        s+=1
    results.write(':: ' + str(k) + ' \n')
    l+=1
    if l == len(rep_names):
        l = 0
results.close()