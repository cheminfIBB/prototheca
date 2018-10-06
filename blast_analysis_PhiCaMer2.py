from __future__ import unicode_literals
# -*- coding: utf-8 -*-
import numpy as np
from copy import deepcopy
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

#słownik "kodów" białek organizmów - napotykając białko w pliku wynikowym z blasta program dzięki temu słownikowi
#potrafi przyporządkować białko do konkretnego organizmu

#wczytuje kody z pliku
results = open("codes","r")
codes = dict()
for line in results:
    if line[0] != '>':
        current_base = line[0:-2]
        codes[current_base] = []
    else:
        codes[current_base].append(line[1:6])
results.close()

codes_org = dict()
codes_org["Chlorella variabilis"] = ["XP_005"]
codes_org['Auxenochlorella protothecoides'] = ['XP_011']
codes_org['Candida albicans'] = ['XP_436','XP_713','XP_019','XP_711','XP_710','XP_714','XP_717','XP_722','XP_709','XP_715','XP_721','XP_712','XP_443','XP_719','XP_441','XP_716','XP_718','XP_720','XP_723','XP_442']
codes_org['Helicosporidium sp.'] = ['KDD']
codes_org['Trichophyton rubrum'] = ['XP_003']

#nazwy plików wynikowych blast z bazą phiscan
names_phi = ["prototheca.phiscan"]
#, "chlorella_variabilis.phiscan","auxenochlorella_protothecoides.phiscan","helicosporidium_sp.phiscan","trichophyton_rubrum.phiscan","candida_albicans.phiscan"]

#nazwy plików wynikowych blast z bazą merops
names_merops = ["prototheca.merops"]
#,"chlorella_variabilis.merops","auxenochlorella_protothecoides.merops","helicosporidium_sp.merops","trichophyton_rubrum.merops","candida_albicans.merops"]

#nazwy plików wynikowych blast z bazą cazy
names_cazy = ["prototheca.cazy"]
#,"chlorella_variabilis.cazy","auxenochlorella_protothecoides.cazy","helicosporidium_sp.cazy","trichophyton_rubrum.cazy","candida_albicans.cazy"]

names = names_phi+names_merops+names_cazy

#Rozmiar proteomów organizmów (w moim przypadku z ncbi genome) w kolejności: Chlorella, Auxenochlorella, Helicosporidum,Trichophyton, Candida
whole_protein_count = [5930.9780, 7014, 6033, 8706, 6030]

#Ilość sekwencji w bazach w kolejności: PhiScan,Merops,Cazy
whole_base_count = [4003,4967,921174]

#nazwy plikow baz danych w formacie fasta
db_names = ["phi_accessions.fasta","merops_scan.fasta","CAZyDB.fasta"]

#nazwy "reprezentatywne" organizmów
rep_names = ["Prototheca","Chlorella variabilis","Auxenochlorella protothecoides","Helicosporidium sp.","Trichophyton rubrum","Candida albicans"]
base_names = ["PhiScan","MEROPS","CAZY"]

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
        self.identities = []    #lista ułamków obrazujących jak podobne jest białko do hitu pod danym indeksem
        self.e_val = []     #lista watości e-value kolejnych hitów
        
        #słownik zawierający pary nazwa rganizmu : lista kodów białek z tego organizmu które były hitami 
        self.taxa = {"PhiScan": [], "MEROPS": [], "CAZY": []}
        for p2 in codes_org:
            for l in codes_org[p2]:
                if self.id[0:len(l)] == l:
                    self.organism = p2

    def __repr__(self):
        return self.id
        
    #dodawanie wartości do list białka
    def add_hit(self,hit):
        self.hits.append(hit)

    def add_ident(self,ident):
        self.identities.append(ident)

    def add_eval(self,val):
        self.e_val.append(val)
        
    #funkcja która usuwa informacje o hitach z e-value gorszym niż wskazane przez nas
    def eval_check(self):
        if len(self.e_val[-1].split("-")) > 1 and int(self.e_val[-1].split("-")[1]) >= 10: #ta linijka wznacza minimalną dopuszczalną wartosć ujemnej potęgi e-value jakie zaakceptujemy
            for p2 in codes:
                for l in codes[p2]:
                    if self.hits[-1][0:len(l)] == l:
                        #kolejne dwie linijki warunkują czy przechowujemy swoje id czy id hitów w słowniku zawierającym organizmy
                        #to znaczy albo wiemy, że to białko ma hit z tego organizmu, albo jakie białka z innych organizmów miały hit do tego białka 
                        self.taxa[p2].append(self.hits[-1])
                        #self.taxa[p2].append(self.id)

        elif len(self.e_val[-1].split("."))>1 and int(self.e_val[-1].split(".")[1]) == 0 and int(self.e_val[-1].split(".")[0]) == 0:
            for p2 in codes:
                    for l in codes[p2]:
                        if self.hits[-1][0:len(l)] == l:
                            self.taxa[p2].append(self.hits[-1])
                            #self.taxa[p2].append(self.id)
        else:   #jeżeli hit ma zbyt niskie e-value to jest usuwany
            del(self.hits[-1])
            del(self.identities[-1])
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
        print(self.id)
        return ""

#kod czytający pliki blastowe i tworzący na ich podstawie obiekty klasy Protein          
for i in range(0,len(names)):
    blast_results = []
    name = names[i] #aktualnie czytany plik
    print(name)
    rep_name = rep_names[i%len(rep_names)] #nazwa organizmu dla którego czytamy wyniki blasta
    file_name = open(name, "r")
    p=0
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
        if k[0:6] == " Ident"  and len(blast_results[-1].hits)>len(blast_results[-1].identities): #spisuje wartość identyczności pomiędzy zapytaniem a ostatnim hitem
            row = k.strip().split()
            blast_results[-1].add_ident(row[2])
            blast_results[-1].eval_check() #sprawdza czy ostatnio dodany hit spełnia próg e-value narzucony przez funkcję eval_check
        p+=1
        if p%10000 == 0:
            print(p)
    org = Organism(blast_results) #tworzy obiekt klasy Organism na podstawie otrzymanej wyżej listy białek
    if name == "candida_albicans.phiscan":
        print(org.all_hits)
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
rep_names = ["\it{Chlorella variabilis}","\it{Auxenochlorella protothecoides}","\it{Helicosporidium sp.}","\it{Trichophyton rubrum}","\it{Candida albicans}"]
sns.set_style("whitegrid")
dicti = dict()
dicti['whole'] = whole_protein_count
dicti2 = dict()
dicti['proteomes'] = rep_names
dicti2['proteomes'] = rep_names
wyniki = []     
s = 0
current = ''
for i in organizmy2:
    #ilość hitow z danej bazy
    l = 0
    if len(i.all_hits) == 1:     
        wyniki.append(len(i.all_hits[list(i.all_hits.keys())[0]]))
        current = list(i.all_hits.keys())[0]
    else:
        wyniki.append(l)
    if len(wyniki) == 1:
        plt.figure(s)
        s+=1
        dicti2['blasthits'] = wyniki
        df = pd.DataFrame(dicti)
        df2 = pd.DataFrame(dicti2)
        sns.set(font_scale = 2.4)
        #sns.barplot(x = rep_names,y = 'whole', color = "skyblue",data = df)
        ax = sns.barplot(x=rep_names,y='blasthits',color="blue",data = df2)
        ax.set_ylabel("Number of hits from " + current,fontsize=24, alpha=0.8)
        ss = 0
        for p in ax.patches[0:int((len(ax.patches)))]:
            if ss<1:
                #percentage = str(round(round(1-(wyniki[ss]/whole_protein_count[ss]),4)*100,2))
                percentage = wyniki[ss]
                height = p.get_height()
                #ax.text(p.get_x()+p.get_width()/2.,height-500,str(percentage) + '\%',ha="center")
                ax.text(p.get_x()+p.get_width()/2.,height+10,str(percentage),ha="center")
                ss+=1
        print(wyniki)
        wyniki = []
        