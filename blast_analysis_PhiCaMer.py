from __future__ import unicode_literals
#-*- coding: utf-8 -*-
#import numpy as np
from copy import deepcopy
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import time

#nazwy "reprezentatywne" organizmów
rep_names = ["Prototheca","Chlorella variabilis","Auxenochlorella protothecoides","Helicosporidium sp.","Trichophyton rubrum","Candida albicans"]
base_names = ["PhiBase","MEROPS","CAZY"]

#słownik "kodów" białek organizmów - napotykając białko w pliku wynikowym z blasta program dzięki temu słownikowi
#potrafi przyporządkować białko do konkretnego organizmu

#wczytuje kody z pliku
results = open("codes","r")
codes = dict()
for line in results:
    if line.strip() in base_names:
        current_base = line.strip()
        codes[current_base] = []
    else:
        codes[current_base].append(line.strip())
results.close()

codes_org = dict()
codes_org["Chlorella variabilis"] = ["XP_005"]
codes_org['Auxenochlorella protothecoides'] = ['XP_011']
codes_org['Candida albicans'] = ['XP_436','XP_713','XP_019','XP_711','XP_710','XP_714','XP_717','XP_722','XP_709','XP_715','XP_721','XP_712','XP_443','XP_719','XP_441','XP_716','XP_718','XP_720','XP_723','XP_442']
codes_org['Helicosporidium sp.'] = ['KDD']
codes_org['Trichophyton rubrum'] = ['XP_003']
codes_org['Prototheca'] = ['genema', 'maker-', 'august']

#nazwy plików wynikowych blast z bazą phibase
names_phi = ["prototheca_proto.phibase","chlorella_proto.phibase","auxenochlorella_proto.phibase","helicosporidium_proto.phibase","trichophyton_proto.phibase","candida_proto.phibase"]

#nazwy plików wynikowych blast z bazą merops
names_merops = ["prototheca_proto.merops","chlorella_proto.merops","auxenochlorella_proto.merops","helicosporidium_proto.merops","trichophyton_proto.merops","candida_proto.merops"]

#nazwy plików wynikowych blast z bazą cazy
names_cazy = ["prototheca_proto.cazy","chlorella_proto.cazy","auxenochlorella_proto.cazy","helicosporidium_proto.cazy","trichophyton_proto.cazy","candida_proto.cazy"]

#łączna ilość linii w plikach: 6 064 393
names = names_phi+names_merops+names_cazy

#Rozmiar proteomów organizmów (w moim przypadku z ncbi genome) w kolejności: Chlorella, Auxenochlorella, Helicosporidum,Trichophyton, Candida
whole_protein_count = [733,2680, 1401, 1125, 4044, 1853]

#Ilość sekwencji w bazach w kolejności: phibase,Merops,Cazy
whole_base_count = [4003,4967,921174]

#nazwy plikow baz danych w formacie fasta
db_names = ["phi_accessions.fasta","merops_scan.fasta","CAZyDB.fasta"]

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
        #self.identities = []    #lista ułamków obrazujących jak podobne jest białko do hitu pod danym indeksem
        self.e_val = []     #lista watości e-value kolejnych hitów
        
        #słownik zawierający pary nazwa rganizmu : lista kodów białek z tego organizmu które były hitami 
        self.taxa = {"PhiBase": [], "MEROPS": [], "CAZY": []}
        for p2 in codes_org:
            for l in codes_org[p2]:
                if self.id[0:len(l)] == l:
                    self.organism = p2

    def __repr__(self):
        return self.id
        
    #dodawanie wartości do list białka
    def add_hit(self,hit):
        self.hits.append(hit)

    #def add_ident(self,ident):
    #   self.identities.append(ident)

    def add_eval(self,val):
        self.e_val.append(val)
        
    #funkcja która usuwa informacje o hitach z e-value gorszym niż wskazane przez nas
    def eval_check(self):
        if len(self.e_val[-1].split("-")) > 1 and int(self.e_val[-1].split("-")[1]) >= 10: #ta linijka wyznacza minimalną dopuszczalną wartosć ujemnej potęgi e-value jakie zaakceptujemy
            for p2 in codes:
                for l in codes[p2]:
                    if self.hits[-1][0:len(l)] == l:
                        #kolejne dwie linijki warunkują czy przechowujemy swoje id czy id hitów w słowniku zawierającym organizmy
                        #to znaczy albo wiemy, że to białko ma hit z tego organizmu, albo jakie białka z innych organizmów miały hit do tego białka 
                        #self.taxa[p2].append(self.hits[-1])
                        self.taxa[p2].append(self.id)

        elif len(self.e_val[-1].split("."))>1 and int(self.e_val[-1].split(".")[1]) == 0 and int(self.e_val[-1].split(".")[0]) == 0:
            for p2 in codes:
                    for l in codes[p2]:
                        if self.hits[-1][0:len(l)] == l:
                            #self.taxa[p2].append(self.hits[-1])
                            self.taxa[p2].append(self.id)
        else:   #jeżeli hit ma zbyt niskie e-value to jest usuwany
            del(self.hits[-1])
            #del(self.identities[-1])
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
    tic = time.time()
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
        p+=1
        if p%100000 == 0:
            toc = time.time() - tic
            print(p,toc)
            tic = time.time()
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
rep_names = ["\it{Prototheca}","\it{Chlorella} \n \it{variabilis}","\it{Auxenochlorella} \n \it{protothecoides}","\it{Helicosporidium sp.}","\it{Trichophyton} \n \it{rubrum}","\it{Candida} \n \it{albicans}"]
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
        #sns.barplot(x = rep_names,y = 'whole', color = "skyblue",data = df)
        ax = sns.barplot(x=rep_names,y='blasthits',color="blue",data = df2)
        ax.set_ylabel("Number of organism proteins with hits from " + current,fontsize=24, alpha=0.8)
        ss = 0
        for p in ax.patches[0:int((len(ax.patches)))]:
            if ss<6:
                height = p.get_height()
                ax.text(p.get_x()+p.get_width()/2.,height+3,str(temp[ss]),ha="center")
                ss+=1      
        s+=1
        wyniki = wyniki + temp
        temp = []
        plt.gcf().subplots_adjust(bottom=0.1)
        plt.show()
        #fig = ax.get_figure()
        #fig.savefig(str(s) + "_proto_uniq.png")
        
            
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