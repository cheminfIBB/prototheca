from __future__ import unicode_literals
#-*- coding: utf-8 -*-
#import numpy as np
from copy import deepcopy
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import time

#nazwy "reprezentatywne" organizmow
rep_names = ["Prototheca","Chlorella variabilis","Auxenochlorella protothecoides","Helicosporidium sp.","Trichophyton rubrum","Candida albicans"]
base_names = ["PhiBase","MEROPS","CAZY"]

#slownik "kodow" bialek organizmow - napotykajac bialko w pliku wynikowym z blasta program dzieki temu slownikowi
#potrafi przyporzadkowac bialko do konkretnego organizmu

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

#nazwy plikow wynikowych blast z baza phibase
names_phi = ["prototheca.phibase","chlorella_variabilis.phibase","auxenochlorella_protothecoides.phibase","helicosporidium_sp.phibase","trichophyton_rubrum.phibase","candida_albicans.phibase"]

#nazwy plikow wynikowych blast z baza merops
names_merops = ["prototheca.merops","chlorella_variabilis.merops","auxenochlorella_protothecoides.merops","helicosporidium_sp.merops","trichophyton_rubrum.merops","candida_albicans.merops"]

#nazwy plikow wynikowych blast z baza cazy
names_cazy = ["prototheca.cazy","chlorella_variabilis.cazy","auxenochlorella_protothecoides.cazy","helicosporidium_sp.cazy","trichophyton_rubrum.cazy","candida_albicans.cazy"]

#laczna ilosc linii w plikach: 36 830 685
names = names_phi+names_merops+names_cazy

#Rozmiar proteomow organizmow (w moim przypadku z ncbi genome) w kolejnosci: Chlorella, Auxenochlorella, Helicosporidum,Trichophyton, Candida
whole_protein_count = [5930, 9780, 7014, 6033, 8706, 6030]

#Ilosc sekwencji w bazach w kolejnosci: phibase,Merops,Cazy
whole_base_count = [4003,4967,921174]

#nazwy plikow baz danych w formacie fasta
db_names = ["phi_accessions.fasta","merops_scan.fasta","CAZyDB.fasta"]

#lista zawierajaca obiekty klasy organizm - czyli obrobione wyniki z blasta
organizmy = []

#lista zawierajaca zbior id bialek unikalnych dla kazdego organizmu
unikalne = []

#w pozniejszym etapie zawiera liste wygenerowana funkcja make_merge na podstawie listy organizmy
#zawiera wszystkie mozliwe kombinacje czesci wspolnych zbiorow kodow z kazdego organizmu
merged = []  
    
class Protein:
    def __init__(self,ids):
        self.id = ids       #identyfikator bialka (jego kod z pliku blast)
        self.hits = []      #lista kodow bialek ktore blast zaliczyl jako hit
        #self.identities = []    #lista ulamkow obrazujacych jak podobne jest bialko do hitu pod danym indeksem
        self.e_val = []     #lista watosci e-value kolejnych hitow
        
        #slownik zawierajacy pary nazwa rganizmu : lista kodow bialek z tego organizmu ktore byly hitami 
        self.taxa = {"PhiBase": [], "MEROPS": [], "CAZY": []}
        for p2 in codes_org:
            for l in codes_org[p2]:
                if self.id[0:len(l)] == l:
                    self.organism = p2

    def __repr__(self):
        return self.id
        
    #dodawanie wartosci do list bialka
    def add_hit(self,hit):
        self.hits.append(hit)

    #def add_ident(self,ident):
    #   self.identities.append(ident)

    def add_eval(self,val):
        self.e_val.append(val)
        
    #funkcja ktora usuwa informacje o hitach z e-value gorszym niz wskazane przez nas
    def eval_check(self):
        if len(self.e_val[-1].split("-")) > 1 and int(self.e_val[-1].split("-")[1]) >= 10: #ta linijka wyznacza minimalna dopuszczalna wartosc ujemnej potegi e-value jakie zaakceptujemy
            for p2 in codes:
                for l in codes[p2]:
                    if self.hits[-1][0:len(l)] == l:
                        #kolejne dwie linijki warunkuja czy przechowujemy swoje id czy id hitow w slowniku zawierajacym organizmy
                        #to znaczy albo wiemy, ze to bialko ma hit z tego organizmu, albo jakie bialka z innych organizmow mialy hit do tego bialka 
                        #self.taxa[p2].append(self.hits[-1])
                        self.taxa[p2].append(self.id)

        elif len(self.e_val[-1].split("."))>1 and int(self.e_val[-1].split(".")[1]) == 0 and int(self.e_val[-1].split(".")[0]) == 0:
            for p2 in codes:
                    for l in codes[p2]:
                        if self.hits[-1][0:len(l)] == l:
                            #self.taxa[p2].append(self.hits[-1])
                            self.taxa[p2].append(self.id)
        else:   #jezeli hit ma zbyt niskie e-value to jest usuwany
            del(self.hits[-1])
            #del(self.identities[-1])
            del(self.e_val[-1])           
         
#uogolnienie klasy Protein - przechowuje slownik zawierajacy:
#sume analogicznych slownikow ze wszystkich bialek aktualnie badanego organizmu
#pozwala sprawdzic co wlasciwie dal nam blast - ile lacznie hitow bylo w pozostalych organizmach i jak sie one maja do siebie
#mozna wyliczyc czesci wspolne
#pracuje na liscie obiektow klasy Protein
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

#kod czytajacy pliki blastowe i tworzacy na ich podstawie obiekty klasy Protein          
for i in range(0,len(names)):
    blast_results = []
    name = names[i] #aktualnie czytany plik
    print(name)
    rep_name = rep_names[i%len(rep_names)] #nazwa organizmu dla ktorego czytamy wyniki blasta
    file_name = open(name, "r")
    p=0
    tic = time.time()
    for k in file_name:
        if k[0:6] == "Query=":  #rozpoznaje wystapienie nowego zapytania
            row = k.strip().split()
            blast_results.append(Protein(row[1]))
        if k[0:2] == "> ":  #wykrywa pojawienie sie nowego hitu
            row = k.strip().split()
            blast_results[-1].add_hit(row[1])
        if k[0:6] == " Score" and len(blast_results[-1].hits)>len(blast_results[-1].e_val): #spisuje wartosc e-value ostatniego hitu
            row = k.strip().split()
            blast_results[-1].add_eval(row[7].strip(","))
            blast_results[-1].eval_check() #sprawdza czy ostatnio dodany hit spelnia prog e-value narzucony przez funkcje eval_check
            
        #if k[0:6] == " Ident"  and len(blast_results[-1].hits)>len(blast_results[-1].identities): #spisuje wartosc identycznosci pomiedzy zapytaniem a ostatnim hitem
            #row = k.strip().split()
            #blast_results[-1].add_ident(row[2])
        p+=1
        if p%100000 == 0:
            toc = time.time() - tic
            print(p,toc)
            tic = time.time()
    file_name.close()
    org = Organism(blast_results) #tworzy obiekt klasy Organism na podstawie otrzymanej wyzej listy bialek
    organizmy.append(org)   #zapisuje obiekt Organism w liscie globalnej

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
#dicti = dict()
#dicti['whole'] = whole_protein_count
dicti2 = dict()
#dicti['proteomes'] = rep_names
dicti2['proteomes'] = rep_names
wyniki = []    
temp = [] 
s = 0
current = ''
phi_uni = [0,36,2,0,345,137]
mer_uni = [0,65,1,3,62,13]
caz_uni = [0,148,3,9,140,54]
r = 0
for i in organizmy2:

    #ilosc hitow z danej bazy
    l = 0
    if len(i.all_hits) == 1:     
        temp.append(len(i.all_hits[list(i.all_hits.keys())[0]]))
        current = list(i.all_hits.keys())[0]
    else:
        temp.append(l)
    if len(temp) == len(rep_names):
        dicti = dict()
        dicti['proteomes'] = rep_names
        sns.set()
        plt.figure(s)
        s+=1
        if r ==0:
            dicti['blasthits2'] =  phi_uni
            searcher = temp + phi_uni
        if r ==1:
            dicti['blasthits2'] = mer_uni  
            searcher = temp + mer_uni
        else:
            dicti['blasthits2'] = caz_uni
            searcher = temp + caz_uni
        r+=1
        dicti2['blasthits'] = temp
        df = pd.DataFrame(dicti)
        df2 = pd.DataFrame(dicti2)
        sns.set(font_scale = 2.4)
        #sns.barplot(x = rep_names,y = 'whole', color = "skyblue",data = df)
        ap = sns.barplot(x=rep_names,y='blasthits',color="skyblue",data = df2)
        ax = sns.barplot(x=rep_names,y='blasthits2',color="blue",data = df)
        ax.set_ylabel("Number of organism proteins with hits from " + current,fontsize=24, alpha=0.8)
        ss = 0
        for p in ax.patches:
            height = p.get_height()
            ax.text(p.get_x()+p.get_width()/2.,p.get_height()+1,str(searcher[ss]),ha="center")
            ss+=1
        wyniki = wyniki + temp
        temp = []
        plt.gcf().subplots_adjust(bottom=0.1)
        plt.show()
        #fig = ax.get_figure()
        #fig.savefig(str(s) + "_proto_uniq.png")
        
'''           
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
'''
dicti = dict()
dicti['proteomes'] = rep_names
dicti2 = dict()
dicti2['proteomes'] = rep_names

wyniki = []    
temp = [] 

current = ''
phi = [977,1653,1205,962,2087,1552]
phi_uni = [0,36,2,0,345,137]
mer = [211,361,242,222,266,188]
mer_uni = [0,65,1,3,62,13]
caz = [648,1378,826,641,1017,744]
caz_uni = [0,148,3,9,140,54]
for i in range(3):
        sns.set()
        plt.figure(i)
        if i == 0:
            dicti['blasthits'] =  phi_uni
            dicti2['blasthits'] = phi
            searcher = phi + phi_uni
            current = 'PhiBase'
        elif i == 1:
            dicti['blasthits'] = mer_uni
            dicti2['blasthits'] = mer
            searcher = mer + mer_uni
            current = 'MEROPS'
        else:
            dicti['blasthits'] = caz_uni
            dicti2['blasthits'] = caz
            searcher = caz + caz_uni
            current = 'CAZY'
        df = pd.DataFrame(dicti)
        df2 = pd.DataFrame(dicti2)
        sns.set(font_scale = 2.4)
        ap = sns.barplot(x=rep_names,y='blasthits',color="skyblue",data = df2)
        ax = sns.barplot(x=rep_names,y='blasthits',color="blue",data = df)
        ax.set_ylabel("Number of organism proteins with hits from " + current,fontsize=24, alpha=0.8)
        print(len(ax.patches))
        ss = 0
        for p in ax.patches:
            height = p.get_height()
            ax.text(p.get_x()+p.get_width()/2.,p.get_height()+1,str(searcher[ss]),ha="center")
            ss+=1
        plt.gcf().subplots_adjust(bottom=0.1)
        plt.show()