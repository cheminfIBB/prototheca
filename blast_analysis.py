# -*- coding: utf-8 -*-
#import numpy as np
from itertools import combinations
import pandas

#słownik "kodów" białek organizmów - napotykając białko w pliku wynikowym z blasta program dzięki temu słownikowi
#potrafi przyporządkować białko do konkretnego organizmu
codes = dict()
codes["Chlorella variabilis"] = ["XP_005"]
codes['Auxenochlorella protothecoides'] = ['XP_011']
codes['Candida albicans'] = ['XP_436','XP_713','XP_019','XP_711','XP_710','XP_714','XP_717','XP_722','XP_709','XP_715','XP_721','XP_712','XP_443','XP_719','XP_441','XP_716','XP_718','XP_720','XP_723','XP_442']
codes['Helicosporidium sp.'] = ['KDD']
codes['Trichophyton rubrum'] = ['XP_003']

#nazwy plików fasta zawierających proteomy 
names_faa = ["chlorella_variabilis.faa","trichophyton_rubrum.faa","auxenochlorella_protothecoides.faa","helicosporidium_sp.faa","candida_albicans.faa"]

#Rozmiar proteomów organizmów (w moim przypadku z ncbi genome) w kolejności: Chlorella, Auxenochlorella, Helicosporidum,Trichophyton, Candida
whole_protein_count = [9780, 7014, 6033, 8706, 6030]

#nazwy plików wynikowych blasta, też muszą być w tej samej kolejności
names = ["no_variabilis.out","no_rubrum.out","no_prototheco.out","no_helico.out","no_albicans.out"]

#nazwy "reprezentatywne" organizmów
rep_names = ["Chlorella variabilis","Trichophyton rubrum","Auxenochlorella protothecoides","Helicosporidium sp.","Candida albicans"]
rep_names2 = ["Chlorella variabilis","Auxenochlorella protothecoides","Helicosporidium sp."]

#lista zawierająca obiekty klasy organizm - czyli obrobione wyniki z blasta
organizmy = []

#lista zawiera dla każdego białka z kolejnych organizmów obiekty klasy Matrix odpowiadające tym białkom
organizmy_macierze = []

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
        self.taxa = {"Chlorella variabilis": [], "Trichophyton rubrum": [], "Auxenochlorella protothecoides": [],
                     "Helicosporidium sp.": [], "Candida albicans": []}
        for p2 in codes:
            for l in codes[p2]:
                if self.id[0:len(l)] == l:
                    self.organism = p2
                    del(self.taxa[p2])

    def __repr__(self):
        return self.id
        
    #dodawanie wartości do list białka
    def add_hit(self,hit):
        self.hits.append(hit)

    def add_ident(self,ident):
        self.identities.append(ident)

    def add_eval(self,val):
        self.e_val.append(val)
        
    #funkcja która usuwa informacje o hitach z e-value gorszym niż wskazanye przez nas
    def eval_check(self):
        if len(self.e_val[-1].split("-")) > 1 and int(self.e_val[-1].split("-")[1]) >= 10: #ta linijka wznacza minimalną dopuszczalną wartosć ujemnej potęgi e-value jakie zaakceptujemy
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
            del(self.identities[-1])
            del(self.e_val[-1])           

#funkcja potrafi przedstawić w formie DataFrame'a swoje atrybuty hits/identitys/e_val
#potrafi także zapisywać takiego dataframe do wskazanego przez nas pliku w formie CSV, wartości poprzedzielane są spacjami
#pracuje na obiekcie klasy Protein
class Matrix:
    def __init__(self,prot):
        self.id = prot.id
        for x,v in list(prot.taxa.items()):
            if v == 0:
                del(prot.taxa[x])
        self.taxa = prot.taxa
        data_dict = {"Identity": prot.identities, "E-Value": prot.e_val}
        self.matrix = pandas.DataFrame(data_dict,index = prot.hits)
        
    def __repr__(self):
        return str(self.id)
        
    def write_results(self,file):
        file.write("> " + self.id +" \n")
        for i2 in list(self.matrix.index.values):
            temp = i2 + " "
            for l in self.matrix.loc[i2]:
                temp += l+ " "
            temp += "\n"
            file.write(temp)
            
#uogólnienie klasy Protein - przechowuję słownik zawierający:
#sumę analogicznych słowników ze wszystkich białek aktualnie badanego organizmu
#pozwala sprawdzić co właściwie dał nam blast - ile łącznie hitów było w pozostałych organizmach i jak się one mają do siebie
#można wyliczyć części wspólne
#pracuje na liście obiektów klasy Protein
class Organism:
    def __init__(self,prot_list):   
        self.id = list(set(rep_names) - set(prot_list[0].taxa.keys()))[0]
        self.all_hits = prot_list[0].taxa
        for prot in prot_list:
            if prot.organism == self.id:
                for each in prot.taxa:
                    self.all_hits[each] = self.all_hits[each] + prot.taxa[each]
        name = names_faa[0]
        del(names_faa[0])
        file_name = open(name, "r")
        current = set()
        for line in file_name:
            new = line.strip().split(">")
            if len(new) == 2:
                new = new[1].split(" ")
                current.add(new[0])
        self.all_hits[self.id] = current
        file_name.close()
    def __repr__(self):
        print(self.id)
        return ""
            
#pracuje na liście obiektów typu Organism
#ta funkcja działa tylko i wyłącznie dla 5 organizmów - zwraca listę:
#zwrócona lista zawiera w sobie cztery listy - dla 2,3,4,5 części współnych 
#każda z czterech list ma w sobie pięc list - po jednej dla każdego organizmu w porządku takim samym jak w liście jaką mu podaliśmy
#w liście organizmów mamy kolejne listy przedstawiające części wspólne 
def make_merge(organism_list):
    two_hits = []
    #dla każdego organizmu w liście kod generuje listę wszystkich możliwych części wspólnych 2 organizmów
    #kolejne pętle funkcji są analogiczne tylko dla przecięć wszystkich możliwych kombinacji 3,4 i 5 organizmów.
    for i in organism_list:
        k = [0,1,2,3,4]
        temp_k = list(combinations(k, 2))
        temp = []
        item_list = list(i.all_hits.items())
        for i2 in temp_k:
            temp_id = [item_list[i2[0]][0],item_list[i2[1]][0]]
            temp_set = item_list[i2[0]][1] & item_list[i2[1]][1]
            temp_id.append(temp_set)
            temp_id.append(len(list(temp_set)))
            temp.append(temp_id)
        two_hits.append(temp)
        
    three_hits = []
    
    for i in organism_list:
        k = [0,1,2,3,4]
        temp_k = list(combinations(k, 3))
        temp = []
        item_list = list(i.all_hits.items())
        for i2 in temp_k:
            temp_id = [item_list[i2[0]][0],item_list[i2[1]][0],item_list[i2[2]][0]]
            temp_set = item_list[i2[0]][1] & item_list[i2[1]][1] & item_list[i2[2]][1]
            temp_id.append(temp_set)
            temp_id.append(len(list(temp_set)))
            temp.append(temp_id)
        three_hits.append(temp)

    four_hits = []
    
    for i in organism_list:
        k = [0,1,2,3,4]
        temp_k = list(combinations(k, 4))
        temp = []
        item_list = list(i.all_hits.items())
        for i2 in temp_k:
            temp_id = [item_list[i2[0]][0],item_list[i2[1]][0],item_list[i2[2]][0],item_list[i2[3]][0]]
            temp_set = item_list[i2[0]][1] & item_list[i2[1]][1] & item_list[i2[2]][1] & item_list[i2[3]][1]
            temp_id.append(temp_set)
            temp_id.append(len(list(temp_set)))
            temp.append(temp_id)
        four_hits.append(temp)
    
    five_hits = []
    for i in organism_list:
        k = [0,1,2,3,4]
        temp_k = list(combinations(k, 5))
        temp = []
        item_list = list(i.all_hits.items())
        for i2 in temp_k:
            temp_id = [item_list[i2[0]][0],item_list[i2[1]][0],item_list[i2[2]][0],item_list[i2[3]][0],item_list[i2[4]][0]]
            temp_set = item_list[i2[0]][1] & item_list[i2[1]][1] & item_list[i2[2]][1] & item_list[i2[3]][1] & item_list[i2[4]][1]
            temp_id.append(temp_set)
            temp_id.append(len(list(temp_set)))
            temp.append(temp_id)
        five_hits.append(temp)
    results = [two_hits,three_hits,four_hits,five_hits]
    return results

#funkcja służy do 
def venn_prepare(organizmy):
    for org in organizmy:
        for org2 in org.all_hits:
            k = org.id+ "_" + org2 + ".list"
            save_file = open(k,"w")
            for i in list(org.all_hits[org2]):
                save_file.write(i + "\n")
            save_file.close()
    
#kod czytający pliki blastowe i tworzący na ich podstawie obiekty klasy Protein           
for i in range(0,len(names)):
    blast_results = []
    name = names[i] #aktualnie czytany plik
    rep_name = rep_names[i] #nazwa organizmu dla którego czytamy wyniki blasta
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
    org = Organism(blast_results) #tworzy obiekt klasy Organism na podstawie otrzymanej wyżej listy białek
    organizmy.append(org)   #zapisuje obiekt Organism w liście globalnej
    macierze = []
    for i in blast_results: #wytwarza dla każdego białka reprezentację jego wyników w formie obiektu klasy Matrix i zapisuje taką listę
        k = Matrix(i)
        macierze.append(k)
    organizmy_macierze.append(macierze)

for i in organizmy:
    for i2 in i.all_hits:
        i.all_hits[i2] = set(i.all_hits[i2])
        
merged = make_merge(organizmy)            
for i in merged:
    for k in i:
        print("new organism")
        for i2 in k:
            print(i2[0:len(i2)-2],i2[-1])