# -*- coding: utf-8 -*-
#import numpy as np
from itertools import combinations
import time
          
class Protein:
    def __init__(self, ids, rep_names, codes):
        self.id = ids       #identyfikator białka (jego kod z pliku blast)
        self.hits = []      #lista kodów białek które blast zaliczył jako hit
        self.e_val = []     #lista watości e-value kolejnych hitów
        
        #słownik zawierający pary nazwa organizmu : lista kodów białek z tego organizmu które były hitami 
        self.taxa = {}
        for name in rep_names:
            self.taxa[name] = []
        for p2 in codes:
            for l in codes[p2]:
                if self.id[0:len(l)] == l:
                    self.organism = p2
                    del(self.taxa[p2])

    def __repr__(self):
        return self.id
        
    #dodawanie wartości do list białka
    def add_hit(self , hit):
        self.hits.append(hit)

    def add_eval(self , val):
        self.e_val.append(val)
        
    #funkcja która usuwa informacje o hitach z e-value gorszym niż wskazanye przez nas
    def eval_check(self, codes, cutoff):
        if len(self.e_val[-1].split("-")) > 1 and int(self.e_val[-1].split("-")[1]) >= cutoff: #ta linijka wznacza minimalną dopuszczalną wartosć ujemnej potęgi e-value jakie zaakceptujemy
                for p2 in codes:
                    if p2 != self.organism:
                        for l in codes[p2]:
                            if self.hits[-1][0:len(l)] == l:
                                #kolejne dwie linijki warunkują czy przechowujemy swoje id czy id hitów w słowniku zawierającym organizmy
                                #to znaczy albo wiemy , że to białko ma hit z tego organizmu , albo jakie białka z innych organizmów miały hit do tego białka 
                                self.taxa[p2].append(self.hits[-1])
        elif len(self.e_val[-1].split("."))>1 and int(self.e_val[-1].split(".")[1]) == 0 and int(self.e_val[-1].split(".")[0]) == 0:
            for p2 in codes:
                    for l in codes[p2]:
                        if self.hits[-1][0:len(l)] == l:
                            self.taxa[p2].append(self.hits[-1])
        else:   #jeżeli hit ma zbyt niskie e-value to jest usuwany
            del(self.hits[-1])
            del(self.e_val[-1])           

class Organism:
    """
    uogólnienie klasy Protein - przechowuję słownik zawierający:
    sumę analogicznych słowników ze wszystkich białek aktualnie badanego organizmu
    pozwala sprawdzić co właściwie dał nam blast - ile łącznie hitów było w pozostałych organizmach i jak się one mają do siebie
    można wyliczyć części wspólne
    pracuje na liście obiektów klasy Protein
    """
    def __init__(self , prot_list, names_faa, rep_names):
        self.id = list(set(rep_names) - set(prot_list[0].taxa.keys()))[0]
        self.all_hits = prot_list[0].taxa
        self.prots = []
        for prot in prot_list:
            if prot.organism == self.id:
                self.prots.append(prot)
                for each in prot.taxa:
                    self.all_hits[each] = self.all_hits[each] + prot.taxa[each]
        name = names_faa[0]
        del(names_faa[0])
        file_name = open(name , "r")
        current = set()
        for line in file_name:
            new = line.strip().split(">")
            if len(new) == 2:
                new = new[1].split(" ")
                current.add(new[0])
        self.all_hits[self.id] = current
        file_name.close()
    def __repr__(self):
        return str(self.id)

#1
def load_org_codes():
    results = open("codes" ,"r")
    #słownik "kodów" białek organizmów - napotykając białko w pliku wynikowym z blasta program dzięki temu słownikowi
    #potrafi przyporządkować białko do konkretnego organizmu
    codes = dict()
    current_org = None
    for line in results:
        if line.strip()[0] == ">":
            current_org = line.strip()[1:]
            codes[current_org] = []
        elif current_org in list(codes.keys()):
            codes[current_org].append(line.strip())
    results.close()
    print("Codes loaded")
    return codes

#2
def read_proteins(blast_names, names_faa, rep_names, codes, cutoff):                       
    #kod czytający pliki blastowe i tworzący na ich podstawie obiekty klasy Protein
    organisms = []
    for file in blast_names:
        blast_results = []
        file_name = open(file , "r")
        
        for k in file_name:
            if k[0:6] == "Query=":  #rozpoznaje wystąpienie nowego zapytania
                row = k.strip().split()
                blast_results.append(Protein(row[1], rep_names, codes))
            if k[0:2] == "> ":  #wykrywa pojawienie się nowego hitu
                row = k.strip().split()
                blast_results[-1].add_hit(row[1])
            if k[0:6] == " Score" and len(blast_results[-1].hits)>len(blast_results[-1].e_val): #spisuje wartość e-value ostatniego hitu
                row = k.strip().split()
                blast_results[-1].add_eval(row[7].strip(" ,"))
                blast_results[-1].eval_check(codes, cutoff) #sprawdza czy ostatnio dodany hit spełnia próg e-value narzucony przez funkcję eval_check
        
        org = Organism(blast_results, names_faa, rep_names) #tworzy obiekt klasy Organism na podstawie otrzymanej wyżej listy białek
        organisms.append(org)   #zapisuje obiekt Organism w liście
        print(str(file) + " loaded")
    print("All organisms loaded")         
    return organisms

#pracuje na liście obiektów typu Organism
#w zwróconej liście mamy kolejne listy przedstawiające części wspólne
def make_merge(organism_list, rep_names):
    #dla każdego organizmu w liście kod generuje listę wszystkich możliwych części wspólnych 2 organizmów
    #kolejne pętle funkcji są analogiczne tylko dla przecięć wszystkich możliwych kombinacji coraz większej liczby organizmów( 3 ,4 ,itd.).
    results = []
    k = list(range(len(rep_names)))
    if len(k) > 1:
        for x in range(2 ,len(rep_names)):
            hits = []
            for i in organism_list:      
                temp_k = list(combinations(k , x))
                temp = []
                item_list = list(i.all_hits.items())
                for i2 in temp_k:
                    temp_id = [item_list[i2[0]][0]]
                    temp_set = set(item_list[i2[0]][1])
                    for s2 in range(1,len(i2)):
                        temp_id.append(item_list[i2[s2]][0])
                        temp_set = temp_set & set(item_list[i2[s2]][1])
                    temp_id.append(temp_set)
                    temp_id.append(len(list(temp_set)))
                    temp.append(temp_id)
                hits.append(temp)
            results.append(hits)
        return results
    else:
        return("Not enough organisms provided.")

def make_absent_merge(organism_list, rep_names):
    #dla każdego organizmu w liście kod generuje listę wszystkich możliwych części wspólnych 2 organizmów
    #kolejne pętle funkcji są analogiczne tylko dla przecięć wszystkich możliwych kombinacji 3 ,4 i 5 organizmów.
    new_organism_list = organism_list
    for i in new_organism_list:
        for l in i.all_hits:
            if l != i.id:
                i.all_hits[l] = i.all_hits[i.id] - i.all_hits[l]
        i.all_hits[i.id] = i.all_hits[i.id] - i.all_hits[i.id]
    results = []
    k = list(range(len(rep_names)))
    if len(k) > 1:
        for x in range(1 ,len(rep_names)):
            hits = []
            for i in organism_list:      
                temp_k = list(combinations(k , x))
                temp = []
                item_list = list(i.all_hits.items())
                for i2 in temp_k:
                    temp_id = [item_list[i2[0]][0] ,item_list[i2[1]][0]]
                    temp_set = set(item_list[i2[0]][1]) & set(item_list[i2[1]][1])
                    temp_id.append(temp_set)
                    temp_id.append(len(list(temp_set)))
                    temp.append(temp_id)
                hits.append(temp)
            results.append(hits)
        return results
    else:
        return("Not enough organisms provided.")
    return results

#funkcja służy do zapisania w pliku białek z konkretnego organizmu mający hit do danego organizmu
def venn_prepare(organisms):
    for org in organisms:
        for org2 in org.all_hits:
            k = org.id+ "_" + org2 + ".list"
            save_file = open(k ,"w")
            for i in list(org.all_hits[org2]):
                save_file.write(i + "\n")
            save_file.close()

#wyszukuje w pliku proteomu białka unique_prots i zapisuje je do osobnego pliku 
def unique_save(unique_list,names_faa):
    for i in range(len(unique_list)):
        k = names_faa[i]+ ".unique"
        save_file = open(k ,"w")
        open_file = open(names_faa[i] ,"r")
        save = False
        p2 = 0
        for line in open_file:
            if line[0] == ">":
                if line.split(">")[1].split(" ")[0] in unique_list[i]:
                    save = True
                    p2+=1
                else:
                    save = False
            if save == True:
                save_file.write(line)
           
def cross_validation(organism_list):
    validated = {}
    for i in organism_list: #badamy kolejne organizmy
        print("Starting analysis for " + str(i.id))
        validated[i.id] = {}
        for another in organism_list:   #dobieramy drugi organizm , powstaje para dwoch badanych organizmow
            if another.id != i.id:
                validated[i.id][another.id] = set()
                tic = time.time()
                for protein in i.prots: #badamy białka organizmu pierwszego
                    main_id = protein.id
                    for hit in protein.hits:    #sprawdzamy czy hity są obustronne - bierzemy listę trafień z pierwszego organizmu do drugiego organizmu 
                        for protein_other in another.prots: #i badamy czy drugi organizm także posiada trafienie do danego białka z pierwszego organizmu
                            if protein_other.id == hit:
                                if main_id in protein_other.hits:
                                    validated[i.id][another.id].add(protein) #zapisujemy obustronny hit w liście
                                    
                toc = time.time() - tic
                print(str(toc) + " sec for computing pair " + str(i) + " -> " + str(another))
                print(len(validated[i.id][another.id]))
        print("Finished analysis for " + str(i.id))
    return validated

def protein_list_to_set(organisms):   
    for org in organisms:
        for prot in org.all_hits:
            org.all_hits[prot] = list(set(org.all_hits[prot]))

def make_validation(organisms,rep_names):
    print("Starting two-wayed blast hit validation")
    time_killer = cross_validation(organisms)
    print("Validation finished")
    for org1 in range(len(rep_names)):
        for org2 in range(len(rep_names)):
            if org1 != org2:
                first_name = rep_names[org1]
                second_name = rep_names[org2]
                new_list = []
                for protein in time_killer[rep_names[org1]][rep_names[org2]]:
                    new_list.append(protein.id)
                for x in organisms:
                    if x.id == first_name:
                        x.all_hits[second_name] = set(new_list)
                    
def save_validated_and_unique(organisms, names_faa, rep_names):
    unique_prots = []
    merged = make_merge(organisms,rep_names) 
    cross_results = open("cross_results.txt" ,'w')     
    for i in merged:
        l = 0
        for k in i:        
            for i2 in k:
                pp = i2[0:len(i2)-2]
                if rep_names[l] in pp:
                    pp.remove(rep_names[l])
                cross_results.write(str(rep_names[l]) + " " + str(pp) + " " + str(i2[-1]) + " \n")
            l+=1  
            if l == len(rep_names):
                l = 0
                    
    for i in range(len(rep_names)):
        unique = set(organisms[i].all_hits[rep_names[i]])
        for i2 in range(1 ,len(rep_names)):
            unique = set(unique) - set(organisms[i].all_hits[rep_names[i-i2]])
        unique_prots.append(unique)
    unique_save(unique_prots,names_faa)