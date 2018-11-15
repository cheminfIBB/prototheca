from __future__ import unicode_literals
from copy import deepcopy
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

class Protein:
    def __init__(self, ids, rep_names, base_names, base_codes, codes_org):
        self.id = ids       #identyfikator białka (jego kod z pliku blast)
        self.hits = []      #lista kodów białek które blast zaliczył jako hit
        self.e_val = []     #lista watości e-value kolejnych hitów
        #słownik zawierający pary nazwa organizmu : lista kodów białek z tego organizmu które były hitami 
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
    def add_hit(self, hit):
        self.hits.append(hit)

    def add_eval(self, val):
        self.e_val.append(val)
        
    #funkcja która usuwa informacje o hitach z e-value gorszym niż wskazanye przez nas
    def eval_check(self, base_codes, cutoff):
        if len(self.e_val[-1].split("-")) > 1 and int(self.e_val[-1].split("-")[1]) >= int(cutoff): #ta linijka wyznacza minimalną dopuszczalną wartosć ujemnej potęgi e-value jakie zaakceptujemy
            for p2 in base_codes:
                for l in base_codes[p2]:
                    if self.hits[-1][0:len(l)] == l:
                        self.taxa[p2].append(self.id)

        elif len(self.e_val[-1].split("."))>1 and int(self.e_val[-1].split(".")[1]) == 0 and int(self.e_val[-1].split(".")[0]) == 0:
            for p2 in base_codes:
                    for l in base_codes[p2]:
                        if self.hits[-1][0:len(l)] == l:
                            self.taxa[p2].append(self.id)
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
    def __init__(self , prot_list):
        self.id = prot_list[0].organism
        self.all_hits = prot_list[0].taxa
        for prot in prot_list:
            for each in prot.taxa:
                self.all_hits[each] = self.all_hits[each] + prot.taxa[each]
    def __repr__(self):
        return str(self.id)

#1
def load_org_and_base_codes(org_fasta, base_names, rep_names):
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
    base_codes = dict()
    for line in results:
        if line.strip().strip(">") in base_names:
            current_base = line.strip().strip(">")
            base_codes[current_base] = []
        else:
            base_codes[current_base].append(line.strip().strip(">"))
    results.close()
    return base_codes, codes_org
#2
def read_proteins(names_faa, rep_names, base_names, base_codes, codes_org, cutoff):                       
    #kod czytający pliki blastowe i tworzący na ich podstawie obiekty klasy Protein
    organisms = []
    names = []
    for base_name in base_names:
        for org_name in names_faa:
            file_name = org_name.split(".")[0]
            file_name = file_name + "." + str(base_name).lower()
            names.append(file_name)
            
    for i in range(0,len(names)):
        blast_results = []
        name = names[i] #aktualnie czytany plik
        rep_name = rep_names[i%len(rep_names)] #nazwa organizmu dla którego czytamy wyniki blasta
        file_name = open(name, "r")
        for k in file_name:
            if k[0:6] == "Query=":  #rozpoznaje wystąpienie nowego zapytania
                row = k.strip().split()
                blast_results.append(Protein(row[1], rep_names, base_names, base_codes, codes_org))
            if k[0:2] == "> ":  #wykrywa pojawienie się nowego hitu
                row = k.strip().split()
                blast_results[-1].add_hit(row[1])
            if k[0:6] == " Score" and len(blast_results[-1].hits)>len(blast_results[-1].e_val): #spisuje wartość e-value ostatniego hitu
                row = k.strip().split()
                blast_results[-1].add_eval(row[7].strip(","))
                blast_results[-1].eval_check(base_codes, cutoff) #sprawdza czy ostatnio dodany hit spełnia próg e-value narzucony przez funkcję eval_check
        file_name.close()
        org = Organism(blast_results) #tworzy obiekt klasy Organism na podstawie otrzymanej wyżej listy białek
        organisms.append(org)   #zapisuje obiekt Organism w liście globalnej 
        print(name+" loaded")
    return organisms
    
def protein_list_to_set(organisms):   
    for org in organisms:
        for prot in org.all_hits:
            org.all_hits[prot] = list(set(org.all_hits[prot]))

def make_plots(organisms, rep_names, tex):
    
    for i in organisms:
        for k, v in list(i.all_hits.items()):
             if len(v) == 0:
                 del i.all_hits[k]
                 
    if tex is True:
        plt.rc('text', usetex=True)
    else:
        plt.rc('text', usetex=False)
        
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
        
    sns.set_style("whitegrid")
    dicti2 = dict()
    dicti2['proteomes'] = rep_names
    
    results = []     
    temp = []
    s = 0
    current = ''
    for i in organisms:
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
            plt.figure(figsize=(20,10))
            dicti2['blasthits'] = temp
            df2 = pd.DataFrame(dicti2)
            sns.set(font_scale = 2.4)
            plot_names = deepcopy(rep_names)
            if tex is True:
                for i in range(len(plot_names)):
                    plot_names[i] = r'\textit{' + plot_names[i] + '}'
            ax = sns.barplot(x=plot_names,y='blasthits',color="blue",data = df2)
            ax.set_ylabel("Number of organism proteins with hits from " + current,fontsize=24, alpha=0.8)
            ss = 0
            for p in ax.patches:
                ax.text(p.get_x()+p.get_width()/2.,p.get_height()+1,str(temp[ss]),ha="center")
                ss+=1      
            s+=1
            results = results + temp
            temp = []
            plt.gcf().subplots_adjust(bottom=0.1)
            plt.savefig(str(current)+"_plot",bbox_inches="tight",pad_inches=0.5)
            plt.close()
    return results

def save_results(results, rep_names, base_names):   
    results_file = open('base_results.txt','w')
    l=0
    s=0
    for k in results:
        if l == 0:
            results_file.write("> " + base_names[s] + " \n")
            s+=1
        results_file.write(rep_names[l] + ':: ' + str(k) + ' \n')
        l+=1
        if l == len(rep_names):
            l = 0
    results_file.close()