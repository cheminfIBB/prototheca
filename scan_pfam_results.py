from __future__ import unicode_literals
import pandas as pd
#import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import copy

#nazwy plików do wczytania
names = ["prototheca.fasta.tsv", "chlorella_variabilis.faa.tsv","auxenochlorella_protothecoides.faa.tsv","helicosporidium_sp.faa.tsv",
"trichophyton_rubrum.faa.tsv","candida_albicans.faa.tsv"]
#nazwy organizmów
rep_names = ["Prototheca", "Chlorella variabilis","Auxenochlorella protothecoides","Helicosporidium sp","Trichophyton rubrum","Candida albicans"]
#nazwy organizmów kursywą
rep_names2 = ["\it{Prototheca}","\it{Chlorella} \n \it{variabilis}","\it{Auxenochlorella} \n \it{protothecoides}","\it{Helicosporidium sp.}","\it{Trichophyton} \n \it{rubrum}","\it{Candida} \n \it{albicans}"]
#nazwy plików z proteomami organizmów
names_faa = ["prototheca.fasta.unique","chlorella_variabilis.faa.unique","auxenochlorella_protothecoides.faa.unique","helicosporidium_sp.faa.unique","trichophyton_rubrum.faa.unique","candida_albicans.faa.unique"]

#["prototheca.fasta.unique","chlorella_variabilis.faa.unique","auxenochlorella_protothecoides.faa.unique","helicosporidium_sp.faa.unique","trichophyton_rubrum.faa.unique","candida_albicans.faa.unique"]
#["prototheca.fasta","chlorella_variabilis.faa","auxenochlorella_protothecoides.faa","helicosporidium_sp.faa","trichophyton_rubrum.faa","candida_albicans.faa"]
#names_tsv = ["prototheca.fasta.tsv_pfam.csv","chlorella_variabilis.faa.tsv_pfam.csv","auxenochlorella_protothecoides.faa.tsv_pfam.csv","helicosporidium_sp.faa.tsv_pfam.csv","trichophyton_rubrum.faa.tsv_pfam.csv","candida_albicans.faa.tsv_pfam.csv"]

# przeprowadzone analizy
analysis = ['TMHMM','SignalP_EUK','Coils','Phobius','MobiDBLite','PIRSF', 'TIGRFAM', 'TMHMM', 'PANTHER', 'SMART', 'SUPERFAMILY', 'Phobius', 'Pfam', 'Gene3D',
            'ProSitePatterns', 'MobiDBLite', 'PRINTS', 'ProSiteProfiles', 'SFLD', 'SignalP_EUK', 'Coils', 'Hamap',
            'ProDom', 'CDD']

# legenda do wyników z interproscan
headers = ["Sequence ID", "Sequence MD5 digest", "Sequence Length",
                   "Analysis", "Signature Accession", "Signature Description",
                   "Start Location", "Stop Location", "Score (e-value)", "Status",
                   "Run Date", "InterPro annotations description", "GO annotations",
                   "Pathways annotations"]

whole_protein_count = []
analysis_protein_count = []
unique = []
pfams = []
dictus = {}
for i in analysis:
    dictus[i] = set()
    
maine_dictus = {}

for i in range(len(rep_names)):
    maine_dictus[rep_names[i]] = copy.deepcopy(dictus)
    name = names_faa[i]
    file_name = open(name, "r")
    protein_names = []
    for line in file_name:      #robi listę białek których szukamy (mogą być to białka unikalne, może być cały proteom) na podstawie pliku fasta
        new = line.strip().split(">")
        if len(new) == 2:
            protein_names.append(new[1].split(" ")[0])
    file_name.close()
    whole_protein_count.append(len(protein_names))
    
    k = str(names[i])
    p = pd.read_csv(k, sep='\t', header = None, names=headers) #wczytuje wyniki interproscan i sprawdza linia po linii w poszukiwaniu
    for i2 in range(1, p.shape[0]):                            #interesujących nas białek i analiz
            row = p.iloc[[i2]].values[0]
            if row[0] in protein_names and row[3] in analysis:
                maine_dictus[rep_names[i]][row[3]].add(row[0])
    print('Loading ' + rep_names[i] + ' done')

plt.rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
rr = 0
for i in analysis:    # pętla robiąca wykres słupkowy dla każdej analizy, każdy organizm 
    plt.figure(rr)    # na wykresie zawiera liczbę wszystkich białek oraz liczbę białek z jego proteomu mające hit w czasie analizy 
    rr += 1
    #whole_protein_count = [733,2512, 1065, 1051, 4039, 1845]
    analysis_protein_count = []
    for name in rep_names:
        analysis_protein_count.append(len(maine_dictus[name][i]))
    sns.set()
    sns.set_style("whitegrid")
    dicti = dict()
    dicti["proteomes"] = rep_names2
    dicti["whole"] = whole_protein_count
    dicti2 = dict()
    dicti2["proteomes"] = rep_names2
    dicti2[i] = analysis_protein_count
    df = pd.DataFrame(dicti)
    df2 = pd.DataFrame(dicti2)
    #sns.barplot(data = df)
    sns.set(font_scale = 2.4)
    ap = sns.barplot(x = rep_names2,y = 'whole', color = "skyblue",data = df)
    ax = sns.barplot(x=rep_names2,y=i,color="blue",data=df2)
    ax.set_ylabel("Number of proteins with hit from " + i,fontsize=24, alpha=0.8)
    searcher = whole_protein_count + analysis_protein_count
    s = 0
    for p in ax.patches:
        height = p.get_height()
        ax.text(p.get_x()+p.get_width()/2.,p.get_height()+1,str(searcher[s]),ha="center")
        #percentage = str(round(round(1-(analysis_protein_count[s]/whole_protein_count[s]),4)*100,2))  <- czy wywietlać wartosć 
        #ax.text(p.get_x()+p.get_width()/2.,0.98*height - 100,str(percentage) + '\%',ha="center")      <- procentową na wykresie
        s+=1
    #plt.gcf().subplots_adjust(bottom=0.35) <- dodanie na dole wykresu trochę miejsca
    #plt.tight_layout()  <- dopasowanie wykresu do całej powierzchni okna
    #fig = ax.get_figure()
    #fig.savefig("pfam.png")
    print(i , analysis_protein_count)
    plt.show()