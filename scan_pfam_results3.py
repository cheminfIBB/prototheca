from __future__ import unicode_literals
import pandas as pd
#import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import copy

names = ["prototheca.fasta.tsv", "chlorella_variabilis.faa.tsv","auxenochlorella_protothecoides.faa.tsv","helicosporidium_sp.faa.tsv",
"trichophyton_rubrum.faa.tsv","candida_albicans.faa.tsv"]
rep_names = ["Prototheca", "Chlorella variabilis","Auxenochlorella protothecoides","Helicosporidium sp","Trichophyton rubrum","Candida albicans"]
rep_names2 = ["\it{Prototheca}","\it{Chlorella} \n \it{variabilis}","\it{Auxenochlorella} \n \it{protothecoides}","\it{Helicosporidium sp.}","\it{Trichophyton} \n \it{rubrum}","\it{Candida} \n \it{albicans}"]
names_faa = ["prototheca.fasta.unique","chlorella_variabilis.faa.unique","auxenochlorella_protothecoides.faa.unique","helicosporidium_sp.faa.unique","trichophyton_rubrum.faa.unique","candida_albicans.faa.unique"]

#["prototheca.fasta.unique","chlorella_variabilis.faa.unique","auxenochlorella_protothecoides.faa.unique","helicosporidium_sp.faa.unique","trichophyton_rubrum.faa.unique","candida_albicans.faa.unique"]
#["prototheca.fasta","chlorella_variabilis.faa","auxenochlorella_protothecoides.faa","helicosporidium_sp.faa","trichophyton_rubrum.faa","candida_albicans.faa"]
#names_tsv = ["prototheca.fasta.tsv_pfam.csv","chlorella_variabilis.faa.tsv_pfam.csv","auxenochlorella_protothecoides.faa.tsv_pfam.csv","helicosporidium_sp.faa.tsv_pfam.csv","trichophyton_rubrum.faa.tsv_pfam.csv","candida_albicans.faa.tsv_pfam.csv"]

# przeprowadzone analizy
analysis = ['TMHMM','SignalP_EUK','Coils','Phobius','MobiDBLite']
#['PIRSF', 'TIGRFAM', 'TMHMM', 'PANTHER', 'SMART', 'SUPERFAMILY', 'Phobius', 'Pfam', 'Gene3D',
#            'ProSitePatterns', 'MobiDBLite', 'PRINTS', 'ProSiteProfiles', 'SFLD', 'SignalP_EUK', 'Coils', 'Hamap',
#            'ProDom', 'CDD']

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
    #ref seq białka
    protein_names = []
    for line in file_name:
        new = line.strip().split(">")
        if len(new) == 2:
            protein_names.append(new[1].split(" ")[0])
    file_name.close()
    whole_protein_count.append(len(protein_names))
    
    k = str(names[i])
    p = pd.read_csv(k, sep='\t', header = None, names=headers)
    for i2 in range(1, p.shape[0]):
            row = p.iloc[[i2]].values[0]
            if row[0] in protein_names and row[3] in analysis:
                maine_dictus[rep_names[i]][row[3]].add(row[0])
                
    #for x in analysis:      
    #    ll = str(k)+ "_" + str(x) + ".tsv"
    #    results = open(ll, "w")
    #    results.write(rep_names[i] + "\n")
    #    for pp in headers:
    #        results.write(str(pp) + "\t")
    #    p = pd.read_csv(k, sep='\t', header=None, names=headers)
    #    for iii in range(0, p.shape[0]):
    #        row = p.iloc[[iii]].values[0]
    #        if row[3] == x and row[0] in protein_names:
    #            for val in row:
    #                results.write(str(val) + "\t")
    #            results.write("\n")
    #    results.close()      
    print('Loading ' + rep_names[i] + ' done')


plt.rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
rr = 0
for i in analysis:    
    plt.figure(rr)
    rr += 1
    #whole_protein_count = [733,2512, 1065, 1051, 4039, 1845]
    analysis_protein_count = []
    for name in rep_names:
        analysis_protein_count.append(len(maine_dictus[name][i]))
    sns.set()
    #nazwy kursywą
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
        #percentage = str(round(round(1-(analysis_protein_count[s]/whole_protein_count[s]),4)*100,2))
        #ax.text(p.get_x()+p.get_width()/2.,0.98*height - 100,str(percentage) + '\%',ha="center")
        s+=1
    #plt.gcf().subplots_adjust(bottom=0.35)
    #plt.tight_layout()
    #fig = ax.get_figure()
    #fig.savefig("pfam.png")
    print(i , analysis_protein_count)
    plt.show()