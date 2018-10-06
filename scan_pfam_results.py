from __future__ import unicode_literals
import pandas as pd
#import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
"""
Opisy analiz:




"""
names = ["chlorella_variabilis.faa.tsv","auxenochlorella_protothecoides.faa.tsv","helicosporidium_sp.faa.tsv",
"trichophyton_rubrum.faa.tsv","candida_albicans.faa.tsv"]
rep_names = ["Chlorella Variabilis","Auxenochlorella Protothecoides","Helicosporidium Sp","Trichophyton Rubrum","Candida Albicans"]

# przeprowadzone analizy
analysis = ['PIRSF', 'TIGRFAM', 'TMHMM', 'PANTHER', 'SMART', 'SUPERFAMILY', 'Phobius', 'Pfam', 'Gene3D',
                'ProSitePatterns', 'MobiDBLite', 'PRINTS', 'ProSiteProfiles', 'SFLD', 'SignalP_EUK', 'Coils', 'Hamap',
                'ProDom', 'CDD']
    # legenda do wyników z interproscan
headers = ["Sequence ID", "Sequence MD5 digest", "Sequence Length",
                   "Analysis", "Signature Accession", "Signature Description",
                   "Start Location", "Stop Location", "Score (e-value)", "Status",
                   "Run Date", "InterPro annotations description", "GO annotations",
                   "Pathways annotations"]

#whole_protein_count = []
#pfam_protein_count = []
"""
for each in names:
    k = str(each)+ "_" + 'Pfam' + ".csv"
    p = pd.read_csv(k, sep='\t', header = None, names=headers)
    seen_protein = 0
    name = None
    for i in range(1, p.shape[0]):
            row = p.iloc[[i]].values[0]            
            if row[0] != name:
                name = row[0]
                seen_protein += 1
    pfam_protein_count.append(seen_protein)
    p = pd.read_csv(each, sep='\t', header = None, names=headers)
    seen_protein = 0
    name = None
    for i in range(0, p.shape[0]):
            row = p.iloc[[i]].values[0]            
            if row[0] != name:
                name = row[0]
                seen_protein += 1
    whole_protein_count.append(seen_protein)
"""
plt.rc('text', usetex=True)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
whole_protein_count = [5930,9780, 7014, 6033, 8706, 6030]
pfam_protein_count = [4124,7357, 4896, 3923, 5899, 4823]
#nazwy kursywą
#nazwy kursywą
rep_names = ["\it{Prototheca}","\it{Chlorella} \n \it{variabilis}","\it{Auxenochlorella} \n \it{protothecoides}","\it{Helicosporidium sp.}","\it{Trichophyton} \n \it{rubrum}","\it{Candida} \n \it{albicans}"]
sns.set_style("whitegrid")
dicti = dict()
dicti["proteomes"] = rep_names
dicti["whole"] = whole_protein_count
dicti2 = dict()
dicti2["proteomes"] = rep_names
dicti2["pfam"] = pfam_protein_count
df = pd.DataFrame(dicti)
df2 = pd.DataFrame(dicti2)
print(df)
print(df2)
#sns.barplot(data = df)
sns.set(font_scale = 2.4)
sns.barplot(x = rep_names,y = 'whole', color = "skyblue",data = df)
ax = sns.barplot(x=rep_names,y='pfam',color="blue",data=df2)
ax.set_ylabel("Number of proteins",fontsize=24, alpha=0.8)
s = 0
for p in ax.patches[0:int((len(ax.patches)/2))]:
    if s<6:
        percentage = str(round(round(1-(pfam_protein_count[s]/whole_protein_count[s]),4)*100,2))
        height = p.get_height()
        ax.text(p.get_x()+p.get_width()/2.,height-500,str(percentage) + '\%',ha="center")
        s+=1
plt.gcf().subplots_adjust(bottom=0.35)
plt.tight_layout()
fig = ax.get_figure()
fig.savefig("pfam.png")
