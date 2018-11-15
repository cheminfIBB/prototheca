from __future__ import unicode_literals
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from copy import deepcopy
import re

def tex_escape(text):
    """
        :param text: a plain text message
        :return: the message escaped to appear correctly in LaTeX
    """
    conv = {
        '&': r'\&',
        '%': r'\%',
        '$': r'\$',
        '#': r'\#',
        '_': r'\_',
        '{': r'\{',
        '}': r'\}',
        '~': r'\textasciitilde{}',
        '^': r'\^{}',
        '\\': r'\textbackslash{}',
        '<': r'\textless{}',
        '>': r'\textgreater{}',
    }
    regex = re.compile('|'.join(re.escape(str(key)) for key in sorted(conv.keys(), key = lambda item: - len(item))))
    return regex.sub(lambda match: conv[match.group()], text)

def make_analysis_dict(interpro_names, analysis):
    # legenda do wyników z interproscan
    headers = ["Sequence ID", "Sequence MD5 digest", "Sequence Length",
                       "Analysis", "Signature Accession", "Signature Description",
                       "Start Location", "Stop Location", "Score (e-value)", "Status",
                       "Run Date", "InterPro annotations description", "GO annotations",
                       "Pathways annotations"]
    
    analysis_dict = {}
    if analysis == []:
        for org in interpro_names:
            k = str(org)
            print("Creating analysis list from file " + str(k))
            p = pd.read_csv(k, sep='\t', header = None, names=headers)
            for i2 in range(1, p.shape[0]):           
                    row = p.iloc[[i2]].values[0]
                    if row[3] not in analysis:
                        analysis.append(row[3])    
    for i in analysis:
        analysis_dict[i] = set()
    
    return analysis_dict, analysis, headers

def load_interpro(names_faa, rep_names, interpro_names, analysis_dict, analysis, headers ):
    main_dict = {}
    whole_protein_count = []
    
    for i in range(len(rep_names)):
        main_dict[rep_names[i]] = deepcopy(analysis_dict)
        name = names_faa[i]
        file_name = open(name, "r")
        protein_names = []
        print("Loading fasta file for " + str(name))
        for line in file_name:      #robi listę białek których szukamy (mogą być to białka unikalne, może być cały proteom) na podstawie pliku fasta
            new = line.strip().split(">")
            if len(new) == 2:
                protein_names.append(new[1].split(" ")[0])
        file_name.close()
        whole_protein_count.append(len(protein_names))
        
        k = str(interpro_names[i])
        p = pd.read_csv(k, sep='\t', header = None, names=headers) 
        print("Loading InterProScan results for " + str(k))        #wczytuje wyniki interproscan i sprawdza linia po linii w poszukiwaniu      
        for i2 in range(1, p.shape[0]):                            #interesujących nas białek i analiz
                row = p.iloc[[i2]].values[0]
                if row[0] in protein_names and row[3] in analysis:
                    main_dict[rep_names[i]][row[3]].add(row[0])
        print('Loading ' + rep_names[i] + ' done')
    return main_dict, whole_protein_count

def make_plots(rep_names, use_tex, analysis, main_dict, whole_protein_count):
    if use_tex == True:
        plt.rc('text', usetex=True)
        rep_names2 = []
        #nazwy organizmów kursywą
        for i in rep_names:
            rep_names2.append(r'\textit{' + tex_escape(i) + '}')
    else:
        plt.rc('text', usetex=False)
        rep_names2 = deepcopy(rep_names)
        
    print("Generating plots...")    
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    rr = 0
    for a_name in analysis:    # pętla robiąca wykres słupkowy dla każdej analizy, każdy organizm 
        plt.figure(rr)    # na wykresie zawiera liczbę wszystkich białek oraz liczbę białek z jego proteomu mające hit w czasie analizy 
        plt.figure(figsize=(20,10))
        rr += 1
        analysis_protein_count = []
        for name in rep_names:
            analysis_protein_count.append(len(main_dict[name][a_name]))
        sns.set()
        sns.set_style("whitegrid")
        dicti = dict()
        dicti["proteomes"] = rep_names2
        dicti["whole"] = whole_protein_count
        dicti2 = dict()
        dicti2["proteomes"] = rep_names2
        dicti2[a_name] = analysis_protein_count
        df = pd.DataFrame(dicti)
        df2 = pd.DataFrame(dicti2)
        sns.set(font_scale = 2.4)
        ap = sns.barplot(x = rep_names2,y = 'whole', color = "skyblue",data = df)
        ax = sns.barplot(x = rep_names2,y = a_name,color="blue",data = df2)
        ax.set_ylabel("Number of proteins with hit from " + tex_escape(a_name), fontsize=24, alpha = 0.8)
        searcher = whole_protein_count + analysis_protein_count
        s = 0
        for p in ax.patches:
            ax.text(p.get_x()+p.get_width()/2.,p.get_height()+1,str(searcher[s]),ha="center")
            #percentage = str(round(round(1-(analysis_protein_count[s]/whole_protein_count[s]),4)*100,2))  <- czy wywietlać wartosć 
            #ax.text(p.get_x()+p.get_width()/2.,0.98*height - 100,str(percentage) + '\%',ha="center")      <- procentową na wykresie
            s+=1
        plt.gcf().subplots_adjust(bottom=0.1)
        plt.savefig(str(a_name)+"_plot",bbox_inches="tight",pad_inches=0.5)
        plt.close()   