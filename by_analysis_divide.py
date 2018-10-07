import pandas as pd

names = ["chlorella_variabilis.faa.tsv","trichophyton_rubrum.faa.tsv","auxenochlorella_protothecoides.faa.tsv",
"helicosporidium_sp.faa.tsv","candida_albicans.faa.tsv"]
rep_names = ["Chlorella Variabilis", "Trichophyton Rubrum", "Auxenochlorella Protothecoides", "Helicosporidium Sp",
             "Candida Albicans"]
             
for l in names:
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
    file_name = open(l, "r")
    for x in analysis:
        k = str(l)+ "_" + str(x) + ".csv"
        results = open(k, "w")
        results.write(rep_names[0] + "\n")
        p = pd.read_csv(l, sep='\t', header=None, names=headers)
        for i in range(0, p.shape[0]):
            row = p.iloc[[i]].values[0]
            if row[3] == x:
                for val in row:
                    results.write(str(val) + "\t")
                results.write("\n")
        results.close()      
    del (rep_names[0])
