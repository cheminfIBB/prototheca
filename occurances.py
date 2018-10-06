import pandas as pd

names = ["chlorella_variabilis.faa.tsv"]
#,"trichophyton_rubrum.faa.tsv","auxenochlorella_protothecoides.faa.tsv","helicosporidium_sp.faa.tsv","candida_albicans.faa.tsv"]
rep_names = ["Chlorella Variabilis","Trichophyton Rubrum","Auxenochlorella Protothecoides","Helicosporidium Sp","Candida Albicans"]
results = open("interpro_stats.txt","w")
results2  = open("interpro_hits.txt","w")
for i in names:
    file_name = open(i, "r")
    results.write("> " + rep_names[0] + "\n")
    results2.write("> " + rep_names[0] + "\n")
    del(rep_names[0])
    #legenda do wyników z interproscan
    
    headers = ["Sequence ID", "Sequence MD5 digest", "Sequence Length",
               "Analysis", "Signature Accession", "Signature Description",
               "Start Location", "Stop Location", "Score (e-value)", "Status",
               "Run Date", "InterPro annotations description","GO annotations",
               "Pathways annotations"]
    
    #przeprowadzone analizy
    
    analysis = {'PIRSF', 'TIGRFAM', 'TMHMM', 'PANTHER', 'SMART', 'SUPERFAMILY', 'Phobius', 'Pfam', 'Gene3D',
                'ProSitePatterns', 'MobiDBLite', 'PRINTS', 'ProSiteProfiles', 'SFLD', 'SignalP_EUK', 'Coils', 'Hamap',
                'ProDom', 'CDD'}
    
    p = pd.read_csv(i, sep='\t', header = None, names=headers)
    occurances = dict()  #np do średniej ilości trafień na jedno białko
    occurances_unique = dict()
    unique_occurances = dict()  # to samo tylko unikalnych trafień
    high_score = 0
    for i in analysis:
        occurances[i] = 0
        occurances_unique[i] = 0
        unique_occurances[i] = 0
    
    length = 0  # np do średniej długości hitu
    length_list = []  #np do mediany, najdłuższego i najkrótszego hitu
    
    seen_proteins = set() # lista wszystkich widzianych w pliku białek
    
    hits = []  #lista trafień (id,analiza, start,stop, e-value)
    longer_hits = []
    for i in range(0,p.shape[0]):
        row = p.iloc[[i]].values[0]
        if row[0] not in seen_proteins:         #KOLUMNY
            seen_analysis = set()
            seen_proteins.add(row[0])           #1 przeszliśmy do nowego białka
            unique_occurances[row[3]] += 1      #4 dodanie do konkretnej analizy jej unikatowego wystąpenia
            seen_analysis.add(row[3])           # zapisanie, że widzieliśmy daną analizę dla tego białka
            length += row[2]                    #3 dodaje długosć białka do sumy długości wszystkich białek
            length_list.append(row[2])
        else:
            if row[3] not in seen_analysis:
                unique_occurances[row[3]] += 1
                seen_analysis.add(row[3])
        occurances[row[3]] += 1
        if row[8] == "-":  # Bada ile trafień ma małe e-value (na podstawie potęgi)
            high_score += 1
        else:
            e_power = row[8].split('-')
            if len(e_power)>1:
                if int(e_power[1]) >= 5:
                    high_score += 1
        k = (row[0],row[3],row[4])
        k2 = (row[6],row[7],row[8])    #zapisuje ref_seq,analizę,start,stop
        if k not in hits or k != hits[-1]:
            occurances_unique[row[3]] += 1
            results.write(row[0] + " " + row[3] + " " + row[4] + "\n")
        results.write(str(row[6]) + " " + str(row[7]) + " " + str(row[8]) + "\n")
        hits.append(k)
        longer_hits.append(k+k2)
        results2.write(str(row[0]) + " " + str(row[3]) + " " + str(row[4]) + " " + str(row[6]) + " " + str(row[7]) + " " + str(row[8]) + "\n")
    print("Ilosć trafień w poszczególnych bazach: ",occurances)
    print("Ilosć unikalnych trafień w poszczególnych bazach: ",unique_occurances)
    print("Łączna długość wszystkich trafień: ",length,"Ilość białek w których zanotowano trafienie: ",len(length_list))
    print("Średnia długość na trafienie: ",round(length/len(length_list)),"Najdłuższe/najkrótsze trafienie: ",max(length_list),min(length_list))
    #wykres ile trafień znajduje się w jakich przedziałach?
    
    print("Baza z której pochodzi najwięcej trafień oraz liczba tych trafień: ",max(occurances, key=occurances.get),max(occurances.values()))
    print("Baza z której pochodzi najwięcej unikalnych trafień oraz liczba tych trafień: ",max(occurances_unique, key=occurances_unique.get),max(occurances_unique.values()))
    
    #wykres z tym ile trafień z której bazy?
    
    print("Hits with e-value above",5,":",high_score)
    print(headers)
    unique_hits = set(hits)
    print("Unique Hits", len(unique_hits))
    #print(p.loc[p['Analysis'].isin(['PANTHER'])]) - zwraca wszystkie hity o wskazanej wartości atrybutu
    
    # p [[kolumna][wiersz]]
    #p.iloc[ albo zakres np 1:10 albo lista[1,2,3...], i tak i dla kolumn i, po przecinku, dla wierszy
    # print(p.iloc[[100,101], 2:4].values)

results.close()