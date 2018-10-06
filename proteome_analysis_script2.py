# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

# sprawdzić cborn

# nazwy plików:
# chlorella_variabilis.faa
# trichophyton_rubrum.faa
# auxenochlorella_protothecoides.faa
# helicosporidium_sp.faa
# candida_albicans.faa
names = ["chlorella_variabilis.faa","auxenochlorella_protothecoides.faa","helicosporidium_sp.faa",
"trichophyton_rubrum.faa","candida_albicans.faa"]
rep_names = ["Chlorella Variabilis","Auxenochlorella Protothecoides","Helicosporidium Sp","Trichophyton Rubrum","Candida Albicans"]
results = open("proteomes_stats.txt","w")
proteomes_lengths = []
proteomes_average = []
s = 1
for i in names:
    name = i
    results.write("> " + i + "\n")
    file_name = open(name, "r")
    #ilość białek w proteomie
    protein_count = 0
    #długości kolejnych białek proteomu
    protein_lengths = []
    #łączna długość białek w proteomie
    prot_len = 0
    #ref seq białka
    protein_names = []
    for line in file_name:
        new = line.strip().split(">")
        if len(new) == 2:
            if protein_count != 0:
                #spisywanie długości białka
                protein_lengths.append(prot_len)
                results.write(str(prot_len) + "\n")
            #odczytanie nowego białka w pliku
            protein_count += 1
            prot_len = 0
            protein_names.append(new[1].split(" ")[0])
            results.write(protein_names[-1] + " ")
        if len(new) == 1:
            prot_len += len(new[0])
        if len(new) == 0:
            raise Exception("Blank line in a file - please remove it.")
    results.write(str(prot_len) + "\n")
    file_name.close()            
    protein_lengths.append(prot_len)        
    arr = np.array(protein_lengths)
    results.write("Proteins: " + str(protein_count) + "\n")
    results.write("Average_lenght: " + str(round(sum(protein_lengths)/protein_count)) + "\n")    
    #print(protein_lengths.index(max(protein_lengths)))
    print(round(sum(protein_lengths)/protein_count), i)
    proteomes_lengths.append(protein_count)
    proteomes_average.append(round(sum(protein_lengths)/protein_count))

    #histogram przedstawiający liczbę białek proteomu 
    #w zależoności od długości sekwencji tych białek
    plt.figure(s)
    plt.hist(arr, bins=np.arange(0, 7500))
    plt.ylabel('Occurances')
    plt.xlabel('Length')
    plt.show()
    s+=1
    plt.figure(s)
    plt.hist(arr, bins=np.arange(0, 2000))
    plt.ylabel('Occurances')
    plt.xlabel('Length')
    plt.show()
    s+=1
    protein_lengths.sort()
    print(protein_lengths[-1:-10:-1])
print(proteomes_lengths)
plt.figure(s)
plt.bar([1,2,3,4,5],proteomes_lengths)
plt.ylabel('Number of proteins in proteome')
plt.xticks([1,2,3,4,5],rep_names,rotation=60)
plt.show()
results.close()