# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

# nazwy plików:
# chlorella_variabilis.faa
# trichophyton_rubrum.faa
# auxenochlorella_protothecoides.faa
# helicosporidium_sp.faa
# candida_albicans.faa
names = ["CAZyDB.fasta","merops_scan.fasta","phi_accessions.fasta"]
rep_names = ["CAZY","MEROPS","PhiBase"]

kody = []
for i in names:
    name = i
    file_name = open(name, "r")
    kod = set()
    protein_count = 0
    for line in file_name:
        new = line.strip().split(">")
        if len(new) == 2:
            kod.add(line[1:6])
            protein_count += 1
    file_name.close()
    kody.append(kod)
#for i in kody:
#    print(i)
      
results = open("codes","w")
for i in range(len(kody)):
    results.write(rep_names[i] + " \n")
    for each in kody[i]:
        results.write(each + " \n")
results.close()

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
    print(protein_count)
   
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

plt.figure(s)
plt.bar([1,2,3],proteomes_lengths)
plt.ylabel('Number of sequences in base')
plt.xticks([1,2,3],rep_names,rotation=60)
plt.show()
results.close()