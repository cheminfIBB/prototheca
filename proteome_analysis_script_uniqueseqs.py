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