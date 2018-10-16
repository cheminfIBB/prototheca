# -*- coding: utf-8 -*-
import numpy as np
import itertools
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f" ,"--fasta" , nargs = "+" , help = "FASTA files with organisms proteoms.")
parser.add_argument("-r" ,"--rep_names_faa" , nargs = '+' ,help = "Representative names of examined organisms.")
parser.add_argument("-b" ,"--blast_results" , nargs = '+' , help = "Files with BLAST results for each organism.")
parser.add_argument("-t", "--type", nargs = '?',default = "org", help = "org - Analyze BLAST results from blasting organisms against each other\n base - Analyze BLAST results from blasting organisms against sequence databases")
parser.add_argument("-db", "--databases", nargs = '*', help = "If analyzis type is base than here goes FASTA files of databases.")
parser.add_argument("-db_names", "--databases_names", nargs = '*', help = "Representative names of databases - these will be shown ")
parser.add_argument("-e" ,"--e_cutoff" , nargs = '?',default = 10 , help = "For example for e = 10, the cut off will be e^-10")
args = parser.parse_args()

names_faa = args.fasta
rep_names_faa = args.rep_names_faa
bl_names = args.blast_results
cutoff = float(args.e_cutoff)
kody = []
code_len = 2
while kody == []:
    for i in names_faa:
        name = i
        file_name = open(name, "r")
        kod = set()
        for line in file_name:
            new = line.strip().split(">")
            if len(new) == 2:
                kod.add(line[1:code_len])
        file_name.close()
        kody.append(kod)
    check = []
    for k in kody:
        check.append(set(k))
    check = list(itertools.combinations(check,2))
    for x in check:
        if len(x[0].intersection(x[1])) > 0:
            kody = []
            code_len += 1
            break
        
results = open("codes","w")
for i in range(len(kody)):
    results.write(">"+rep_names_faa[i] + " \n")
    for each in kody[i]:
        results.write(each + " \n")
results.close()

if args.type == "base":
    base_files = args.databases
    base_names = args.databases_names
    kody = []
    code_len = 2
    while kody == []:
        for i in base_files:
            name = i
            file_name = open(name, "r")
            kod = set()
            for line in file_name:
                new = line.strip().split(">")
                if len(new) == 2:
                    kod.add(line[1:code_len])
            file_name.close()
            kody.append(kod)
        check = []
        for k in kody:
            check.append(set(k))
        check = list(itertools.combinations(check,2))
        for x in check:
            if len(x[0].intersection(x[1])) > 0:
                kody = []
                code_len += 1
                break
    results = open("codes_bases","w")
    for i in range(len(kody)):
        results.write(">"+base_names[i] + " \n")
        for each in kody[i]:
            results.write(each + " \n")
    results.close()
            
            
results = open("proteomes_stats.txt","w")
proteomes_lengths = []
proteomes_average = []
s = 1

for i in names_faa:
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
    protein_names_faa = []
    for line in file_name:
        new = line.strip().split(">")
        if len(new) == 2:
            if protein_count != 0:
                #spisywanie długości białka
                protein_lengths.append(prot_len)
            #odczytanie nowego białka w pliku
            protein_count += 1
            prot_len = 0
            protein_names_faa.append(new[1].split(" ")[0])
        if len(new) == 1:
            prot_len += len(new[0])
            
    file_name.close()            
    protein_lengths.append(prot_len)        
    arr = np.array(protein_lengths)
    results.write("Proteins: " + str(protein_count) + "\n")
    results.write("Average_lenght: " + str(round(sum(protein_lengths)/protein_count)) + "\n")    
    proteomes_lengths.append(protein_count)
    proteomes_average.append(round(sum(protein_lengths)/protein_count))
results.close()

results = open("proteomes_stats.txt","r")
for line in results:
    print(line)
results.close()

fasta_names = ""
rep_names = ""
blast_names = ""
proteomes_count = ""

for i in range(len(names_faa)):
    fasta_names += str(names_faa[i]) + " "
    rep_names += "\"" + str(rep_names_faa[i]) + "\" "
    blast_names += str(bl_names[i]) + " "
    proteomes_count += str(proteomes_lengths[i]) + " "
    
if args.type == "org":
    print("Moving to blast_analysis script")
    os.system("py -3 blast_analysis.py -f {0} -r {1} -b {2} -e {3}".format(fasta_names,rep_names,blast_names,cutoff))
    
elif args.type == "base":
    base_string = ""
    base_filestring = ""
    for i in range(len(base_names)):
        base_string += str(base_names[i]) + " "
        base_filestring += str(base_files[i]) + " "
    os.system("py -3 base_blast_analysis.py -f {0} -r {1} -pc {2} -e {3} -db {4} -db_names {5}".format(fasta_names,rep_names,proteomes_count,cutoff,base_filestring,base_string))