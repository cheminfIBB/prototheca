# -*- coding: utf-8 -*-
import itertools

def make_protein_codes(names_faa,rep_names,a_type="",base_files="",base_names=""):
    codes = []
    code_len = 2
    while codes == []:
        for name in names_faa:
            file_name = open(name, "r")
            kod = set()
            for line in file_name:
                new = line.strip().split(">")
                if len(new) == 2:
                    kod.add(line[1:code_len])
            file_name.close()
            codes.append(kod)
        check = []
        for k in codes:
            check.append(set(k))
        check = list(itertools.combinations(check, 2))
        for x in check:
            if len(x[0].intersection(x[1])) > 0:
                codes = []
                code_len += 1
                break
            
    results = open("codes", "w")
    for i in range(len(codes)):
        results.write(">"+rep_names[i] + " \n")
        for each in codes[i]:
            results.write(each + " \n")
    results.close()
    
    if a_type == "base":
        codes = []
        code_len = 2
        while codes == []:
            for name in base_files:
                file_name = open(name, "r")
                kod = set()
                for line in file_name:
                    new = line.strip().split(">")
                    if len(new) == 2:
                        kod.add(line[1:code_len])
                file_name.close()
                codes.append(kod)
            check = []
            for k in codes:
                check.append(set(k))
            check = list(itertools.combinations(check, 2))
            for x in check:
                if len(x[0].intersection(x[1])) > 0:
                    codes = []
                    code_len += 1
                    break
        results = open("codes_bases", "w")
        for i in range(len(codes)):
            results.write(">"+base_names[i] + " \n")
            for each in codes[i]:
                results.write(each + " \n")
        results.close()

def make_proteomes_stats(names_faa):
    results = open("proteomes_stats.txt", "r+")
    proteomes_lengths = []
    proteomes_average = []
    for name in names_faa:
        results.write("> " + name + "\n")
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
        results.write("Proteins: " + str(protein_count) + "\n")
        results.write("Average_lenght: " + str(round(sum(protein_lengths)/protein_count)) + "\n")    
        proteomes_lengths.append(protein_count)
        proteomes_average.append(round(sum(protein_lengths)/protein_count))
    
    results.seek(0)
    for line in results:
        print(line)   
    results.close()
    return proteomes_lengths