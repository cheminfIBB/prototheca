# -*- coding: utf-8 -*-

import argparse
from copy import deepcopy
import blast_analysis, base_blast_analysis, interproscan_analysis, proteome_analysis

parser = argparse.ArgumentParser()
parser.add_argument("-f" ,"--fasta" , nargs = "+" , 
                    help = "FASTA files with organisms proteoms.")

parser.add_argument("-r" ,"--rep_names" , nargs = '+' ,
                    help = "Representative names of examined organisms.")

parser.add_argument("-b" ,"--blast_results" , nargs = '*',default = [] ,
                    help = "Files with BLAST results for each organism.")

parser.add_argument("-t", "--type", nargs = '?',default = "org", 
                    help = "org - Analyze BLAST results from blasting organisms against each other \n base - Analyze BLAST results from blasting organisms against sequence databases\n inter - Analyze InterProScan results for given organisms")

parser.add_argument("-ips" ,"--interproscan" , nargs = '*' ,default = [], 
                    help = "InterProScan results files for organisms.")

parser.add_argument("-a" ,"--interproscan_analysis" , nargs = '*' ,default = [], 
                    help = "Analysis to be considered from InterProScan results.")

parser.add_argument("-db", "--databases", nargs = '*', 
                    help = "If analysis type is base than here goes FASTA files of databases.")

parser.add_argument("-db_names", "--databases_names", nargs = '*', 
                    help = "Representative names of databases - these will be shown ")

parser.add_argument("-e" ,"--e_cutoff" , nargs = '?',default = 10 , 
                    help = "For example for e = 10, the cut off will be e^-10")

parser.add_argument("-tex" ,"--use_tex",default = False,  type = lambda x: (str(x).lower() == 'true'), 
                    help = "By default False, if set to True, organisms names on plots will be italicized with LaTeX.")
args = parser.parse_args()

names_faa = args.fasta #FASTA files
rep_names = args.rep_names #organism names
blast_names = args.blast_results #files with blast results
cutoff = float(args.e_cutoff) #cut-off value for e-value
use_tex = args.use_tex #using LaTeX on plots
a_type = args.type #type of analysis to execute

whole_protein_count = proteome_analysis.make_proteomes_stats(names_faa) #preotomes statistics generation

if a_type == 'org':
    proteome_analysis.make_protein_codes(names_faa,rep_names)
    codes = blast_analysis.load_org_codes()
    organisms = blast_analysis.read_proteins(blast_names, deepcopy(names_faa), rep_names, codes, cutoff)
    blast_analysis.protein_list_to_set(organisms)
    blast_analysis.make_validation(organisms, rep_names)
    blast_analysis.save_validated_and_unique(organisms, names_faa, rep_names)
    
if a_type == 'base':
    base_files = args.databases
    base_names = args.databases_names
    proteome_analysis.make_protein_codes(names_faa,rep_names,a_type,base_files,base_names)
    base_codes, codes_org = base_blast_analysis.load_org_and_base_codes(names_faa, base_names, rep_names)
    organisms = base_blast_analysis.read_proteins(
             deepcopy(names_faa), rep_names, 
             base_names, base_codes, codes_org, cutoff )
    base_blast_analysis.protein_list_to_set(organisms)
    results = base_blast_analysis.make_plots(organisms, rep_names, use_tex)
    base_blast_analysis.save_results(results, rep_names, base_names)
    
if a_type == 'inter':
    inter_names = args.interproscan
    analysis = args.interproscan_analysis
    analysis_dict, analysis, headers = interproscan_analysis.make_analysis_dict(inter_names, analysis)
    main_dict, whole_protein_count = interproscan_analysis.load_interpro(names_faa, rep_names, inter_names, analysis_dict, analysis, headers)
    interproscan_analysis.make_plots(rep_names, use_tex, analysis, main_dict, whole_protein_count)