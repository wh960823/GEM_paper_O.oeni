"""
This script was used to generate 20 combinations of metabolome data for mixOmics analyses

This script was for paper 'Integrative omics and genome-scale metabolic modeling approaches
reveal mechanisms of Oenococcus oeni combating acidic stress' 

"""

from itertools import combinations

# Input file 'metabolome_data.txt' cotains the full metabolome data (all samples)
with open("metabolome_data.txt", "r") as f1:
    signal_dist = {}
    for each_line in f1.readlines()[1:]:
        metabolite_id = each_line.strip().split("\t")[0]
        samples = len(each_line.strip().split("\t"))
        for i in range(1, samples):
            if metabolite_id not in signal_dist:
                signal_dist[metabolite_id] = []
                signal_dist[metabolite_id].append(each_line.strip().split("\t")[i])
            else:
                signal_dist[metabolite_id].append(each_line.strip().split("\t")[i])

# Now generate C(6,3) combinations for test MixOmics integration
# 1-6: pH4.8_0h; 7-12: pH4.8_1h; 13-18: pH4.8_3h; 19_24: pH3.0 1h; 25-30: pH3.0 3h

a_list = [1, 2, 3, 4, 5, 6]
q = 1
for each_combination in list(combinations(a_list, 3)):
    with open("metabolite_samples_{0}.txt".format(str(q)), "w+") as outfile: # outfile: 20 combinations of metabolome data
        for each_metabolite in signal_dist.keys():
            outfile.writelines(each_metabolite + "\t")
            select_list = []
            for k in range(len(each_combination)):
                sample_id = each_combination[k]
                for m in range(0, 30, 6):
                    select_list.append(sample_id + m)
                    select_list.sort()
            add_line = ""
            for l in select_list:
                add_line += signal_dist[each_metabolite][l - 1]
                add_line += "\t"
            add_line = add_line.strip()
            outfile.writelines(add_line + "\n")
    q += 1
