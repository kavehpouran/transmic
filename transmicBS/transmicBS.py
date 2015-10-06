# transmicBS.py
#
# A program for hierarchical clustering of genetic sequences
# on phylogenetic trees. The code is customised for the identification of 
# putative transmission clusters of viruses. In particular, it allows 
# for multiple follow-up viral sequences extracted from the same individual.
#
# Kaveh Pouran Yousef, Sylvana Gromoeller, 2015
 
from transmiclib import *
import os
    
#Parse the control file of the program 
input_filename =  "controlfile.tbsc"
paramlist = parse_input_file(input_filename)
outgroupSeqId = paramlist[7]
 
#Read the phylogenetic tree
print("Reading the tree and rooting it using the outgroup sequence.")
filename1 = paramlist[0]
f = open(filename1, "r") 
s = f.read()
f.close()
tree=dendropy.Tree.get_from_string(s, schema='newick')

tree.is_rooted = True
tree.to_outgroup_position(tree.find_node_with_taxon_label(string.replace(outgroupSeqId,"_"," ")))

pdm=treecalc.PatristicDistanceMatrix(tree)

#number of sequences in the tree
taxonlist = tree.infer_taxa()
taxon_labels = []
for t in taxonlist:
    taxon_labels.append(string.replace(t.label," ","_"))
nr_leaves = len(taxon_labels)    

#Read the cut-offs for bootstrap support and mean patristic distance from standard input
patrDist_cutoff = paramlist[1]
bootstrap_cutoff = paramlist[2]

#Set pointer to the root-node
start_node = tree.seed_node
 
#determine significant nodes using the bootstrap edge labels
is_node_significant_raxml(start_node, bootstrap_cutoff)
         	
#Compute the clusters
patWiseClust = paramlist[3]

print("Hierarchical clustering.")
if (patWiseClust):
	#Inter-patient sequence-wise cluster computation
    clusters_meandist = flatten(get_clusters_meandist_interPatient(start_node, patrDist_cutoff, pdm, outgroupSeqId))
else:
    clusters_meandist = flatten(get_clusters_meandist(start_node, patrDist_cutoff))

#Set output file name for text output
outdir1="text_out";
outfilename = outdir1+"/TransmicBS_clusters_out"+str(patrDist_cutoff)+"_bootstrap"+str(bootstrap_cutoff)+".txt"

#generate output directory
if not os.path.exists(outdir1):
    os.makedirs(outdir1)

#Write clusters to the text file
clusters_meandist = write_clusters_to_file(clusters_meandist, outfilename, nr_leaves, patrDist_cutoff, bootstrap_cutoff)
nr_clust = len(clusters_meandist)  
print('Text summary written to'), outfilename+"." 
 
#Generate an excel table with clusters (requires another excel file with sequence infos)
#RKI database file
write_out_table = paramlist[4]
if (write_out_table):
    outdir2="table_out";
    if not os.path.exists(outdir2):
        os.makedirs(outdir2)
            	
    xls_table_infilename = paramlist[5] 
    enforce_RKI_format = paramlist[6]

    xls_outfilename = outdir2+"/combined_output_clustering_mppd"+str(patrDist_cutoff)+"_bootstrap"+str(bootstrap_cutoff)+"_all.tsv"
    write_clusters_to_table(clusters_meandist, xls_table_infilename, xls_outfilename,taxon_labels, outgroupSeqId, enforce_RKI_format)
    print('Table summary of clustering results written to'), xls_outfilename +"." 
