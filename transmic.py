# transmicBS.py
#
# A program for hierarchical clustering of genetic sequences
# on phylogenetic trees. The code is customised for the identification of 
# putative transmission clusters of viruses. In particular, it allows 
# for multiple follow-up viral sequences extracted from the same individual.
#
# Kaveh Pouran Yousef, Silvana Gromoeller, 2015
 
from transmiclib import *

#from paperReviewLib import * #PAPER REVIEW FUNCTIONS 
import sys
import os
    
#Parse the control file of the program 
input_filename =  "controlfile.tbsc"
paramlist = parse_input_file(input_filename)

#Read the phylogenetic tree
filename1 = paramlist[0]
tree_format = paramlist[1]
node_support = paramlist[2]

#distance cut-off
patrDist_cutoff_str = paramlist[3]
sep_idx = [pos for pos, char in enumerate(patrDist_cutoff_str) if char == ":"]
if len(sep_idx)==2:
    try:
        patrDist_cutoff1 = float(patrDist_cutoff_str[0:sep_idx[0]])
        patrDist_step = float(patrDist_cutoff_str[sep_idx[0]+1:sep_idx[1]])
        patrDist_cutoff2 = float(patrDist_cutoff_str[sep_idx[1]+1:])
        patrDistList=arange(patrDist_cutoff1, patrDist_cutoff2, patrDist_step)
    except:
        raise ValueError('controlfile error: distance cutoff range probably not a numeric value.')
else:
    try:
        patrDist_cutoff = float(patrDist_cutoff_str)
        patrDistList=[patrDist_cutoff]
    except:      
        raise ValueError('controlfile error: distance cutoff probably not a numeric value.')

#support cut-off  
support_cutoff_str = paramlist[4]
sep_idx = [pos for pos, char in enumerate(support_cutoff_str) if char == ":"]
if len(sep_idx)==2:
    try:
        support_cutoff1 = float(support_cutoff_str[0:sep_idx[0]])
        support_step = float(support_cutoff_str[sep_idx[0]+1:sep_idx[1]])
        support_cutoff2 = float(support_cutoff_str[sep_idx[1]+1:])
        supportList=arange(support_cutoff1, support_cutoff2, support_step)
    except:
        raise ValueError('controlfile error: support cutoff range probably not a numeric value.')
else:
    try:
        support_cutoff = float(support_cutoff_str)
        supportList=[support_cutoff]
    except:       
        raise ValueError('controlfile error: support cutoff probably not a numeric value.')


patWiseClust = True if paramlist[5]=="true" else False
xls_table_infilename = paramlist[6]
outgroupSeqId = paramlist[7]
figure_distances= True if paramlist[8]=="true" else False
figure_nr_clusters= True if paramlist[9]=="true" else False
figure_sizeof_clusters= True if paramlist[10]=="true" else False

print("Reading the tree and rooting it.")
f = open(filename1, "r") 
s = f.read()
f.close()
tree=dendropy.Tree.get_from_string(s, schema=tree_format)

tree.is_rooted = True
tree.update_bipartitions()
#pdm=treecalc.PatristicDistanceMatrix(tree) V3.12
pdm = dendropy.PhylogeneticDistanceMatrix.from_tree(tree) # until here

#number of sequences in the tree
taxonlist = tree.update_taxon_namespace()
taxon_labels = []
for t in taxonlist:
    taxon_labels.append(string.replace(t.label," ","_"))
nr_leaves = len(taxon_labels)    

#Set pointer to the root-node
start_node = tree.seed_node

#determine significant nodes using the bootstrap edge labels
#is_node_significant_raxml(start_node, bootstrap_cutoff)

#generate output directories
outdir1="text_out";
if not os.path.exists(outdir1):
    os.makedirs(outdir1)

outdir2="table_out";
if not os.path.exists(outdir2):
    os.makedirs(outdir2)

#defines how many columns of the sequence database file xls_table_infilename should be skipped 
#when writing the data
columns_skip=0;

#generate sequence database file if it is not given as parameter
if not xls_table_infilename:
	columns_skip=1
	xls_table_infilename="__auxiliary__seqID__file.xls"
	wb = xlwt.Workbook()
	ws = wb.add_sheet('Sheet 1')
	ws.write(0, 0, "SeqID")
	for j in range(len(taxon_labels)):
		idstr = taxon_labels[j]
		ws.write(j+1, 0, idstr)
	wb.save(xls_table_infilename)

#conduct clustering and write out the results
cluster_size_dist_mean=[];
nr_clusters=[];
print('Conducting clustering for given thresholds.')
for support_cutoff in supportList:
    cluster_size_dist_mean_tmp=[];
    nr_clusters_tmp=[];
    sys.stdout.write(node_support+' support cut-off: ' + str(support_cutoff)+'%\n')
    sys.stdout.write('distance cut-off: ')  
    for patrDist_cutoff in patrDistList:
        sys.stdout.write(str(patrDist_cutoff)+' ')
        assess_node_support(start_node, support_cutoff, node_support) #------------------------------------------
		#Compute the clusters
        if (patWiseClust):
            #Inter-patient sequence-wise cluster computation
            clusters_meandist= flatten(get_clusters_meandist_interPatient (start_node, patrDist_cutoff, outgroupSeqId, pdm))	
        else:
            clusters_meandist = flatten(get_clusters_meandist(start_node, patrDist_cutoff)) 
		    #clusters_meandist = flatten(get_clusters_meandist(start_node, patrDist_cutoff)) # --> this works with dendropy 4
		    #Set output file name for text output
        text_outfilename = outdir1+"/TransmicBS_clusters_out"+str(patrDist_cutoff)+"_bootstrap"+str(support_cutoff)+".txt"
        xls_outfilename = outdir2+"/combined_output_clustering_mppd"+str(patrDist_cutoff)+"_bootstrap"+str(support_cutoff)+"_all.tsv"
        taxon_labels_copy=  list(taxon_labels)
        write_clusters_to_file(clusters_meandist, text_outfilename, nr_leaves, patrDist_cutoff, support_cutoff)
        write_clusters_to_table(clusters_meandist, xls_table_infilename, xls_outfilename,taxon_labels_copy, outgroupSeqId, columns_skip)
        cluster_size_dist_mean_tmp.append(mean(array(map(len, clusters_meandist))))
        nr_clusters_tmp.append(len(clusters_meandist))
    cluster_size_dist_mean.append(array(cluster_size_dist_mean_tmp));
    nr_clusters.append(array(nr_clusters_tmp));    
        
        
        
    sys.stdout.write('\n')

print('Text-based clustering results written to directory:'), outdir1+"." 
print('Table-based clustering summary written to'), outdir2 +"." 

if xls_table_infilename=="__auxiliary__seqID__file.xls":
    os.remove(xls_table_infilename)

if (figure_distances):
	#visualization 1
	p1=plt.figure("Distribution of patristic distances in the phylogenetic tree")
	distvec=[]
	for i in range(len(taxonlist)):
		for j in range(i+1,len(taxonlist)):
			distvec.append(pdm(taxonlist[i],taxonlist[j]))
			
	#num_bins = ceil(max(distvec))
	# the histogram of the data
	num_bins=(30)
	n, bins, patches = plt.hist(distvec, num_bins, normed=1,  facecolor=[0,0.2,1], alpha=0.5)
	# add a 'best fit' line
	
	plt.xlabel('Pairwise patristic distance between taxa')
	plt.ylabel('Probability')
	#plt.title(r'Distribution of patristic distances in the phylogenetic tree.')
	
	# Tweak spacing to prevent clipping of ylabel
	p1.subplots_adjust(left=0.15)
	
#visualization 2
if (figure_nr_clusters):
	p2 = plt.figure("Number of clusters")
	xax = array([patrDistList,]*len(supportList)).transpose()
	yax = array(nr_clusters).transpose()
	plt.plot(xax,yax,marker='o')
	plt.xlabel('Patristic distance threshold')
	plt.ylabel('Number of clusters')
	plt.legend(["support=" + s + "%" for s in array(supportList).astype(str) ]) 
 
##visualization 3
if (figure_sizeof_clusters):
	p3 = plt.figure("Mean clusters sizes")
	xax = array([patrDistList,]*len(supportList)).transpose()
	yax = array(cluster_size_dist_mean).transpose()
	plt.plot(xax,yax,marker='o')
	plt.xlabel('Patristic distance threshold')
	plt.ylabel('Mean size of clusters')
	plt.legend(["support=" + s + "%" for s in array(supportList).astype(str) ])
	
if (figure_distances or figure_nr_clusters or figure_sizeof_clusters):
    plt.show()
