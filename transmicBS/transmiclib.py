# LIBRARY OF TRANSMICBS FUNCTIONS


from __future__ import with_statement
import math
import dendropy
from dendropy import treecalc
import sys 
import csv
import string

from numpy import *
from Bio import motifs
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import Alphabet
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from xlrd import open_workbook
from xlwt import Workbook
from xlutils.copy import copy

#Parser for the control file 
def parse_input_file(filename):
    try:
        fstream = open(filename, "r")
        rlines = fstream.readlines()
        fstring = "";
        for s in rlines:
            if s[0]!= '#': 			
                fstring=fstring+s; 					 
        fstream.close()
    except IOError:
        print 'Input error: the name of the control file must be "controlfile.tbsc" \n\n '    				
    paramlist = []              
    #1. tree file name     
    str1  = fstring    
    idx1 = str1.find('"',1)  
    idx2 = str1.find('"',idx1+1)      
    treefile_name = str1[idx1+1:idx2]
    paramlist.append(treefile_name)
    
    #2. clustering threshold
    idx3 = str1.find('"',idx2+1) 
    idx4 = str1.find('"',idx3+1)
    cutoff_mppd = float(str1[idx3+1:idx4])
    paramlist.append(cutoff_mppd)
      
    #3. bootstrap support threshold
    idx5 = str1.find('"',idx4+1) 
    idx6 = str1.find('"',idx5+1) 
    cutoff_bootstrap = float(str1[idx5+1:idx6])
    paramlist.append(cutoff_bootstrap) 
       
    #4. patient-wise clustering
    idx7 = str1.find('"',idx6+1) 
    idx8 = str1.find('"',idx7+1)
    if str1[idx7+1:idx8] == "true":
        paramlist.append(1)
    else:
        paramlist.append(0)
    
	#5. write out clusters to a tabular file?
    idx9 = str1.find('"',idx8+1) 
    idx10 = str1.find('"',idx9+1)
    if str1[idx9+1:idx10] == "true":
        paramlist.append(1)
    else:
        paramlist.append(0)
                 
    #6. name of .xls-file with sequence IDs/patient IDs and
    #potentially meta information for each sequence
    idx11 = str1.find('"',idx10+1) 
    idx12 = str1.find('"',idx11+1)
    sequence_table_name = str1[idx11+1:idx12]
    paramlist.append(sequence_table_name)
           				
    #7. set true if string format of Robert Koch Institute is required 
    idx13 = str1.find('"',idx12+1) 
    idx14 = str1.find('"',idx13+1)
    if str1[idx13+1:idx14] == "true":
		paramlist.append(1)
    else:
		paramlist.append(0)           				
           					
    #8. sequence id of the outgroup    
    idx15 = str1.find('"',idx14+1) 
    idx16 = str1.find('"',idx15+1)
    outgroupSeqId = str1[idx15+1:idx16]
    paramlist.append(outgroupSeqId)
                                                                       
    return paramlist

    
#Read in bootstrap supports for each node
def is_node_significant_raxml(node, bootstrap_cutoff):
    """Checks whether a given tree node has sufficient bootstrap support"""
    if node.is_leaf():
        node.support = 1
        return node.support      
 
    #recursive function call
    #support_of_downstream_subtree = []        
    for n in node.child_nodes():
        #support_of_downstream_subtree.append(is_node_significant_raxml(n, bootstrap_cutoff))
        is_node_significant_raxml(n, bootstrap_cutoff)
    #Algorithm:
    #compare the taxa (leaves) of the node to the list of clades from the consensus tree
    #and if the exact list of the taxa matches one of the clades from the consensus tree
    #then check if the corresponding support value is greater-equal than the bootstrap cutoff. 
    #If both queries are positive, recursively check if all these queries are also positive
    #for all child nodes of the node in question. If the latter is true then mark the node in
    #question as significant.
    tn_set=set(taxon_names(node))
    
    #the top-level node is always considered as significant
    if node.level() == 0:
        support_val = 1
    else:
        if(len(node.get_node_str()))==0:
			support_val=0
        else:                		
            support_val = float(node.get_node_str())/100
		
    if support_val >= float(bootstrap_cutoff)/100:
		bootstrap_support = True
    else:
		bootstrap_support = False
		   
    #node.support = support_val
    node.support = bootstrap_support
    node.support_val = support_val
    return node.support    

#This function computes the mean of all pairwise patristic distances of the taxa descending from a given node.
#If the mean of patristic distances of the node is less than cutoff and the node has sufficient bootstrap support,
#then the taxa descending of this node are identified as a transmission cluster. 
def get_clusters_meandist(node, cutoff):
    cluster_list = []
    if node.is_leaf():
        return cluster_list 
    nr_leafs = get_num_leaf_nodes(node, 0)    
    avg=get_average_length2(node, nr_leafs)
    sum_dist = sum(avg)
    
    #devide the sum of patristic distances by the possible number of pairs (binomial)
    mean_dist = sum_dist/nchoosek(nr_leafs,2)
    
    current_clade = taxon_names(node)
    bootstrap_accept = node.support
    
    #if node is root assign weight zero
    nodeweight = node.support_val
	    
    if mean_dist <= float(cutoff)/100 and bootstrap_accept:
        current_clade += [mean_dist, nodeweight]
        cluster_list += [current_clade]
    else:
        for n in node.child_nodes():
            cluster_list += [get_clusters_meandist(n, cutoff)]
    return cluster_list
    
#the same as  get_clusters_mediandist() but the distance
#is only computed between individuals (no distances between follow-up
#sequences of the same individual are computed)
def get_clusters_meandist_interPatient(node, cutoff, pdm, outgroupSeqId):
    cluster_list = []
    if node.is_leaf():
        return cluster_list 
           
    #compute sequenc-wise inter-patient distance
    mean_dist = get_average_length3(node, pdm, outgroupSeqId)
    
    current_clade = taxon_names(node)     
    bootstrap_accept = node.support
    
    #if node is root assign weight zero
    #if node.level()==0:
    nodeweight = node.support_val		
    #else:
	#	nodeweight = node.support_val
	    
    if mean_dist == -1:
      return cluster_list
    else:       
      if mean_dist <= float(cutoff)/100 and  bootstrap_accept:
          current_clade += [mean_dist, nodeweight]
          cluster_list += [current_clade]   
      else:
          for n in node.child_nodes():
              cluster_list += [get_clusters_meandist_interPatient(n, cutoff, pdm, outgroupSeqId)]
    return cluster_list	

#A function running recursively through the subtree descending from 'node' and calculating the length 
#of the branches for each node, and weighting the branch lengths according to the number of times 
#they appear in the set of pairwise patristic distances.
def get_average_length2(node, nrLeafsCluster):
  own_lengths=[]
  if node.is_leaf():
    return [node.edge.length*(nrLeafsCluster-1)] # nrLeafsCluster should always be >= 2, since this function is only called for inner nodes
  if node.edge.length is not None:
    own_lengths += [(nrLeafsCluster-get_num_leaf_nodes(node, 0))*get_num_leaf_nodes(node, 0)*node.edge.length]
  for n in node.child_nodes():
    own_lengths += get_average_length2(n, nrLeafsCluster)
  return own_lengths
  
#Same as get_average_length2, but computes only inter-individual distances
def get_average_length3(node, pdm, outgroupSeqId):
  patIdList=[] #list of patIDs (the column seSKNr in the RKI table)	
  leaflist = node.leaf_nodes()
  patIDdict = dict()
  for i in range(len(leaflist)):  
    try: 	  
      patID=leaflist[i].get_node_str().strip("'").split(' ')[1]
    except IndexError:		
      if leaflist[i].get_node_str()==outgroupSeqId:
	    continue
      else:
        raise NameError('Sequence ID does not have the proper format.') 		    	    
    if patID not in patIDdict:
      patIDdict[patID]=[]
    patIDdict[patID].append(leaflist[i])     
  keylist = patIDdict.keys()  
  minDistList=[]
  dictsz = len(keylist)  
  if dictsz<2:
    return -1
  else:     
    for i in range(dictsz):
      leafList1=patIDdict[keylist[i]]     		
      for j in range(i+1,dictsz):
        leafList2=patIDdict[keylist[j]]     		  		
        minDistList.append(minDistTwoPatients(leafList1, leafList2, pdm))
  return mean(minDistList)    

#A function running recursively through the tree and emitting the names of corresponding OTUs for each node.
#Note: The Dendropy library puts name strings of leafs in single quotes, if the name contains a hyphen.
#Therefore the function taxon_names also searches for single quotes as a part of the name string and extracts
#the actual name string.
def taxon_names(node):
  if node.is_leaf():
	namestr  = node.get_node_str()
	modstr =  namestr.split("'")
	if len(modstr) > 1:
		namestr = modstr[1]
	return [namestr.replace(' ', '_')]
  taxo = []
  for n in node.child_nodes():
    taxo += taxon_names(n)
  return taxo

#The minimal distance between any two sequences belonging to two individuals
def minDistTwoPatients(leafList1, leafList2, pdm):
  sz1=len(leafList1)
  sz2=len(leafList2)
  minDist=inf
  for i in range(sz1):
    taxon1 = leafList1[i].taxon  
    for j in range(sz2):
	  	taxon2 = leafList2[j].taxon
	  	minDist = min(minDist, pdm(taxon1, taxon2))
  return minDist	  	

#Auxiliary function      		 	  	     
def nchoosek(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)
 	
#Auxiliary function for unfolding nested lists with multiple levels, to lists of elementary lists
def flatten(S):
    if S == []:
        return S
    if not isinstance(S[0], list):
        return [S]    
    if isinstance(S[0], list):
        return flatten(S[0]) + flatten(S[1:])
    return S[:1] + flatten(S[1:])

#A function running recursively through the tree and calculating the number of leaves of a node 
def get_num_leaf_nodes(node, num):
  if node.is_leaf():
    return num+1
  for n in node.child_nodes():
    num = get_num_leaf_nodes(n, num)
  return num
      			                       
#This function writes cluster lists to a file.
#It can have either one or two list arguments. 
def write_clusters_to_file(cluster_list1, filename, nr_seqs, patrDist_cutoff,  bootstrap_cutoff):
    fobj = open(filename, "w") 
    
    #sorts by the clustering thresold
    cluster_list1.sort(key=lambda tup: tup[-2])
               
    len1 =len(cluster_list1)
    cc = 0
    for i in range(len1):
		cc= cc + len(cluster_list1[i][0:-2])
    print >> fobj, "This is an output file of TransmicBS.\n"
	
    if len1 > 0:
        print >> fobj, "Transmission clusters computed according to the distance  threshold of", patrDist_cutoff,"% and a bootstrap threshold of", bootstrap_cutoff,"%."
        if nr_seqs > 0:
			print >> fobj, "Number of sequences: ", nr_seqs
        print >> fobj, "Number of identified clusters: ", len1
        print >> fobj, "Average number of sequences in a cluster: ", "{0:.2f}".format(float(cc)/float(len1))
        if nr_seqs > 0:					
            print >> fobj, "Number of sequences not included in clusters: ", int(round(nr_seqs - len1* float(cc)/float(len1)))         
    else:
        print >> fobj, "No transmission clusters were found. Try different cut-off values."    
    print >> fobj, ""
    
    for i in range(len1):
        print >> fobj,  str(i+1)+"."
        print >> fobj, "Mean value of pairwise patristic distances (%):", "{0:.4f}".format(cluster_list1[i][-2]*100)
        print >> fobj, "Bootstrap support (%):", "{0:1.1f}".format(cluster_list1[i][-1]*100) 
        print >> fobj, cluster_list1[i][0:-2]
        print >> fobj, ""   
    fobj.close()
    return cluster_list1

#Auxiliary function
def RKI_format_seqID(idstr):
	str1 = idstr.split('_',1)[0].split('-')[0]
	str2 = "0"+idstr.split('_',1)[0].split('-')[1]
	seq_id = str1+"-"+str2 
	return(seq_id)

#Write clusters to excel             
def write_clusters_to_table(cluster_list1, RKIxls_infilename, xls_outfilename, taxon_labels, outgroupSeqId, enforce_RKI_format):
	book_orig = open_workbook(RKIxls_infilename)
	sheet_orig = book_orig.sheet_by_index(0)
	col_len = len(sheet_orig.col_values(0))
	row_len = len(sheet_orig.row_values(0))   

	#book_new = Workbook()
	#sheet_new = book_new.add_sheet('Clustering results')
	#sheet_new.write(0,0,'ClusterNr')   
	#sheet_new.write(0,1,'SeqId')
    
    #write the header of the new xls table
	with open(xls_outfilename, 'wb') as csvfile:
		wr = csv.writer(csvfile, delimiter='\t')
		wr.writerow(['ClusterNr'] + ['SeqId'] + [cell.value.encode('utf8') if isinstance(cell.value, unicode) else cell.value for cell in sheet_orig.row(0)] + ['Mppd'] + ['Bootstrap'])
		#wr.writerow(['ClusterNr'] + ['SeqId'])
		#for i in range(row_len):
		#	sheet_new.write(0,i+2,sheet_orig.cell(0,i).value)
		
		scount_list_orig = sheet_orig.col_values(0)[1:]
		len1 =len(cluster_list1)
		seq_idx=0
		for i in range(len1):
			for j in range(len(cluster_list1[i])-2):
				seq_idx = seq_idx + 1 
				idstr = cluster_list1[i][j]
				#mppd and bootstrap value per cluster
				mppd = cluster_list1[i][-2]*100
				bootstrap = cluster_list1[i][-1]*100
				if not outgroupSeqId == idstr: 	# make sure it is nor the the outgroup sequence					
					str1 = idstr.split('_',1)[0].split('-')[0]
					str2 = "0"+idstr.split('_',1)[0].split('-')[1]
					seq_id = str1+"-"+str2
					if enforce_RKI_format:        
						seq_id  = RKI_format_seqID(idstr)                        																							
					else:
						seq_id  = idstr.split('_',1)[0]							
					idx = [k for k,x in enumerate(scount_list_orig) if x == seq_id]
					if len(idx) != 1: 
						print idx
						print seq_id
						print 'Error: Seq.Id. found multiple times or not found'    
					wr.writerow([i+1] + [idstr] + [cell.value.encode('utf8') if isinstance(cell.value, unicode) else cell.value for cell in sheet_orig.row(idx[0]+1)]+ [mppd] + [bootstrap])
					#wr.writerow([i+1] + [idstr])
					#sheet_new.write(seq_idx,0,i+1) # index of the cluster
					#sheet_new.write(seq_idx,1,idstr)
					#for ii in range(row_len):
					#	sheet_new.write(seq_idx,ii+2,sheet_orig.cell(idx[0]+1,ii).value)
				taxon_labels.remove(idstr)
                	            
	 # print out sequences not included in clusters 
		i=0
		for j in range(len(taxon_labels)):
			idstr = taxon_labels[j]
			if not outgroupSeqId in idstr:
				i=i+1
				str1 = idstr.split('_',1)[0].split('-')[0]    
				str2 = "0"+idstr.split('_',1)[0].split('-')[1]
				seq_id = str1+"-"+str2
				idx = [k for k,x in enumerate(scount_list_orig) if x == seq_id]
				if len(idx) != 1:
					print 'Error: Seq.Id. found multiple times or not found'
				wr.writerow([0] + [idstr]+[cell.value.encode('utf8') if isinstance(cell.value, unicode) else cell.value for cell in sheet_orig.row(idx[0]+1)]+ [''] + [''])
   
