# LIBRARY OF TRANSMIC FUNCTIONS

from __future__ import with_statement
import math
import dendropy
from dendropy import calculate
import sys 
import csv
import string
import xlwt
import re
import matplotlib.pyplot as plt
 

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
from xlutils.copy import copy
from numpy import *

#Parser for the control file 
def parse_input_file(filename):
    try:
        fstream = open(filename, "r")
        rlines = fstream.readlines()
        fstring = "";
        for s in rlines:
            if s[0]!= '#' and s: 			
                fstring=fstring+s; 					 
        fstream.close()
    except IOError:
        print 'Input error: the name of the control file must be "controlfile.tbsc" \n\n '            
    rowlist = fstring.split("\n")
    paramlist=[]
    for item in rowlist:
		if not(item == '' or re.match('\s+$', item)):
			paramlist.append(parser_aux(item))
    return paramlist
        				    
def parser_aux(parstring):
	idx1 = parstring.find('"',0)
	idx2 = parstring.find('"',idx1+1)
	return parstring[idx1+1:idx2]
	
#Read in bootstrap supports for each node
def assess_node_support(node, cutoff, typeOfSupport):
    """Checks whether a given tree node has sufficient bootstrap support"""
    if node.is_leaf():
        node.support = 1
        return node.support      
 
    #recursive function call
    #support_of_downstream_subtree = []        
    for n in node.child_nodes():
        #support_of_downstream_subtree.append(is_node_significant_raxml(n, bootstrap_cutoff))
        assess_node_support(n, cutoff, typeOfSupport)
    
    #the top-level node (e.g. root) is always considered as significant
    if node.level() == 0:
        support_val = 1
    else:
        if typeOfSupport=="bootstrap":
            #print node
            #print node.child_nodes()
            #print node.is_leaf()    #---> include try-statement
            try:                    			
                label_len = len(node.label)
            except TypeError:
				label_len = 0 
            if label_len == 0:
            #if node.label==None:
			    support_val=0
            else: 			    			                		
                support_val = float(node.label)/100
        elif typeOfSupport=="posterior":
            try:
                support_val = float(node.annotations.get_value("posterior"))
            except TypeError:
				support_val =  0                				                
        else: raise NameError("type of support should be either bootstrap or posterior.")           
                        			                    			                                		
    if support_val >= float(cutoff)/100:
		node_support = True
    else:
		node_support = False
	
    #node.support = support_val
    node.support = node_support
    node.support_val = support_val
    return node.support   

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
    tn_set=set(taxon_names(node))
    
    #the top-level node is always considered as significant
    if node.level() == 0:
        support_val = 1
    else:
        if(len(node.label))==0:
			support_val=0
        else:                		
            support_val = float(node.label)/100
		
    if support_val >= float(bootstrap_cutoff)/100:
		bootstrap_support = True
    else:
		bootstrap_support = False
		   
    #node.support = support_val
    node.support = bootstrap_support
    node.support_val = support_val
    return node.support    

#sequence-wise clustering
def get_clusters_meandist(node, dist_cutoff):
    cluster_list = []
    if node.is_leaf():
        return cluster_list 
    nr_leafs = get_num_leaf_nodes(node, 0)    
    avg=get_average_length2(node, nr_leafs)
    sum_dist = sum(avg)
    
    #devide the sum of patristic distances by the possible number of pairs (binomial)
    mean_dist = sum_dist/nchoosek(nr_leafs,2)   
    current_clade = taxon_names(node)
 	    
    if mean_dist <= float(dist_cutoff) and node.support: #node.posterior >= posterior_cutoff:
        current_clade += [mean_dist, node.support_val]
        cluster_list += [current_clade]
    else:
        for n in node.child_nodes():
            cluster_list += [get_clusters_meandist(n, dist_cutoff)]
    return cluster_list    


#the same as  get_clusters_mediandist() but the distance
#is only computed between individuals (no distances between follow-up
#sequences of the same individual are computed)
def get_clusters_meandist_interPatient(node, cutoff, outgroupSeqId, pdm):
    cluster_list = []
    if node.is_leaf():
        return cluster_list 
           
    #compute sequenc-wise inter-patient distance
    mean_dist = get_average_length3(node, outgroupSeqId, pdm)
    
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
      if mean_dist <= float(cutoff) and  bootstrap_accept:
          current_clade += [mean_dist, nodeweight]
          cluster_list += [current_clade]   
      else:
          for n in node.child_nodes():
              cluster_list += [get_clusters_meandist_interPatient(n, cutoff, outgroupSeqId, pdm)]
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
def get_average_length3(node, outgroupSeqId, pdm):
  patIdList=[] #list of patIDs 
  leaflist = node.leaf_nodes()
  patIDdict = dict()
  for i in range(len(leaflist)):  
      #patID=leaflist[i].get_node_str().strip("'").split(' ')[1] #V3.12
    try: 	  
      patID=leaflist[i].taxon.label.strip("'").split(' ')[1]             
    except IndexError:		
      #if leaflist[i].get_node_str()==outgroupSeqId:  #V3.12
      if leaflist[i].taxon.label==outgroupSeqId:		  
	    continue
      else:
        raise NameError('Sequence ID does not have the proper format.') 		    	    
    if patID not in patIDdict:
      patIDdict[patID]=[]
    patIDdict[patID].append(leaflist[i])     
  keylist = patIDdict.keys()  
  minDistList=[]
  dictsz = len(keylist)  
  if dictsz<2: #at least two patients required for distance computation
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
	namestr  = node.taxon.label
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
        print >> fobj, "Transmission clusters computed according to the distance  threshold of", patrDist_cutoff,"% and a support threshold of", bootstrap_cutoff,"%."
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
        print >> fobj, "Node support (%):", "{0:1.1f}".format(cluster_list1[i][-1]*100) 
        print >> fobj, cluster_list1[i][0:-2]
        print >> fobj, ""   
    fobj.close()

#Auxiliary function
def RKI_format_seqID(idstr):
	str1 = idstr.split('_',1)[0].split('-')[0]
	str2 = "0"+idstr.split('_',1)[0].split('-')[1]
	seq_id = str1+"-"+str2 
	return(seq_id)

#Write clusters to excel             
def write_clusters_to_table(cluster_list1, RKIxls_infilename, xls_outfilename, taxon_labels, outgroupSeqId, skipColumns):
	book_orig = open_workbook(RKIxls_infilename)
	sheet_orig = book_orig.sheet_by_index(0)
	col_len = len(sheet_orig.col_values(0))
	row_len = len(sheet_orig.row_values(0))   

    #write the header of the new xls table
	with open(xls_outfilename, 'wb') as csvfile:
		wr = csv.writer(csvfile, delimiter='\t')
		wr.writerow(['ClusterNr'] + ['SeqId'] + [cell.value.encode('utf8') if isinstance(cell.value, unicode) else cell.value for cell in sheet_orig.row(0)[skipColumns:]] + ['Mppd'] + ['Support, %'])		
		scount_list_orig = sheet_orig.col_values(0)[1:]
		len1 =len(cluster_list1)
		seq_idx=0
		for i in range(len1):
			for j in range(len(cluster_list1[i])-2):
				seq_idx = seq_idx + 1 
				idstr = cluster_list1[i][j]
				#mppd and bootstrap value per cluster
				mppd = cluster_list1[i][-2]
				bootstrap = cluster_list1[i][-1]*100
				if not outgroupSeqId == idstr: 	# make sure it is not the the outgroup sequence
					seq_id = idstr
					#seq_id = __extractScountFromIdStr__(idstr)    #---------------------------- specific for the TC-prject pouranyousef et al. (required for the RKI excel files)   
					#if enforce_RKI_format:        
					#	seq_id  = RKI_format_seqID(idstr)                        																							
					#else:
					#	seq_id  = idstr.split('_',1)[0]							
					idx = [k for k,x in enumerate(scount_list_orig) if x == seq_id]
					if len(idx) != 1: 
						print idx
						print seq_id
						print 'Error: Seq.Id. found multiple times or not found'
					#print i+1
					#print idstr   
					#print idx
					#print sheet_orig.row(idx[0]+1)
					#print sheet_orig.row(idx[0]+1)[skipColumns:]
					wr.writerow([i+1] + [idstr] + [cell.value.encode('utf8') if isinstance(cell.value, unicode) else cell.value for cell in sheet_orig.row(idx[0]+1)[skipColumns:] ]+ [mppd] + [bootstrap])
					taxon_labels.remove(idstr)
                	            
	 # print out sequences not included in clusters 
		for j in range(len(taxon_labels)):			
			idstr = taxon_labels[j]
			if not outgroupSeqId == idstr:				
				seq_id = idstr
				#seq_id = __extractScountFromIdStr__(idstr) #---------------------------- specific for the TC-prject pouranyousef et al. (required for the RKI excel files)  
				idx = [k for k,x in enumerate(scount_list_orig) if x == seq_id]
				if len(idx) != 1:
					print 'Error: Seq.Id. found multiple times or not found'
				rowstr = [cell.value.encode('utf8') if isinstance(cell.value, unicode) else cell.value for cell in sheet_orig.row(idx[0]+1)[skipColumns:] ]
				wr.writerow([0] + [idstr] + rowstr + [''] + [''])				
   
def __extractScountFromIdStr__(idstr):
    str1 = idstr.split('_',1)[0].split('-')[0]
    str2 = "0"+idstr.split('_',1)[0].split('-')[1]
    seq_id = str1+"-"+str2
    return seq_id
    
 











