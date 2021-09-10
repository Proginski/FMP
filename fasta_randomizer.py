# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 15:28:08 2021

@author: Snoopybreton
"""

def randomize(fasta_path):
    
    from Bio import SeqIO
    import random
    import sys
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    
    entries = SeqIO.parse(open(fasta_path),'fasta')
    output = []
    
    i = 0
    
    with open(fasta_path[0:-7]+"_randomized.nfasta",
              'a') as f_out:
    
        # Clean up the file
        f_out.truncate(0)
        
        for entrie in entries :
            
            # Randomize the seq
            seq = ''.join(random.sample(str(entrie.seq), len(str(entrie.seq))))
            entrie.description += "_randomized"
            
            # write new fasta file
            r=SeqIO.write(SeqRecord(Seq(seq),id=entrie.id,description=entrie.description), f_out, 'fasta')
            
            i += 1
            
            sys.stdout.write("\r"+str(i))
            sys.stdout.flush()
        
            if r!=1: print('Error while writing sequence:  ' + entrie.id)
            
            output.append(entrie)

    print("\n")
    
    return(output)
    

path_start = "/home1/paul.roginski/WD"
# path_start = "C:/Users/Snoopybreton/Desktop/BIM"

# ages = range(0,11)
# for age in ages : 
#     randomize(path_start + "/Levure/Abondance/Scer_CDS_age_"+str(age)+".nfasta")

# randomize(path_start + "/Levure/Abondance/Scer_CDS.nfasta")

# import os
# dir_path = "/home1/paul.roginski/WD/Archaea_Data/codingNfastaNP/"
# directory = os.listdir(dir_path)

# for file in directory :
#     randomize(dir_path+file)

randomize(path_start+ "/Levure/Abondance/Scer_CDS.nfasta")