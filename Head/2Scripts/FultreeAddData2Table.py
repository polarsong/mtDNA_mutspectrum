#!/Users/xtinaushakova/opt/anaconda3/bin/python

"""

REQUIREMENTS: This script requires biopython to run.

This simple script appends data on aminoacid substitutions to fulltree.csv table. 
It determines ancestral/derived aa in accordance to mitochondrial translation table.
Adds the following columns: 
	"pos_in_codon" - position of substitution in codon (1, 2, 3)
	"synonymous" - is the mutation synonymous or not 
	"ancestral_aa", "derived_aa" - self explanatory 
	"note" - any pecularities that should be noted

If the substitution ocurred in a non-coding sequence the script appends NAs to all columns

"""


import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# Open our table. This is the way around relative paths in Python
dirname = os.path.dirname(os.path.abspath(__file__))
sourceTablePath = os.path.join(dirname, "../../Body/3Results/fulltree.csv")
sourceTable = open(sourceTablePath)

# Aminoacid reference
aminoacids = {'C' : 'Cys', 'D' : 'Asp', 'S' : 'Ser', 'Q' : 'Gln', 'K' : 'Lys',
     		  'I' : 'Ile', 'P' : 'Pro', 'T' : 'Thr', 'F' : 'Phe', 'N' : 'Asn', 
              'G' : 'Gly', 'H' : 'His', 'L' : 'Leu', 'R' : 'Arg', 'W' : 'Trp', 
              'A' : 'Ala', 'V' : 'Val', 'E' : 'Glu', 'Y' : 'Tyr', 'M' : 'Met',
              '*' : 'Stop', 'X' : 'Ambiguous', 'J' : 'Leu/Ile', 'B' : 'Asn/Asp',
              'Z' : 'Gln/Glu'}

# Open new table fulltreeCodons for writing
newTablePath = os.path.join(dirname, "../../Body/3Results/fulltreeCodons.csv")
newTable = open(newTablePath, 'wt')

# Read header in old file append and store in new
header = sourceTable.readline().strip('\n')
newHeader = ";".join([header, "pos_in_codon", "synonymous", "ancestral_aa", "derived_aa", "note"]) 
newTable.write(newHeader + "\n")

for line in sourceTable:
	line = line.strip('\n').split(';')
	name = line[-1]
	# Only take mRNA
	if 'mRNA' in name:
		posTable = int(line[3]) # get mutation pos in ref from table
		
		# DISREPANCY CHECK DONE
		# if posTable - startRecord < 0:
		#	 print("No good")

		ancestralSeq = line[4]
		derivedSeq = line[5]
		# determine position in codon
		posInCodon = 3 - (posTable % 3)
		# Get codon from table
		codonStart = 3 - posInCodon
		codonEnd = codonStart + 3
		ancestralCodon = ancestralSeq[codonStart:codonEnd]
		derivedCodon = derivedSeq[codonStart:codonEnd]

		# If gaps in either sequences then ignore
		if '-' in ancestralCodon or '-' in derivedCodon:
			line.extend([posInCodon, "NA", "NA", "NA", "gaps"])
			next
		else:
			ancestralCodon = Seq(ancestralSeq[codonStart:codonEnd], generic_dna)
			derivedCodon = Seq(derivedSeq[codonStart:codonEnd], generic_dna)
			ancestralAa = aminoacids[ancestralCodon.translate(table=2)]
			derivedAa = aminoacids[derivedCodon.translate(table=2)]
			synonymous = "non-synonymous"
			#synonymous or not
			if ancestralAa == derivedAa:
				synonymous = "synonymous"
			note = "normal"
			line.extend([posInCodon, synonymous, ancestralAa, derivedAa, note])
	else:
		line.extend(["NA", "NA", "NA", "NA", "NA"])
	print(line)
	line = ";".join(str(element) for element in line)
	newTable.write(line + "\n")

sourceTable.close()
newTable.close()
