import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd

# PATH
dirname = os.path.dirname(os.path.abspath(__file__))

# SOURCE
sourceTablePath = os.path.join(dirname, "../../Body/1Raw/Wei2019HumanDuos.txt")
sourceTable = open(sourceTablePath)

# HARD REFERENCE
aminoacids = {'C' : 'Cys', 'D' : 'Asp', 'S' : 'Ser', 'Q' : 'Gln', 'K' : 'Lys',
     		  'I' : 'Ile', 'P' : 'Pro', 'T' : 'Thr', 'F' : 'Phe', 'N' : 'Asn', 
              'G' : 'Gly', 'H' : 'His', 'L' : 'Leu', 'R' : 'Arg', 'W' : 'Trp', 
              'A' : 'Ala', 'V' : 'Val', 'E' : 'Glu', 'Y' : 'Tyr', 'M' : 'Met',
              '*' : 'Stop', 'X' : 'Ambiguous', 'J' : 'Leu/Ile', 'B' : 'Asn/Asp',
              'Z' : 'Gln/Glu'}

genes = {'ND1' : "ND1", 'ND2': "ND2", "CO1" : 'COX1', 
		 "CO2" : 'COX2', 'ATP8': "ATP8", 'ATP6': "ATP6", 
		 "CO3" : 'COX3', 'ND3': "ND3", "ND4L" : 'ND4L', 
		 'ND4': "ND4", 'ND5': "ND5", 'ND6': "ND6", "CYB" : 'CYTB'}

position = {1 : 1, 2 : 2, 0 : 3}

# REFERENCE
# Parse genome reference and get gene refs
refGenes = {"geneName" : [],
			"startEnd" : [],
			"sequence" : []}

refGenbank = os.path.join(dirname, "../../Body/3Results/NC_012920.1.gb")
refGenbank = open(refGenbank)

for rec in SeqIO.parse(refGenbank, "genbank"):
    if rec.features:
        for feature in rec.features:
            if feature.type == "CDS":
            	# Extract gene name and append to dict
            	featureName = feature.qualifiers['gene'][0]
            	refGenes["geneName"].append(featureName)
            	# Extract gene location as range and append to dict
            	featureLocation = [int(feature.location.start), int(feature.location.end)]
            	refGenes["startEnd"].append(featureLocation)
            	# Extract gene sequence string and append to dict
            	featureSequence = str(feature.location.extract(rec).seq)
            	refGenes["sequence"].append(featureSequence)

# Turn into data frame
refGenes = pd.DataFrame.from_dict(refGenes)

# WRITING
newTablePath = os.path.join(dirname, "../../Body/3Results/Wei2019HumanDuosDerived.txt")
newTable = open(newTablePath, 'wt')

# Read header in old file append and store in new
header = sourceTable.readline().strip('\n')
header = sourceTable.readline().strip('\n').split('\t')
header.extend(i for i in ["fromTo", "heavyNotation", "codingOrNot", "SynonymousOrNot", "posInCodon","fromToAA", "note"])
newHeader = "\t".join(header)

newTable.write(newHeader + "\n")

# Line for line in Wei
for line in sourceTable:
	line = line.strip('\n').split('\t')

	# Extract mutation from table
	mutation = line[2]
	mutationType = mutation[0] + mutation[-1]
	mutationRef = int(mutation[1:-1])
	
	# Flip to heavy notation
	mutationHeavy = Seq(mutationType, generic_dna).complement()

	# Add > to mutation string
	mutationType = mutationType[0] + ">" + mutationType[-1]
	mutationHeavy = mutationHeavy[0] + ">" + mutationHeavy[-1]

	# Determine if coding
	gene = line[-2]
	codingOrNot = "coding"

	# if non-coding we stop here
	if gene in ["T", "DLOOP", "RNR1", "RNR2", "OLR", "OHR"]:
		codingOrNot = "non-coding"
		line.extend([mutationType, mutationHeavy, "NA", "NA", "NA", "NA", "non-coding"])
		next

	# if coding we go further
	else:
		# Get position in gene
		geneName = genes[gene]
		geneLocationList = list(refGenes[refGenes.values  == geneName]["startEnd"])[0]
		geneLocation = range(geneLocationList[0], geneLocationList[1])
		geneStart = geneLocationList[0]
	
		posInGene = mutationRef - geneStart

		# CHECK IF NEGATIVE
		if posInGene < 0:
			line.extend([codingOrNot, "NA", "NA", "NA", "out of gene boundaries"])
		else:
			# get position in codon
			posInCodon = position[posInGene % 3]

			# get codon
			codonStart = posInGene - posInCodon
			codonEnd = codonStart + 3
			# Get codon from ref sequence
			codonRef = list(refGenes[refGenes.values  == geneName]["sequence"])[0][codonStart:codonEnd]
			codonRefList = list(codonRef)
			# Get previous codon
			codonAnc = codonRefList.copy()
			codonAnc[posInCodon - 1] = mutationType[0]
			codonAnc = Seq("".join(codonAnc), generic_dna)
			# And mutated
			codonDer = codonRefList.copy()
			codonDer[posInCodon - 1] = mutationType[-1]
			codonDer = Seq("".join(codonDer), generic_dna)
			# get aminoacids
			ancestralAa = aminoacids[codonAnc.translate(table=2)]
			derivedAa = aminoacids[codonDer.translate(table=2)]

			# synonymous or not
			synonymousOrNot = "synonymous"
			if ancestralAa != derivedAa:
				synonymousOrNot = "non-synonymous"
			fromToAA = ancestralAa + ">" + derivedAa
			print(fromToAA, synonymousOrNot, posInCodon)

			# Append all results to line and write to file
			line.extend([mutationType, mutationHeavy, codingOrNot, synonymousOrNot, str(posInCodon), fromToAA, "normal"])
			newTable.write("\t".join(line) + "\n")

sourceTable.close()
newTable.close()