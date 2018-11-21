#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
	This script parses .gb source file pulls out complete genomes and stashes them into one file 		genomes.txt

"""

from Bio import SeqIO
from collections import Counter

repeating = ['Bradypus_variegatus', 'Rattus_norvegicus', 'Lateolabrax_maculatus', 'Coilia_nasus']

names = {
			'NADH dehydrogenase subunit 6' : 'ND6',
			'ATP synthase F0 subunit 6' : 'ATP6',
			'NADH dehydrogenase subunit 2' : 'ND2',
			'cytochrome c oxidase subunit II' : 'COX2',
			'cytochrome b' : 'CytB',
			'NADH dehydrogenase subunit 4L' : 'ND4L',
			'cytochrome c oxidase subunit III' : 'COX3',
			'NADH dehydrogenase subunit 1' : 'ND1',
			'NADH dehydrogenase subunit 3' : 'ND3',
			'NADH dehydrogenase subunit 4' : 'ND4',
			'ATP synthase F0 subunit 8' : 'ATP8',
			'cytochrome c oxidase subunit I' : 'COX1',
			'NADH dehydrogenase subunit 5' : 'ND5',
			'NADH denydrogenase subunit 1' : 'ND1',
			'cytochrome c oxidase subunit 1' : 'COX1',
			'cytochrome c oxidase subunit 2' : 'COX2',
			'cytochrome c oxidase subunit 3' : 'COX3'
		}

infile = open("/home/kristina/WB/tRNA/1_raw_data/0_source_file/source.gb")

outAla = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Ala.tRNA", "wt")
outArg = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Arg.tRNA", "wt")
outAsn = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Asn.tRNA", "wt")
outAsp = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Asp.tRNA", "wt")
outCys = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Cys.tRNA", "wt")
outGln = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Gln.tRNA", "wt")
outGlu = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Glu.tRNA", "wt")
outGly = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Gly.tRNA", "wt")
outHis = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/His.tRNA", "wt")
outIle = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Ile.tRNA", "wt")
outLeuCUN = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/LeuCUN.tRNA", "wt")
outLeuUUR = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/LeuUUR.tRNA", "wt")
outLeuUn = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/LeuUnknown.tRNA", "wt") 
outLys = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Lys.tRNA", "wt")
outMet = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Met.tRNA", "wt")
outPhe = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Phe.tRNA", "wt")
outPro = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Pro.tRNA", "wt")
outSerAGY = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/SerAGY.tRNA", "wt")
outSerUCN = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/SerUCN.tRNA", "wt")
outSerUn = open("//home/kristina/WB/tRNA/2_derived_data/2_trna/SerUnknown.tRNA", "wt") 
outThr = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Thr.tRNA", "wt")
outTrp = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Trp.tRNA", "wt")
outTyr = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Tyr.tRNA", "wt")
outVal = open("/home/kristina/WB/tRNA/2_derived_data/2_trna/Val.tRNA", "wt")

outfile_genomes = open("/home/kristina/WB/tRNA/2_derived_data/0_genomes/genomes.fasta", "wt")

with open("/home/kristina/WB/tRNA/2_derived_data/outliars.gb", "wt") as output_handle:
	for rec in SeqIO.parse(infile, 'genbank'):
		species = "_".join(rec.annotations['organism'].split())
		if species in repeating:
			species = species + "_" + rec.name
		
		print(species)
		taxonomy = "_".join(rec.annotations["taxonomy"])
		count = Counter(rec.seq)
		(A, C, G, T) = (count["A"], count["C"], count["G"], count["T"])
		X = str(len(rec.seq) - (A + C + G + T))
		outfile_genomes.write(">" + rec.name + "\t" + species + "\n" + str(rec.seq) + "\n")
		
		if len([i for i in rec.features if i.type.upper() == "TRNA"]) > 22:
			SeqIO.write(rec, output_handle, "genbank")
			next
		
		outCDS = open("/home/kristina/WB/tRNA/2_derived_data/1_CDS/%s.fasta" % (species), "wt") 
		for feature in rec.features:
			if feature.type == "CDS":
				nucleotide = feature.extract(rec.seq)
				if len(nucleotide) % 3 != 0:
					nucleotide = nucleotide[0: -(len(nucleotide) % 3)]
				translation = feature.qualifiers['translation']	
				name = names[feature.qualifiers['product'][0]]
				outCDS.write(">%s\t%s\tnucleotide\n%s\n" % (name, species, nucleotide))
				outCDS.write(">%s\t%s\ttranslation\n%s\n" % (name, species, translation[0]))
		

			if feature.type == 'tRNA':
				if feature.qualifiers['product'] == ['tRNA-Ala']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outAla.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Arg']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outArg.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Asn']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outAsn.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Asp']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outAsp.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Cys']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outCys.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Gln']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outGln.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Glu']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outGlu.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Gly']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outGly.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-His']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outHis.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Ile']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outIle.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				# 3 cases for Leu and Ser
				if feature.qualifiers['product'] == ['tRNA-Leu']:
					if 'codon_recognized' in feature.qualifiers and feature.qualifiers['codon_recognized'] == 'CUN':
						seq = str(feature.extract(rec.seq)).upper()
						seq = seq.replace("T", "U")
						outLeuCUN.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
					elif 'codon_recognized' in feature.qualifiers and feature.qualifiers['codon_recognized'] == 'UUR':
						seq = str(feature.extract(rec.seq)).upper()
						seq = seq.replace("T", "U")
						outLeuUUR.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
					else:
						seq = str(feature.extract(rec.seq)).upper()
						seq = seq.replace("T", "U")
						outLeuUn.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Lys']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outLys.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Met']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outMet.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Phe']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outPhe.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Pro']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outPro.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Ser']:
					if 'codon_recognized' in feature.qualifiers and feature.qualifiers['codon_recognized'] == 'AGY':
						seq = str(feature.extract(rec.seq)).upper()
						seq = seq.replace("T", "U")
						outSerAGY.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
					elif 'codon_recognized' in feature.qualifiers and feature.qualifiers['codon_recognized'] == 'UCN':
						seq = str(feature.extract(rec.seq)).upper()
						seq = seq.replace("T", "U")
						outSerUCN.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
					else:
						seq = str(feature.extract(rec.seq)).upper()
						seq = seq.replace("T", "U")
						outSerUn.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Thr']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outThr.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Trp']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outTrp.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Tyr']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outTyr.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
				if feature.qualifiers['product'] == ['tRNA-Val']:
					seq = str(feature.extract(rec.seq)).upper()
					seq = seq.replace("T", "U")
					outVal.write(">%s\t%s\n%s\n" % (species, feature.qualifiers['product'][0], seq))
		outCDS.close()
outAla.close()
outArg.close()
outAsn.close()
outAsp.close()
outCys.close()
outGln.close()
outGlu.close()
outGly.close()
outHis.close()
outIle.close()
outLeuCUN.close()
outLeuUUR.close()
outLeuUn.close()
outLys.close()
outMet.close()
outPhe.close()
outPro.close()
outSerAGY.close()
outSerUCN.close()
outSerUn.close()
outThr.close()
outTrp.close()
outTyr.close()
outVal.close()

outfile_genomes.close()
