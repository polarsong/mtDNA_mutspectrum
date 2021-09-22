#coding=utf8

import time
import timeit
import os
import sqlite3

from random import choice
import random





import numpy as np
import matplotlib.pyplot as plt








def pern(seq,n,s):
	newseq = ""
	for i in range(len(seq)):

		if i == n:
			newseq += s
		else:
			newseq += seq[i]			


	return newseq



def dia(fi,se):
	if fi == se:
		dia = "Syn"
	else:
		dia = "NonSyn"

	return dia



def synornot(dna):

	standart_map = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
	       "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
	       "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
	       "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
	       "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
	       "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
	       "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
	       "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
	       "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
	       "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
	       "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
	       "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
	       "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
	       "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
	       "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
	       "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}


	mito_map = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
	       "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
	       "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
	       "UGU":"C", "UGC":"C", "UGA":"W", "UGG":"W",
	       "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
	       "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
	       "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
	       "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
	       "AUU":"I", "AUC":"I", "AUA":"M", "AUG":"M",
	       "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
	       "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
	       "AGU":"S", "AGC":"S", "AGA":"*", "AGG":"*",
	       "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
	       "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
	       "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
	       "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}


	mito_map_dna = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
	       "TCT":"S", "TCC":"s", "TCA":"S", "TCG":"S",
	       "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
	       "TGT":"C", "TGC":"C", "TGA":"W", "TGG":"W",
	       "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
	       "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
	       "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
	       "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
	       "ATT":"I", "ATC":"I", "ATA":"M", "ATG":"M",
	       "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
	       "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
	       "AGT":"S", "AGC":"S", "AGA":"*", "AGG":"*",
	       "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
	       "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
	       "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
	       "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}


	#start = DNA.find('AUG')
	# start = 0
	# if start!= -1:
	#     while start+2 < len(DNA):
	#         codon = dna[start:start+3]
	#         #if codon == "UAG": break;
	#         print(map[codon])
	#         start+=3

	#print(dna)

	if len(dna) < 3:
		return "-"
	else:
		return mito_map_dna[dna]













cytba = ""


# 1045
for i in range(1045):

	ac = 0
	tc = 0
	gc = 0
	cc = 0

	allo = open("engraulis_Data/cytb_A.phy",'r')

	arr = []

	read = allo.readline().strip()

	while True:

		read = allo.readline().strip()

		if read == "":
			break

		row = read.split("  ")

		#print(row[1])

		if row[1][i] == "A":
			ac += 1
		elif row[1][i] == "T":
			tc += 1
		elif row[1][i] == "G":
			gc += 1
		elif row[1][i] == "C":
			cc += 1


	if ac >= tc and ac >= gc and ac >= cc:
		cytba = cytba + "A"
	elif tc >= ac and tc >= gc and tc >= cc:
		cytba = cytba + "T"
	elif gc >= ac and gc >= tc and gc >= cc:
		cytba = cytba + "G"
	elif cc >= ac and cc >= tc and cc >= gc:
		cytba = cytba + "C"


gens=open("cytba.txt","w")
gens.write("%s\n" % (cytba))




cytbb = ""


# 1045
for i in range(1045):

	ac = 0
	tc = 0
	gc = 0
	cc = 0

	allo = open("engraulis_Data/cytb_B.phy",'r')

	arr = []

	read = allo.readline().strip()

	while True:

		read = allo.readline().strip()

		if read == "":
			break

		row = read.split("  ")

		#print(row[1])

		if row[1][i] == "A":
			ac += 1
		elif row[1][i] == "T":
			tc += 1
		elif row[1][i] == "G":
			gc += 1
		elif row[1][i] == "C":
			cc += 1


	if ac >= tc and ac >= gc and ac >= cc:
		cytbb = cytbb + "A"
	elif tc >= ac and tc >= gc and tc >= cc:
		cytbb = cytbb + "T"
	elif gc >= ac and gc >= tc and gc >= cc:
		cytbb = cytbb + "G"
	elif cc >= ac and cc >= tc and cc >= gc:
		cytbb = cytbb + "C"


gens=open("cytbb.txt","w")
gens.write("%s\n" % (cytbb))










gens=open("cytba.csv","w")

for i in range(1045):




	a2t = 0
	a2g = 0
	a2c = 0
	t2a = 0
	t2g = 0
	t2c = 0
	g2a = 0
	g2t = 0
	g2c = 0
	c2a = 0
	c2t = 0
	c2g = 0

	allo = open("engraulis_Data/cytb_A.phy",'r')

	count = 0

	read = allo.readline().strip()

	while True:

		read = allo.readline().strip()

		if read == "":
			break

		row = read.split("  ")

		if row[1][i] != cytba[i]:

           

			codnum = (i) // 3 
			posincod = (i) % 3

			cs = codnum * 3


			#print(row[1][i-1:i+1],cytba[i-1:i+1],row[1][cs:cs+3],cytba[cs:cs+3])
			
			cytcod = synornot(cytba[cs:cs+3])



			if cytba[i] == "A" and row[1][i] == "T":
				a2t += 1
			elif cytba[i] == "A" and row[1][i] == "G":
				a2g += 1
			elif cytba[i] == "A" and row[1][i] == "C":
				a2c += 1
			elif cytba[i] == "T" and row[1][i] == "A":
				t2a += 1
			elif cytba[i] == "T" and row[1][i] == "G":
				t2g += 1
			elif cytba[i] == "T" and row[1][i] == "C":
				t2c += 1
			elif cytba[i] == "G" and row[1][i] == "A":
				g2a += 1
			elif cytba[i] == "G" and row[1][i] == "T":
				g2t += 1
			elif cytba[i] == "G" and row[1][i] == "C":
				g2c += 1
			elif cytba[i] == "C" and row[1][i] == "A":
				c2a += 1
			elif cytba[i] == "C" and row[1][i] == "T":
				c2t += 1
			elif cytba[i] == "C" and row[1][i] == "G":
				c2g += 1

	if a2t != 0:
		rowcod = pern(cytba[cs:cs+3],posincod,"T")
		res = synornot(rowcod)
		gens.write("%s;A>T;%s;%s;%s;%s;%s;%s\n" % (i, a2t, cytba[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if a2g != 0:
		rowcod = pern(cytba[cs:cs+3],posincod,"G")
		res = synornot(rowcod)
		gens.write("%s;A>G;%s;%s;%s;%s;%s;%s\n" % (i, a2g, cytba[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if a2c != 0:
		rowcod = pern(cytba[cs:cs+3],posincod,"C")
		res = synornot(rowcod)
		gens.write("%s;A>C;%s;%s;%s;%s;%s;%s\n" % (i, a2c, cytba[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if t2a != 0:
		rowcod = pern(cytba[cs:cs+3],posincod,"A")
		res = synornot(rowcod)
		gens.write("%s;T>A;%s;%s;%s;%s;%s;%s\n" % (i, t2a, cytba[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if t2g != 0:
		rowcod = pern(cytba[cs:cs+3],posincod,"G")
		res = synornot(rowcod)
		gens.write("%s;T>G;%s;%s;%s;%s;%s;%s\n" % (i, t2g, cytba[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if t2c != 0:
		rowcod = pern(cytba[cs:cs+3],posincod,"C")
		res = synornot(rowcod)
		gens.write("%s;T>C;%s;%s;%s;%s;%s;%s\n" % (i, t2c, cytba[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if g2a != 0:
		rowcod = pern(cytba[cs:cs+3],posincod,"A")
		res = synornot(rowcod)
		gens.write("%s;G>A;%s;%s;%s;%s;%s;%s\n" % (i, g2a, cytba[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if g2t != 0:
		rowcod = pern(cytba[cs:cs+3],posincod,"T")
		res = synornot(rowcod)
		gens.write("%s;G>T;%s;%s;%s;%s;%s;%s\n" % (i, g2t, cytba[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if g2c != 0:
		rowcod = pern(cytba[cs:cs+3],posincod,"C")
		res = synornot(rowcod)
		gens.write("%s;G>C;%s;%s;%s;%s;%s;%s\n" % (i, g2c, cytba[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if c2a != 0:
		rowcod = pern(cytba[cs:cs+3],posincod,"A")
		res = synornot(rowcod)
		gens.write("%s;C>A;%s;%s;%s;%s;%s;%s\n" % (i, c2a, cytba[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if c2t != 0:
		rowcod = pern(cytba[cs:cs+3],posincod,"T")
		res = synornot(rowcod)
		gens.write("%s;C>T;%s;%s;%s;%s;%s;%s\n" % (i, c2t, cytba[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if c2g != 0:	
		rowcod = pern(cytba[cs:cs+3],posincod,"G")
		res = synornot(rowcod)
		gens.write("%s;C>G;%s;%s;%s;%s;%s;%s\n" % (i, c2g, cytba[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))














gens=open("cytbb.csv","w")

for i in range(1045):





	a2t = 0
	a2g = 0
	a2c = 0
	t2a = 0
	t2g = 0
	t2c = 0
	g2a = 0
	g2t = 0
	g2c = 0
	c2a = 0
	c2t = 0
	c2g = 0

	allo = open("engraulis_Data/cytb_B.phy",'r')

	count = 0

	read = allo.readline().strip()

	while True:

		read = allo.readline().strip()

		if read == "":
			break

		row = read.split("  ")

		if row[1][i] != cytbb[i]:


			codnum = (i) // 3 
			posincod = (i) % 3

			cs = codnum * 3


			#print(row[1][i-1:i+1],cytba[i-1:i+1],row[1][cs:cs+3],cytba[cs:cs+3])
			
			cytcod = synornot(cytbb[cs:cs+3])




			if cytbb[i] == "A" and row[1][i] == "T":
				a2t += 1
			elif cytbb[i] == "A" and row[1][i] == "G":
				a2g += 1
			elif cytbb[i] == "A" and row[1][i] == "C":
				a2c += 1
			elif cytbb[i] == "T" and row[1][i] == "A":
				t2a += 1
			elif cytbb[i] == "T" and row[1][i] == "G":
				t2g += 1
			elif cytbb[i] == "T" and row[1][i] == "C":
				t2c += 1
			elif cytbb[i] == "G" and row[1][i] == "A":
				g2a += 1
			elif cytbb[i] == "G" and row[1][i] == "T":
				g2t += 1
			elif cytbb[i] == "G" and row[1][i] == "C":
				g2c += 1
			elif cytbb[i] == "C" and row[1][i] == "A":
				c2a += 1
			elif cytbb[i] == "C" and row[1][i] == "T":
				c2t += 1
			elif cytbb[i] == "C" and row[1][i] == "G":
				c2g += 1

	if a2t != 0:
		rowcod = pern(cytbb[cs:cs+3],posincod,"T")
		res = synornot(rowcod)
		gens.write("%s;A>T;%s;%s;%s;%s;%s;%s\n" % (i, a2t, cytbb[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if a2g != 0:
		rowcod = pern(cytbb[cs:cs+3],posincod,"G")
		res = synornot(rowcod)
		gens.write("%s;A>G;%s;%s;%s;%s;%s;%s\n" % (i, a2g, cytbb[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if a2c != 0:
		rowcod = pern(cytbb[cs:cs+3],posincod,"C")
		res = synornot(rowcod)
		gens.write("%s;A>C;%s;%s;%s;%s;%s;%s\n" % (i, a2c, cytbb[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if t2a != 0:
		rowcod = pern(cytbb[cs:cs+3],posincod,"A")
		res = synornot(rowcod)
		gens.write("%s;T>A;%s;%s;%s;%s;%s;%s\n" % (i, t2a, cytbb[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if t2g != 0:
		rowcod = pern(cytbb[cs:cs+3],posincod,"G")
		res = synornot(rowcod)
		gens.write("%s;T>G;%s;%s;%s;%s;%s;%s\n" % (i, t2g, cytbb[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if t2c != 0:
		rowcod = pern(cytbb[cs:cs+3],posincod,"C")
		res = synornot(rowcod)
		gens.write("%s;T>C;%s;%s;%s;%s;%s;%s\n" % (i, t2c, cytbb[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if g2a != 0:
		rowcod = pern(cytbb[cs:cs+3],posincod,"A")
		res = synornot(rowcod)
		gens.write("%s;G>A;%s;%s;%s;%s;%s;%s\n" % (i, g2a, cytbb[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if g2t != 0:
		rowcod = pern(cytbb[cs:cs+3],posincod,"T")
		res = synornot(rowcod)
		gens.write("%s;G>T;%s;%s;%s;%s;%s;%s\n" % (i, g2t, cytbb[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if g2c != 0:
		rowcod = pern(cytbb[cs:cs+3],posincod,"C")
		res = synornot(rowcod)
		gens.write("%s;G>C;%s;%s;%s;%s;%s;%s\n" % (i, g2c, cytbb[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if c2a != 0:
		rowcod = pern(cytbb[cs:cs+3],posincod,"A")
		res = synornot(rowcod)
		gens.write("%s;C>A;%s;%s;%s;%s;%s;%s\n" % (i, c2a, cytbb[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if c2t != 0:
		rowcod = pern(cytbb[cs:cs+3],posincod,"T")
		res = synornot(rowcod)
		gens.write("%s;C>T;%s;%s;%s;%s;%s;%s\n" % (i, c2t, cytbb[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
	if c2g != 0:	
		rowcod = pern(cytbb[cs:cs+3],posincod,"G")
		res = synornot(rowcod)
		gens.write("%s;C>G;%s;%s;%s;%s;%s;%s\n" % (i, c2g, cytbb[cs:cs+3], rowcod, cytcod, res, dia(cytcod,res)))
