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
		dia = True
	else:
		dia = False

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
	       "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
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















gens = open("cytba_diff.csv","w")
cytba = open("cytba.txt",'r').readline().strip()
cytbb = open("cytbb.txt",'r').readline().strip()




# a2t = 0
# a2g = 0
# a2c = 0
# t2a = 0
# t2g = 0
# t2c = 0
# g2a = 0
# g2t = 0
# g2c = 0
# c2a = 0
# c2t = 0
# c2g = 0

a = 0
t = 0
g = 0
c = 0

ao = 0
to = 0
go = 0
co = 0

count = 0

for i in range(len(cytba)):

	if cytba[i] != cytbb[i]:


		codnum = (i) // 3 
		posincod = (i) % 3

		cs = codnum * 3

			
		cytcoda = synornot(cytba[cs:cs+3])
		cytcodb = synornot(cytbb[cs:cs+3])		

		


		gens.write("%s;%s;%s;%s;%s;%s\n" % (i, cytba[cs:cs+3], cytbb[cs:cs+3], cytcoda, cytcodb, dia(cytcoda,cytcodb)))





		if cytba[i] == "A" and cytbb[i] == "T":
			a += 1
			to += 1
		elif cytba[i] == "A" and cytbb[i] == "G":
			a += 1
			go += 1
		elif cytba[i] == "A" and cytbb[i] == "C":
			a += 1
			co += 1
		elif cytba[i] == "T" and cytbb[i] == "A":
			t += 1
			ao += 1
		elif cytba[i] == "T" and cytbb[i] == "G":
			t += 1
			go += 1
		elif cytba[i] == "T" and cytbb[i] == "C":
			t += 1
			co += 1
		elif cytba[i] == "G" and cytbb[i] == "A":
			g += 1
			ao += 1
		elif cytba[i] == "G" and cytbb[i]== "T":
			g += 1
			to += 1
		elif cytba[i] == "G" and cytbb[i] == "C":
			g += 1
			co += 1
		elif cytba[i] == "C" and cytbb[i] == "A":
			c += 1
			ao += 1
		elif cytba[i] == "C" and cytbb[i] == "T":
			c += 1
			to += 1
		elif cytba[i] == "C" and cytbb[i] == "G":
			c += 1
			go += 1



gens.write("%s;%s;%s;%s\n" % (a, t, g, c))
gens.write("%s;%s;%s;%s\n" % (ao, to, go, co))


a = 0
t = 0
g = 0
c = 0

ao = 0
to = 0
go = 0
co = 0

count = 0

for i in range(len(cytba)):

	if cytba[i] != cytbb[i]:
           

		codnum = (i) // 3 
		posincod = (i) % 3

		cs = codnum * 3


		#print(row[1][i-1:i+1],cytba[i-1:i+1],row[1][cs:cs+3],cytba[cs:cs+3])
			
		cytcoda = synornot(cytba[cs:cs+3])
		cytcodb = synornot(cytbb[cs:cs+3])		
		#rowcod = pern(cytba[cs:cs+3],posincod,"T")
		if dia(cytcoda,cytcodb):



			if cytba[i] == "A" and cytbb[i] == "T":
				a += 1
				to += 1
			elif cytba[i] == "A" and cytbb[i] == "G":
				a += 1
				go += 1
			elif cytba[i] == "A" and cytbb[i] == "C":
				a += 1
				co += 1
			elif cytba[i] == "T" and cytbb[i] == "A":
				t += 1
				ao += 1
			elif cytba[i] == "T" and cytbb[i] == "G":
				t += 1
				go += 1
			elif cytba[i] == "T" and cytbb[i] == "C":
				t += 1
				co += 1
			elif cytba[i] == "G" and cytbb[i] == "A":
				g += 1
				ao += 1
			elif cytba[i] == "G" and cytbb[i]== "T":
				g += 1
				to += 1
			elif cytba[i] == "G" and cytbb[i] == "C":
				g += 1
				co += 1
			elif cytba[i] == "C" and cytbb[i] == "A":
				c += 1
				ao += 1
			elif cytba[i] == "C" and cytbb[i] == "T":
				c += 1
				to += 1
			elif cytba[i] == "C" and cytbb[i] == "G":
				c += 1
				go += 1


gens.write("%s;%s;%s;%s\n" % (a, t, g, c))
gens.write("%s;%s;%s;%s\n" % (ao, to, go, co))








a = 0
t = 0
g = 0
c = 0

ao = 0
to = 0
go = 0
co = 0




count = 0

for i in range(len(cytba)):

	if cytba[i] != cytbb[i]:
           

		codnum = (i) // 3 
		posincod = (i) % 3

		cs = codnum * 3


		#print(row[1][i-1:i+1],cytba[i-1:i+1],row[1][cs:cs+3],cytba[cs:cs+3])
			
		cytcoda = synornot(cytba[cs:cs+3])
		cytcodb = synornot(cytbb[cs:cs+3])		
		#rowcod = pern(cytba[cs:cs+3],posincod,"T")




	    # ["S","L","P","R","T","V","A","G"]


	    # "TC":"S",
	    # "CT":"L",
	    # "CC":"P",
	    # "CG":"R",
	    # "AC":"T",
	    # "GT":"V",
	    # "GC":"A",
	    # "GG":"G"

		quadro = "TC", "CT", "CC", "CG", "AC", "GT", "GC", "GG"


		if dia(cytcoda,cytcodb):

			if cytba[cs:cs+2] in quadro: 

				if cytba[i] == "A" and cytbb[i] == "T":
					a += 1
					to += 1
				elif cytba[i] == "A" and cytbb[i] == "G":
					a += 1
					go += 1
				elif cytba[i] == "A" and cytbb[i] == "C":
					a += 1
					co += 1
				elif cytba[i] == "T" and cytbb[i] == "A":
					t += 1
					ao += 1
				elif cytba[i] == "T" and cytbb[i] == "G":
					t += 1
					go += 1
				elif cytba[i] == "T" and cytbb[i] == "C":
					t += 1
					co += 1
				elif cytba[i] == "G" and cytbb[i] == "A":
					g += 1
					ao += 1
				elif cytba[i] == "G" and cytbb[i]== "T":
					g += 1
					to += 1
				elif cytba[i] == "G" and cytbb[i] == "C":
					g += 1
					co += 1
				elif cytba[i] == "C" and cytbb[i] == "A":
					c += 1
					ao += 1
				elif cytba[i] == "C" and cytbb[i] == "T":
					c += 1
					to += 1
				elif cytba[i] == "C" and cytbb[i] == "G":
					c += 1
					go += 1


gens.write("%s;%s;%s;%s\n" % (a, t, g, c))
gens.write("%s;%s;%s;%s\n" % (ao, to, go, co))