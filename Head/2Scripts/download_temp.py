#coding=utf8

import os 


import requests

import urllib2
from bs4 import BeautifulSoup




src = open("list.txt",'r')


path3="TempCount.txt"
reps3 = open(path3, "w")



#for r in range(1,maxid+1):
n=1
while True: 
    print n
    read = src.readline()
    if read=='':
        break
    
    sepread=read.split(';')
    
    url=sepread[0]

    #soup = BeautifulSoup(urllib2.urlopen(url))

    headers = {'User-Agent':'Mozilla/5.0'}


    page = requests.get(url)

    
   

    soup = BeautifulSoup(page.text, "html.parser")


    path="DOWN/%s_%s" % (n,soup.title.string)


    f = open(path, 'w')

    print soup.find_all('table')
    raw=soup.find_all('table')[2].find_all('td')

   
    m = 7 # 8 колонка
    while True: 
        try:
            print raw[m]
            reps3.write ("\"%s\" \"%s\"\n" % (soup.title.string,raw[m]))
            m+=14 # 14 колонок всего
        except:
            #print ("Nope")
            break



       

    f.write(soup.prettify().encode('UTF-8'))

    n+=1







with open ('TempCount.txt', 'r') as f:
    old_data = f.read()

new_data = old_data.replace("<td>","").replace("</td>","").replace("Species Population Growth List  -","")

with open ('Temps.txt', 'w') as f:
    f.write(new_data)