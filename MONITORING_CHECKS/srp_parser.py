import numpy as np
import pandas as pd
import datetime
from collections import OrderedDict

f = open("/home/km3net/analysis/MONITORING_CHECKS/Run10min_second","r")
orderedDOM = [806451575,808981684,808447031,808985194,808971330,806451239,808952022,808967370,808489098,808976266,809537142,808984748,808982228,808980464,808976292,809544159,808996919]
doms=[1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18]

dom_numbers = {orderedDOM[i] : doms[i] for i in range(len(doms))}

ciccio=[]
dom_list_upi = [] #list which contains all dom_i lists for all DOMs NOT USUED
keys=[]
lines=f.readlines()

for index in range(len(lines)):
    if "Data" in lines[index]:
        for DOM in orderedDOM:
            if str(DOM) in lines[index+1]:
                upi=lines[index].split(" ")[1].split("d")
                upifinal=upi[0][:-1]
                dom_list_upi.append(upi[0][:-1])
                ciccio.append(DOM)

        
#dic=OrderedDict(zip(orderedDOM,dom_list_upi))
dic=OrderedDict(zip(ciccio,dom_list_upi))
#for index in range(len(lines)):
#    if "Data" in lines[index]:
#        if dic[806451575] in lines[index]:
#            newkey=lines[index].split("-")[3]
#            keys.append(newkey)

keys=['ahrs_yaw\n', 'ahrs_pitch\n', 'ahrs_roll\n', 'ahrs_a[0]\n', 'ahrs_a[1]\n', 'ahrs_a[2]\n', 'ahrs_g[0]\n', 'ahrs_g[1]\n', 'ahrs_g[2]\n', 'ahrs_h[0]\n', 'ahrs_h[1]\n', 'ahrs_h[2]\n', 'temp\n','humid\n']


total_values=[]
time_list=[]
dom_list=[]
date_list = []
#DOM_list=[]

for key in keys:
    key_values=[]
    for index in range(len(lines)):
        for dom in orderedDOM:
            upi=dic[dom]
            if upi in lines[index] and key in lines[index]:
                ind = 0
                while True:
                    ind +=1
                    if 'Data' in lines[index+ind]:
                        break
                    row = lines[index+ind]
                    splitted_row = row.split(" ")
                    time = splitted_row[1]
                    if splitted_row[2] == 'PM':
                        hour = str(int(splitted_row[1].split(':')[0])+12)
                        time = hour+':'+splitted_row[1].split(':')[1]+':'+splitted_row[1].split(':')[2]
                    date = splitted_row[0]
                    value = splitted_row[3]
                    if key == 'temp\n' or key == 'humid\n':
                        value = float(value)/100
                    if key == 'ahrs_pitch\n':
                        time_list.append(time)
                        date_list.append(date)
                        dom_list.append(dom_numbers[dom])
                    key_values.append(value)
                    print("1",len(key_values))
    total_values.append(key_values)
   # print("total",len(total_values))
#UTC_datetime = str(datetime.datetime.utcnow())
#print(UTC_datetime)               
#print(dom_list[:10])
#print(time_list[:10])
#print(len(total_values[0]))               
                    
list_of_tuples = list(zip(dom_list,date_list,time_list,total_values[0],total_values[1],total_values[2],total_values[3],total_values[4],total_values[5],total_values[6],total_values[7],total_values[8],total_values[9],total_values[10],total_values[11],total_values[12],total_values[13]))            
df = pd.DataFrame(list_of_tuples,columns=['DOM','date','time','yaw','pitch','roll','aX','aY','aZ','gX','gY','gZ','hX','hY','hZ','Temp','Humid'])

df.to_csv("/home/km3net/analysis/MONITORING_CHECKS/datalog_parsed_second.csv",index=0)
f.close()
