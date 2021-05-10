import numpy as np
import io
import os
import sys
import time
import matplotlib
matplotlib.use('Agg')
import km3pipe as kp
import km3modules as km
from km3pipe.io.daq import TMCHData
from km3pipe import Module
from tabulate import tabulate
#from km3pipe.core import Pump
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


class OOSAnalyzer(Module):
    def configure(self):

        #----------------------------------------------------------------                                                                                                                                                                                            
        #conversion dom_id to dom_number for  phase2 @BCI                                                                                                                                                                                                            
        #----------------------------------------------------------------                                                                                                                                                                                            
        self.orderedDOM = [818798894,818845349,818848239,818798906,818848226,818798740,818806250,818785829,818806263,818785853,818798725,818848646,818798728,818798885,818848251,818845350,818785841,818806251]

        self.decimal = list(range(1, 19)) 
        self.Dom_id_name = {str(self.orderedDOM[i]) : self.decimal[i] for i in range(len(self.orderedDOM))}

        self.doms=[]
        self.TOOS=0
        self.numberofactivedom=0
        self.status=0
        self.timetoplot=60
        self.final_delay=0

        self.tshitsdf=pd.DataFrame(columns=['DOM','TSCounter','time']) #contains times of all hits in a TS and the Timeslice timestamp                                                                                                          
        self.delaysdf=pd.DataFrame(columns=['DOM','delay','Time']) #contains all delays of the delay distribution 
        
        self.counter=0
        self.start = round(time.time(),0)
        self.aa = 0

    def process(self,blob):
        info = blob['TimesliceInfo']
        info2=blob['TimesliceFrameInfos']

        self.counter+=1
        TSindex=info['timestamp'][0]
        Time = info['nanoseconds'][0]/1000000000
        TSCounter = TSindex+Time

        self.status=0
        
        tshits = blob['TSHits']

        #to exlude dom 3,5,8... other doms do not give L1 hits so there is no need of excluding them
        tshits = tshits[tshits.dom_id != 818848226]
        tshits = tshits[tshits.dom_id != 818785829]
        tshits = tshits[tshits.dom_id != 818848239]
        
        
        self.doms = np.unique(tshits['dom_id']) 
        self.doms = [self.Dom_id_name[str(i)] for i in self.doms]
        self.doms.sort()
        #print("Doms in the timeslice are ",self.doms)

        if len(self.doms)== 0:
            return blob

        #print('DOMS =', self.doms)
        #print('len of doms ',len(self.doms))
        channels = np.unique(tshits['channel_id'])
        #print('firing PMTs = ',channels)
        #print('number of firing PMTs = ',len(channels))
        print('TSCounter ', TSCounter)
        
        #fixing DOM 1 Large 
        dom1=tshits[tshits.dom_id == 818798894]

        #dom17=tshits[tshits.dom_id == 818785841]
        #print("dom 17", dom17)
        
        pmt1_12 = dom1[dom1.channel_id == 12]

        if len(pmt1_12.time) == 0: #if for some reasons the OP skips the hit injection for that TS, then skip the TS
            print("Hit missing for DOM1 !")
            print("L1 DOMs in this TS = ",self.doms)
            print(pmt1_12)
            return blob

        t = pmt1_12.time[0]
        #print('DOM1 channel 12 hit time  = ',t)

        self.tshitsdf.loc[len(self.tshitsdf)+1] = [1,TSCounter,t]

        #cycle over the remaining  DOMS 
        now = datetime.now()

        #-----------------------------------------------------------------------------------------------
        #take all the chosen-DOM1-PMT hits' times and calculate dt with all chosen DOM-PMT hits' times, doing this for all DOMs
        

        for i in self.orderedDOM[1:]:
            #to skip DOM 3,5,7,8,9,15,17
            if i == 818848239 or i == 818848226 or i == 818806250 or i == 818785829 or i == 818806263 or i == 818848251:
                continue
            dom_number = self.Dom_id_name[str(i)]
            dom2=tshits[tshits.dom_id == i]
            #d = np.unique(dom2['channel_id'])
            #d.sort()
            #print(d)
                
            #fixing for all the LARGE Board
            pmtd_12 = dom2[dom2.channel_id == 12]

            if len(pmtd_12.time) == 0: #if the hit is missing, skip to next DOM
                print("Hit missing for DOM",dom_number,"!")
                print(pmtd_12)
                continue

            #print('-.-.-.-.-.-.-.-.-.-.-.-')
            #print('DOM_',dom_number)

            tt = pmtd_12.time[0]
            self.tshitsdf.loc[len(self.tshitsdf)+1] = [dom_number,TSCounter,tt]

            deltat = tt-t
            self.delaysdf.loc[len(self.delaysdf)+1] = [dom_number,deltat,TSCounter]
        
            #print('Delay:',deltat)
            #print('-.-.-.-.-.-.-.-.-.-.-.-')
           
            #Implementig circular buffer on python dataframe
            # self.over_threshold()
            # print(self.testdf)
               
            # self.plot_delays(str(Dom_numbers),delay_over_TS)    

            #if TSindex%3600==0:
            #    self.plot_delays(str(Dom_numbers),delay_over_TS)    

        #save dataframes to csv every 1 hour
        if TSindex%3600==0:
            print('saving dataframes')
            self.tshitsdf.to_csv("/home/km3net/analysis/Phase2/Dataframes/tshits/tshitsdf_"+str(TSindex)+".csv",mode='a',header=True)
            self.delaysdf.to_csv("/home/km3net/analysis/Phase2/Dataframes/delays/delaysdf_"+str(TSindex)+".csv",mode='a',header=True)
            print("dataframes saved")
            self.tshitsdf = self.tshitsdf.iloc[0:0]
            self.delaysdf = self.delaysdf.iloc[0:0]

        #plot at each hour if OOS not occurred                                                                                                                                                      
        #if TSindex%3600==0:
        #    print('FINALLY PLOTTING!')
        #    self.plotter()
        else:
            return blob
        
    def finish(self):
        print("killed CTRL_C")
    
    
def main():
    pipe=kp.Pipeline()
    pipe.attach(kp.io.ch.CHPump,host='192.168.0.20', port=5553,tags='IO_TSL1',timeout=60 * 5,max_queue=2000)
    pipe.attach(kp.io.daq.TimesliceParser)
    pipe.attach(OOSAnalyzer)
    pipe.drain()

    #/////to store the last 3h dataframe into a hdf5 file/////
    #store = HDFStore('self.SUMMARY_LAST3H.h5')
    #store.put('self.testdf',self.testdf)
    #store.close()

if __name__ == '__main__':
    main()
