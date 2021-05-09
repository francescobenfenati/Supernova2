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

        self.tshitsdf=pd.DataFrame(columns=['DOM','time']) #contains times of all hits in a TS                                                                                                           
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
        pmt1_12 = dom1[dom1.channel_id == 12]
        print(pmt1_12)

        t = pmt1_12.time[0]
        #print('DOM1 channel 12 hit time  = ',t)

        self.tshitsdf.loc[len(self.tshitsdf)+1] = [1,t]

        #cycle over the remaining  DOMS 
        now = datetime.now()

        #-----------------------------------------------------------------------------------------------
        #take all the chosen-DOM1-PMT hits' times and calculate dt with all chosen DOM-PMT hits' times, doing this for all DOMs
        

        for i in self.orderedDOM[1:]:
            #to skip DOM 3,4,5,7,8,9,10,15,17
            if i == 818848239 or i == 818798906 or i == 818848226 or i == 818806250 or i == 818785829 or i == 818785841 or i == 818806263 or i == 818785853 or i == 818848251:
                continue
            dom_number = self.Dom_id_name[str(i)]
            dom2=tshits[tshits.dom_id == i]
            #d = np.unique(dom2['channel_id'])
            #d.sort()
            #print(d)
                
            #fixing for all the LARGE Board
            pmtd_12 = dom2[dom2.channel_id == 12]
            
            #print('-.-.-.-.-.-.-.-.-.-.-.-')
            #print('DOM_',dom_number)

            tt = pmtd_12.time[0]
            self.tshitsdf.loc[len(self.tshitsdf)+1] = [dom_number, tt]

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

        #save zoomdf
        if TSindex%3600==0:
            print('saving dataframes')
            self.tshitsdf.to_csv("/home/km3net/analysis/Phase2/Dataframes/tshits/tshitsdf_"+str(TSindex)+".csv",mode='a',header=True)
            self.delaysdf.to_csv("/home/km3net/analysis/Phase2/Dataframes/delays/delaysdf_"+str(TSindex)+".csv",mode='a',header=True)
            print("dataframes saved")
        
        #plot at each hour if OOS not occurred                                                                                                                                                      
        #if TSindex%3600==0:
        #    print('FINALLY PLOTTING!')
        #    self.plotter()
        else:
            return blob
        
    def over_threshold(self):
        #Threshold delays (1 micro, 100 ns, 50 ns)
        if abs(self.final_delay) > 1000:
            Over1micro = self.testdf.loc[len(self.testdf)].to_csv('Ph2_Over1micro.csv',mode='a',header=False) 
            self.status=1
        elif abs(self.final_delay) > 100:
            Over100nano = self.testdf.loc[len(self.testdf)].to_csv('Ph2_Over100nano.csv',mode='a',header=False) 
            self.status=1
        elif abs(self.final_delay) > 50:
            Over50nano = self.testdf.loc[len(self.testdf)].to_csv('Ph2_Over50nano.csv',mode='a',header=False) 
            self.status=1
       # else:
           # print("Delay is smaller than 50 ns")

    def summary(self,i):
        #le calcoliamo come valore di riferimento sulle prime N TS, poi li possiamo aggiornare per printare                          
        df_dom = self.testdf[(self.testdf.DOMnumber == i)]
        MEAN = df_dom["deltaT"].mean()
        MAX =  df_dom["deltaT"].max()
        MIN =  df_dom["deltaT"].min()
        domname = 'DOM'+str(i+1)
        self.SUMMARY.loc[len(self.SUMMARY)+1] = [now,domname, MEAN, MAX,MIN]

    def plot_delays(self,dom,delays):
        fig,ax = plt.subplots()
        #if dom == '2' or dom == '11':
        #    x_entries,binning,_ = plt.hist(delays,bins=100,range = [-0.2e6,+0.0e6],density=0,alpha=0.5,color='black')
        #elif dom == '5' or dom == '9' or dom == '6' or dom == '7' or dom == '10' or dom == '17':
        #    x_entries,binning,_ = plt.hist(delays,bins=100,range = [0.0e6,+0.1e6],density=0,alpha=0.5,color='black')
        #elif dom == '8' or dom == '18':
        #    x_entries,binning,_ = plt.hist(delays,bins=100,range = [-0.1e6,0.0e6],density=0,alpha=0.5,color='black')
        #elif dom == '14' or dom == '12':
        #    x_entries,binning,_ = plt.hist(delays,bins=100,range = [-0.2e6,-0.1e6],density=0,alpha=0.5,color='black')            
        #elif dom == '16':
         #   x_entries,binning,_ = plt.hist(delays,bins=100,range = [0.1e6,0.2e6],density=0,alpha=0.5,color='black')
 
        x_entries,binning,_ = plt.hist(delays,bins=600,range = [0,+300e3] ,density=0,alpha=0.5,color='black')
        ax.set_xlabel("Delta_t (ns)")
        ax.set_ylabel("Occurrence")
        fig.suptitle("DOM-"+dom)
        namefig = "/home/km3net/analysis/Phase2/Images/Single_TS/DOM_"+dom
        fig.savefig(namefig)
        
    def plotter(self):
        now = datetime.now()
        df_list = []
        missing_dom=18-self.numberofactivedom
        name=[]
        Dom_number=[] #here you can get rid of Dom_number somehow, use self.doms directly instead
        print('doms ',self.doms)
        print('number of active doms ',self.numberofactivedom)
        for j in self.doms:
            if j == 1: #skip DOM 1
                continue
            Dom_number.append(j) #getting the doms in decimal base
        
        Dom_number.sort()
        #Dom_number=Dom_number[2:] #exclude DOM 1
        print('Dom number ',Dom_number)
        for aa in Dom_number:
            partial = self.testdf[(self.testdf.DOMnumber == aa)].tail(60)
            name.append('CLB_'+str(aa))
            # make a list of all dataframes in order                             
            df_list.append(partial)
        if missing_dom!=0:    
            for i in range(missing_dom):
                partial = self.testdf[(self.testdf.DOMnumber == 1)].tail(60)
                df_list.append(partial)
                name.append('CLB_'+str(1))
        #start of the PLOT with 16 subplots
        time = now.strftime("%y_%m_%d_%H_%M")
        
        fig = plt.figure(figsize=(10, 8))
        outer = gridspec.GridSpec(4, 4, wspace=0.3, hspace=0.3)
        for i,j in enumerate(Dom_number):
            inner = gridspec.GridSpecFromSubplotSpec(1, 2,subplot_spec=outer[i], wspace=0.1, hspace=0.1)
            axe = plt.Subplot(fig,inner[0])
            if i<12:
                axe.set_xticks([])
            print('i = ',i)
            print('dom number ', j)
            l=df_list[i].plot(x ='Time', y='deltaT', kind='scatter', ax=axe)
            axe2 = plt.Subplot(fig, inner[1])
            axe2.tick_params(colors='red')
            if i==3 or i == 7 or i == 11 or i == 15:
                axe2.yaxis.tick_right()
                #axe2.tick_params(right=True, labelright=True,colors='red')
            else:
                axe2.set_yticks([])
            
            g=self.testdf[(self.testdf.DOMnumber == j)].hist(column='deltaT',ax=axe2,color = "red")
            fig.add_subplot(axe)
            fig.add_subplot(axe2)
            titlename = name[i]
            axe.set_title(titlename)
            axe.set_xlabel('')
            axe.set_ylabel('')
            axe2.set_title('')
            if i==0:
                axe.set_ylabel('Delay (ns)')
            elif i==15:
                axe2.set_xlabel('TIME (s)')
            
        fig.suptitle('DELAYS_referred_CLB1_at:'+time,y=0.05,x=0.5)
        fig.tight_layout()
        namefigure="/home/km3net/analysis/Phase2/Images/Multi_TS/Delays_DOMs_"+time
        fig.savefig(namefigure)
        self.testdf.sort_values(by = ['DOMnumber','Time']).tail(840).to_csv("/home/km3net/analysis/Phase2/Dataframes/Overtime/Delays_DOMs_"+time,index=False)
        self.testdf = self.testdf.iloc[0:0]

    def finish(self):
        self.plotter()
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
