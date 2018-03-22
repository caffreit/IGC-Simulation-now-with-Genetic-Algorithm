# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 10:28:25 2016

@author: Administrator
"""

import math
import numpy as np
import scipy.stats
from scipy.optimize import curve_fit
#import csv
from  more_itertools import unique_everseen
import matplotlib.pyplot as plt

def func(cov,gam0,gam1,mid,p):
    return ((gam0-gam1)/(1+((cov/mid)**p)))+gam1
popt_list = []
SE_profile_list = []

split_list = [0.01]#, 0.02, 0.05, 0.1, 0.15]#DEfect in Basal, i.e. % of Not Edge that is defective
Temp_list = [25]#, 0.96, 0.955, 0.95, 0.945, 0.94, 0.92, 0.9, 0.85, 0.8, 0.75, 0.7, 0.6, 0.5, 0.4]
for split in split_list:
    for T in Temp_list:
        basal_pop_list = [0.99, 0.97, 0.96, 0.95, 0.92, 0.9, 0.85, 0.8]#, 0.75]#, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001]
        for basal_pop in basal_pop_list:
            Population_3 = round(split*basal_pop,5) # Defect population
            Population_2 = basal_pop-Population_3 #Bulk Population
            Population_1 = 1 - basal_pop # Edge population
            Pop_list =[]
            SE_list = []
            
    #Defect        ####################################################################################################################
            sigma_3 = 0.005
            mu_3 = 0.09
            step_3 = 0.005
            
            SE_array_3 = np.arange(mu_3-(5*sigma_3),5*sigma_3 + mu_3+step_3,step_3)
            for s in SE_array_3:
                if s < 0:
                    SE_array_3 = np.delete(SE_array_3,0) 
            for s in SE_array_3:
                s = round(s,4)
                SE_list.append(s)
            
            Initial_Pop_list_3 = []
            for x in SE_array_3:
                k = scipy.stats.norm(mu_3, sigma_3).pdf(x)
                Initial_Pop_list_3.append(k*step_3)
    
            sum_pop = sum(Initial_Pop_list_3)
            for p in Initial_Pop_list_3:
                Pop_list.append(p*Population_3/sum_pop) # remember to put*Population_3 back in
            
    #Bulk        #################################################################################################################
            sigma_2 = 0.005
            mu_2 = 0.04
            step_2 = 0.005
            
            SE_array_2 = np.arange(mu_2-(5*sigma_2),5*sigma_2 + mu_2+step_2,step_2)
            for s in SE_array_2:
                if s < 0:
                    SE_array_2 = np.delete(SE_array_2,0)
            for s in SE_array_2:
                s = round(s,4)
                SE_list.append(s)
                
            Initial_Pop_list_2 = []
            for x in SE_array_2:
                k = scipy.stats.norm(mu_2, sigma_2).pdf(x)
                Initial_Pop_list_2.append(k*step_2)
            
            sum_pop = sum(Initial_Pop_list_2)
            for p in Initial_Pop_list_2:
                Pop_list.append(p*Population_2/sum_pop)
            
    #Edge        ###############################################################################################################    
            sigma_1 = 0.005
            mu_1 = 0.09
            step_1 = 0.005
                
            SE_array_1 = np.arange(mu_1 - (5*sigma_1), 5*sigma_1 + mu_1 + step_1, step_1)
            for s in SE_array_1:
                if s < 0:
                    SE_array_1 = np.delete(SE_array_1,0)
            for s in SE_array_1:
                s = round(s,4)
                SE_list.append(s)
            
            Initial_Pop_list_1 = []
            for x in SE_array_1:
                k = scipy.stats.norm(mu_1, sigma_1).pdf(x)
                Initial_Pop_list_1.append(k*step_1*Population_1)
                
            sum_pop = sum(Initial_Pop_list_1)
            for p in Initial_Pop_list_1:
                Pop_list.append(p*Population_1/sum_pop)
            dsfg = tuple(Pop_list)
            D_Pop_list = list(dsfg)
                
    #Hexane    ###############################################################################################################
                
            SE_Hexane = 0.0184
            area_Hexane = 5.15E-19
            TinC = T #Temp in Celsius
            kT = 1.38e-23*(273.15 + TinC)
            
            Hexane_P_list = []
            for SE in SE_list:
                free_energy_absorption = -2 * area_Hexane * math.sqrt(SE * SE_Hexane)
                Hexane_P = math.exp(-(free_energy_absorption/kT))
                Hexane_P_list.append(Hexane_P)
            
            Hexane_no_duplicates = list(unique_everseen(Hexane_P_list))
            Hexane_Ptotal = sum(Hexane_no_duplicates) 
            
            Hexane_Energy_prob_list = []
            for Hexane_P in Hexane_P_list:
                Hexane_Energy_prob_list.append(Hexane_P/Hexane_Ptotal)
            
            ##############################################################################################################
            
            granularity = 1000
            inject_size = float(1)/granularity    
            coverage = []
            Hexane_SE_measured_list = []
            Hexane_SE_actual_list = [] 
            injected_amount = 0
            
            for i in range(granularity-1):
                H_PreNorm_probs = [a*b for a,b in zip(Pop_list,Hexane_Energy_prob_list)]
                norm_const = 1/sum(H_PreNorm_probs)
            
                Postnorm_prob_list = []
                for j in H_PreNorm_probs:
                    postnorm_prob = norm_const*j
                    Postnorm_prob_list.append(postnorm_prob)
                    
                Probe_SE_list = [c*d for c,d in zip(Postnorm_prob_list,SE_list)]
                SE_meas = sum(Probe_SE_list)
                Hexane_SE_measured_list.append(SE_meas)
                SE_actual = np.mean(Hexane_SE_measured_list)
                Hexane_SE_actual_list.append(SE_actual)
                
                    
                for Pop_prob, post_norm_prob in zip(Pop_list, Postnorm_prob_list):
                    Pop_prob = ((Pop_prob * (float(1) - injected_amount)) - (inject_size * post_norm_prob)) / (float(1)-(injected_amount+inject_size))
                    Pop_list.pop(0)
                    if Pop_prob < 0:
                        Pop_prob = 0
                    Pop_list.append(Pop_prob)
                    
                injected_amount += inject_size
                injected_amount = round(injected_amount,5)
                coverage.append((float(i)+1)/granularity)
                print(sum(Pop_list))
    
            
    #Nonane  ################################################################################################################
            
            SE_Decane = 0.0227
            area_Decane = 6.9E-19
            
            Decane_P_list = []
            for SE in SE_list:
                D_free_energy_absorption = -2 * area_Decane * math.sqrt(SE * SE_Decane)
                Decane_P = math.exp(-(D_free_energy_absorption/kT))
                Decane_P_list.append(Decane_P)
            
            Decane_no_duplicates = list(unique_everseen(Decane_P_list))
            Decane_Ptotal = sum(Decane_no_duplicates) 
            
            Decane_Energy_prob_list = []
            for Decane_P in Decane_P_list:
                Decane_Energy_prob_list.append(Decane_P/Decane_Ptotal)
            
            ##############################################################################################################
                
            D_coverage = []
            Decane_SE_measured_list = []
            Decane_SE_actual_list = [] 
            D_injected_amount = 0
            
            for i in range(granularity-1):
                D_PreNorm_probs = [a*b for a,b in zip(D_Pop_list,Decane_Energy_prob_list)]
                D_norm_const = 1/sum(D_PreNorm_probs)
            
                D_Postnorm_prob_list = []
                for j in D_PreNorm_probs:
                    D_postnorm_prob = D_norm_const*j
                    D_Postnorm_prob_list.append(D_postnorm_prob)
                    
                D_Probe_SE_list = [c*d for c,d in zip(D_Postnorm_prob_list,SE_list)]
                D_SE_meas = sum(D_Probe_SE_list)
                Decane_SE_measured_list.append(D_SE_meas)
                D_SE_actual = np.mean(Decane_SE_measured_list)
                Decane_SE_actual_list.append(D_SE_actual)
                
                    
                for D_Pop_prob, D_post_norm_prob in zip(D_Pop_list, D_Postnorm_prob_list):
                    D_Pop_prob = ((D_Pop_prob * (float(1) - D_injected_amount)) - (inject_size * D_post_norm_prob)) / (float(1)-(D_injected_amount+inject_size))
                    D_Pop_list.pop(0)
                    if D_Pop_prob < 0:
                        D_Pop_prob = 0
                    D_Pop_list.append(D_Pop_prob)
                    
                D_injected_amount += inject_size
                D_injected_amount = round(D_injected_amount,5)
                D_coverage.append((float(i)+1)/granularity)
                print(sum(D_Pop_list))
                
    ########################################################################################################################            
            SE_profile = []
            RTlnVn_hex = []        
            RTlnVn_dec = []
            
            for hex_experienced_nrj in Hexane_SE_actual_list:
                hex_rtlnvn = 84136*(np.sqrt(hex_experienced_nrj))
                RTlnVn_hex.append(hex_rtlnvn)
                
            for dec_experienced_nrj in Decane_SE_actual_list:
                dec_rtlnvn = 125208*(np.sqrt(dec_experienced_nrj))
                RTlnVn_dec.append(dec_rtlnvn)
                        
            y2_y1 = [f-g for f,g in zip(RTlnVn_dec,RTlnVn_hex)]
            x2_x1 = 3#4.48e-20
            for p in y2_y1:
                m = p/x2_x1
                n = (m*m)/1.86e8
                SE_profile.append(n)
            SE_profile_list.append(SE_profile)
            try:
                popt = list(curve_fit(func,coverage,Hexane_SE_actual_list))
            except:
                pass
                #popt = [np.array((0,0,0,0)),np.array((0,0,0,0))]
            popt_list.append(popt[0])
#        l = zip(*SE_profile_list)
#        with open(str(split)+".csv","wb") as g:
#            writer = csv.writer(g)
#            writer.writerows(l)
            
    popt_array = np.array(popt_list)
    gamma0_list = popt_array[:,0]
    gamma1_list = popt_array[:,1]
    midpoint_list = popt_array[:,2]
    p_list = popt_array[:,3]
    delta_gamma = np.subtract(gamma1_list,gamma0_list)
    
    
    plt.figure(figsize=(15,10))
    plt.suptitle("defect="+str(split)+", sigma_def="+str(sigma_3)+", "+"mu_def="+str(mu_3)+", "+"sigma_bulk="+str(sigma_2)+", "+"mu_bulk="+str(mu_2)+", "+"sigma_edge="+str(sigma_1)+", "+"mu_edge="+str(mu_1))
    plt.subplot(321)
    plt.plot(gamma0_list,gamma1_list,"ro")
    plt.ylabel("gamma1")
    plt.xlabel("gamma0")
    
    plt.subplot(322)
    plt.plot(gamma1_list,delta_gamma, "bo")
    plt.ylabel("delta gamma")
    plt.xlabel("gamma1")
    
    plt.subplot(323)
    plt.plot(gamma0_list,delta_gamma,"ro")
    plt.ylabel("delta gamma")
    plt.xlabel("gamma0")
    
    plt.subplot(324)
    plt.plot(midpoint_list,gamma1_list,"ro")
    plt.ylabel("gamma1")
    plt.xlabel("midpoint_list")
    
    plt.subplot(325)
    plt.plot(midpoint_list,gamma0_list,"ro")
    plt.ylabel("gamma0")
    plt.xlabel("midpoint_list")
    
    plt.subplot(326)
    plt.plot(midpoint_list,delta_gamma,"ro")
    plt.ylabel("delta gamma")
    plt.xlabel("midpoint_list")
    
    #plt.savefig("C:\Users\Administrator\Desktop\plots from igc sim" + "\defect"+str(int(split*100))+"sigma_def"+str(int(sigma_3*1000))+"mu_def"+str(int(mu_3*1000))+"sigma_bulk"+str(int(sigma_2*1000))+"mu_bulk"+str(int(mu_2*1000))+"sigma_edge"+str(int(sigma_1*1000))+"mu_edge"+str(int(mu_1*1000)))
    #plt.show()
