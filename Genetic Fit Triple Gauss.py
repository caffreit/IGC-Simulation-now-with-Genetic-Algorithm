# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 14:08:50 2017

@author: Administrator
"""

import math
import numpy as np
import scipy.stats
#from scipy.optimize import curve_fit
#import csv
from  more_itertools import unique_everseen
import random as rd



def func(cov,gam0,gam1,mid,p):
    return ((gam0-gam1)/(1+((cov/mid)**p)))+gam1

def population(count):
    individual_array = []
    for i in range(count):
        individual = [float(rd.randint(0,50))/100, float(rd.randint(10,150))/1000, float(rd.randint(10,150))/1000, float(rd.randint(10,150))/1000, float(rd.randint(1,30))/1000, float(rd.randint(1,30))/1000, float(rd.randint(1,30))/1000]
        #Basal Defect %, centre; Defect, Bulk, Edge, Width; Defect, Bulk, Edge
        individual_array.append(individual)
    return individual_array
    
    
def AcceptReject(child_array,index):
    while (True):
        number = rd.randint(0,100)
        number2 = rd.uniform(0,fitness_array[par_4_ind])
        if number<85:
            if number%3 ==0:
                if number2>fitness_array[1]:
                    return(child.append(father[i]))
            elif number%3 == 1:
                 if number2>fitness_array[2]:
                     return(child.append(parent3[i]))
            else:
                 if number2>fitness_array[0]:
                     return(child.append(mother[i]))
        else:
            if 0 < i < 4:
                new_gene = float(rd.randint(10,150))/1000
            elif 3 < i < 7:
                new_gene = float(rd.randint(1,30))/1000
            else:
                new_gene = float(rd.randint(0,50))/100
            return(child.append(new_gene))


variableses = population(20)

children = [[0],[0],[0],[0],[0],[0],[0],[0],[0],[0]]
fitness_list = []
best_cost = 99999999
best_parent = []

for generation in range(2000):
    for variables in variableses:
        popt_list = []
        SE_profile_list = []
        split = variables[0]
        Temp = 25
        basal_pop_list = [0.999, 0.997, 0.995, 0.992, 0.99, 0.975, 0.95]#, 0.85, 0.8]#, 0.75]#, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001]
        for basal_pop in basal_pop_list:
            Population_3 = round(split*basal_pop,5) # Defect population
            Population_2 = basal_pop-Population_3 #Bulk Population
            Population_1 = 1 - basal_pop # Edge population
            Pop_list =[]
            SE_list = []
            
        #Defect        ####################################################################################################################
            sigma_3 = variables[1]
            mu_3 = variables[4]
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
            sigma_2 = variables[2]
            mu_2 = variables[5]
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
            sigma_1 = variables[3]
            mu_1 = variables[6]
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
            TinC = Temp #Temp in Celsius
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
                #print(sum(Pop_list))
        
            
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
            SE_profile_list.append(Decane_SE_actual_list)
            midpoint_y = (Decane_SE_actual_list[0]+Decane_SE_actual_list[-1])/2
            for i in range(len(Decane_SE_actual_list)):
                if midpoint_y < Decane_SE_actual_list[i]:
                    pass
                else:
                    midpoint_x = Decane_SE_actual_list[i]
            try:
                popt = [np.array((Decane_SE_actual_list[0],Decane_SE_actual_list[-1],midpoint_x,0)),np.array((0,0,0,0))]                
                #poptfit = list(curve_fit(func,coverage,Decane_SE_actual_list))
            except:
                pass
                popt = [np.array((Decane_SE_actual_list[0],Decane_SE_actual_list[-1],midpoint_x,0)),np.array((0,0,0,0))]
            popt_list.append(popt[0])
        #        l = zip(*SE_profile_list)
        #        with open(str(split)+".csv","wb") as g:
        #            writer = csv.writer(g)
        #            writer.writerows(l)
                    
        popt_array = np.array(popt_list)
        gamma0_list = popt_array[:,0]
        gamma1_list = popt_array[:,1]
        midpoint_list = popt_array[:,2]
        delta_gamma = np.subtract(gamma1_list,gamma0_list)
        
        Real_IGC_gam1 = [0.04725723, 0.04471484, 0.03851335, 0.04082745, 0.0434372, 0.03665321, 0.0319563]
        Real_IGC_deltagam = [-0.02717509, -0.03199348, -0.02774597, -0.0316344, -0.01376077, -0.02359752, -0.03283994]
        Real_IGC_gam0 = [0.07443232, 0.07670832, 0.06625932, 0.07246185, 0.05719797, 0.06025073, 0.06479624]
        Real_IGC_midp = [0.03295, 0.02819, 0.01839, 0.00946, 0.0122, 0.00904, 0.00746]
        
        cost = 0
        for i in range(len(basal_pop_list)):
            distance1 = np.sum((np.subtract(Real_IGC_gam1,gamma1_list))**2)
            distance2 = np.sum((np.subtract(Real_IGC_gam0,gamma0_list))**2)
            distance3 = np.sum((np.subtract(Real_IGC_midp,midpoint_list))**2)
            distance4 = np.sum((np.subtract(Real_IGC_deltagam,delta_gamma))**2)
            
            
            #distance1 = np.sqrt(((gamma0_list[i]-Real_IGC_gam0[i])**2)+((gamma1_list[i]-Real_IGC_gam1[i])**2))
            #distance2 = np.sqrt(((gamma0_list[i]-Real_IGC_gam0[i])**2)+(((delta_gamma[i]-Real_IGC_deltagam[i])/10)**2))
            #distance3 = np.sqrt((((midpoint_list[i]-Real_IGC_midp[i])*10)**2)+((gamma1_list[i]-Real_IGC_gam1[i])**2))
            cost += (distance1+distance2+distance3+distance4)
        fitness_list.append(cost**2) ## to help speed up selection by making cost difference larger, also need to add in
        ## some way of choosing fitter parents more often like an acceptReject function.
    
    fitness_array = np.asarray(fitness_list)
    father_index = fitness_array.argsort()[:4][1]
    mother_index = fitness_array.argsort()[:4][0]
    parent3_index = fitness_array.argsort()[:4][2]
    par_4_ind = fitness_array.argsort()[:4][3]  ## just so that parent3 has a chance to mate.
    father = variableses[father_index]
    mother = variableses[mother_index]
    parent3 = variableses[parent3_index]
    
    if best_cost > fitness_list[mother_index]:
        best_cost = np.sqrt(fitness_list[mother_index])
        best_parent = mother
    
    for j in range(20):
        child = []
        for i in range(len(variables)):
            AcceptReject(child,i)
            """
            number = rd.randint(0,100)
            #number1 = rd.randint(0,2)
            #testnum = np.sum(fitness_array)
            number2 = rd.randint(0,fitness_list[par_4_ind]) ## just so that parent3 has a chance to mate.
            if number<85:
                if number%3 ==0:
                    if number2>fitness_list[1]
                    child.append(father[i])
                elif number%3 == 1:
                     if number2>fitness_list[2]
                    child.append(parent3[i])
                else:
                     if number2>fitness_list[0]
                    child.append(mother[i])
            else:
                if 0 < i < 4:
                    new_gene = float(rd.randint(10,150))/1000
                elif 3 < i < 7:
                    new_gene = float(rd.randint(1,30))/1000
                else:
                    new_gene = float(rd.randint(0,50))/100
                child.append(new_gene)
             """  
        children.append(child)
        children.pop(0)
        fitness_list.pop(0)
    variableses = children
        

        
        
        
        
        