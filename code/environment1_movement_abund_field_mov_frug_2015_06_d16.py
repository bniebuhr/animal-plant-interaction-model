# -*- coding: utf-8 -*-
"""
Spatially-explicit animal-plant interactions model

Bernardo Niebuhr
Jun 2015
"""

# Loading libraries
import sys, os
import numpy as np
from scipy.stats import lognorm, gamma
import time
#import matplotlib.pyplot as plt

# Importing modules
#sys.path.append("/media/windows/Users/ber/Documents/Documentos_ber/Atividades 2014/Colaboracao_Ricardo/simulacoes_interacao/aux")
#from aux_functions import distance

def distance(xy_ind_a, xy_ind_b):
    """
    This function returns the euclidean distance among two individuals.
    Input:
    xy_ind_a: coordinates (x,y) of individual a
    xy_ind_b: coordinates (x,y) of individual b
    Output:
    dist: straight line distance among the individuals
    """
    a_x=xy_ind_a[0]
    a_y=xy_ind_a[1]
    b_x=xy_ind_b[0]
    b_y=xy_ind_b[1]
    
    dist = np.sqrt((a_x-b_x)*(a_x-b_x) + (a_y-b_y)*(a_y-b_y))
    return dist

def regrowth_model(fruits_init = 0, max_fruits = 100, rate = 12):
    '''
    This function performs a regrowth of fruits in a plant individual
    
    Input:
     fruits_init: initial number of fruits
     max_fruits: maximum number of fruits in a plant
     rate: rate of fruit growth, in fruits/day
    '''
    return fruits_init + rate*(1 - fruits_init/max_fruits)
    
def read_abund(file_name, col = 0):    
    
    f = open(file_name, 'r')
    lines = f.readlines()
    lines = lines[1:] # exclude headers
    f.close()
    
    mat_names = []
    matrix = []
    for line in lines:
        mat_names.append(line.strip().split(";")[0])
        matrix.append(map(float, line.strip().split(";")[1:]))
        
    #matrix = np.asarray(matrix, dtype='float32')
    mat_names = np.array(mat_names)
    matrix = np.asarray(matrix)
    
    return mat_names, matrix[:,col], matrix[:,-2:]


def create_environment(env_dim = 10000, plant_sp = 3, individuals = 1500, dist_corr = 200.0, rho = 0.0, stage = 0,
                       file_name = ""):
    '''
    This function populates a continuous landscape of dimension env_dim (in meters) with "individuals" 
    plants of "plant_sp" different species. dist_corr and rho sets the autocorrelation in the position
    of plants.
    
    We still need tox':
    - aggregate species in clusters?
    '''

    # Plant agents
    class Plants:
        idd = np.arange(1, individuals+1)
        coords = np.array([[], []], dtype = float)
        sp = np.zeros((individuals,), dtype=int)
        fruits = np.array([], dtype = int)
        
    plants = Plants()
    pl_sp = plant_sp
    
    #np.random.seed(3)
    
    # Plant location
    x = env_dim*np.random.random(individuals)
    y = env_dim*np.random.random(individuals)
    plants.coords = np.append(x, y).reshape(2, individuals).transpose()
    
    # Number of fruits
    plants.fruits = gamma.rvs(2, scale = 8, size = individuals).astype(int)
    np.putmask(plants.fruits, plants.fruits > 100, 100)
    
    # Defining plant species
    if Parms.p_abund_field:
        names, prop, mat = read_abund(file_name, stage)
        order = np.argsort(prop)[::-1]
        names = names[order]
        prop = np.sort(prop)[::-1]
    else:
        names = np.arange(1, pl_sp+1).astype('S10')
        spid = np.arange(1, pl_sp+1) # species id
        s = 2 # parameter of the lognormal distribution
        prop = lognorm.pdf(spid, s, loc=0, scale=1) # proportions come from a discretized log normal distribution, with parameters (mulog = 0, sdlog = 2)

    prop = prop/sum(prop) # proportion of individuals of each species
    nind = np.rint(prop*individuals).astype(int) # number of individuals of each species
    # removing sps not present
    if Parms.p_remove_absent:    
        names = names[np.where(nind != 0)]
        nind = nind[np.where(nind != 0)]
        pl_sp = len(nind)
    while nind.sum() < individuals:
        x = np.random.choice(np.arange(pl_sp))
        nind[x] = nind[x] + 1
    while nind.sum() > individuals:
        x = np.random.choice(np.arange(pl_sp))
        nind[x] = nind[x] - 1
    prop = nind.astype(float)/nind.sum() # proportion of individuals of each species
    propacum = prop.cumsum(0) # Cumulative probability for each species
    
    # First plant
    x = np.where( plants.sp == 0 )[0]
    index_ind = np.random.choice(x)
        
    nrand = np.random.uniform()
    index_sp = np.amin(np.where( nrand < propacum ))
    
    plants.sp[index_ind] = index_sp+1
    
    # Refreshing
    nind[index_sp] = nind[index_sp]-1
    prop = nind.astype(float)/nind.sum()
    
    if prop[index_sp] > 0.0:
        prop_aux = np.delete(prop, index_sp)
        prop_aux = np.array(map(lambda x: x * (1 - rho), prop_aux))
        prop_aux = np.insert(prop_aux, index_sp, prop[index_sp] + (1 - prop[index_sp]) * rho)
        propacum = prop_aux.cumsum(0)
    else:
        propacum = prop.cumsum(0)
    
    # Other plants
    #while np.any(plants.sp == 0):
    while nind.sum() > 0:
        dists = np.array(map(distance, np.repeat(plants.coords[index_ind].reshape((1,2)), individuals, axis=0), plants.coords ))
        dist_index = np.where( (dists < dist_corr) * (plants.sp == 0) )[0]
        if dist_index.size != 0:
            index_ind = np.random.choice(dist_index)
        else:
            #dists_sort = np.sort(dists) 
            #for i in np.arange(len(dists)):
            #    indice = np.where( dists_sort[i] == dists )
            #    if plants.sp[indice] == 0:
            #        index_ind = indice
            #        break;
            x = np.where( plants.sp == 0 )[0]
            index_ind = np.random.choice(x)
            
        nrand = np.random.uniform()
        index_sp = np.amin(np.where( nrand < propacum ))
        
        plants.sp[index_ind] = index_sp+1
        
        # Refreshing
        nind[index_sp] = nind[index_sp]-1
        prop = nind.astype(float)/nind.sum()
        print nind.sum()
        
        if prop[index_sp] > 0.0:
            prop_aux = np.delete(prop, index_sp)
            prop_aux = np.array(map(lambda x: x * (1 - rho), prop_aux))
            prop_aux = np.insert(prop_aux, index_sp, prop[index_sp] + (1 - prop[index_sp]) * rho)
            propacum = prop_aux.cumsum(0)
        else:
            propacum = prop.cumsum(0)
            
    return plants.idd, plants.sp, plants.fruits, plants.coords, names

# Testing environment
#idd, sp, fruits, coords, names = create_environment(plant_sp = 30, individuals = 500, rho = 0.0, stage = 2,
#                                                    file_name = "plantas_abund_rel3.csv")
#plt.scatter(coords[:,0], coords[:,1], c = sp, s = fruits)
#sum(sp == 1)/500
#
#os.chdir("/media/windows/Users/ber/Documents/Documentos_ber/Atividades 2014/Colaboracao_Ricardo/simulacoes_interacao")
#
#f = open('teste.txt', 'w')
#
#f.write(';'.join(['id', 'x', 'y', 'sp']) + '\n')
#for i in np.arange(Parms.individuals):
#    f.write(';'.join(str(x) for x in [plants.idd[i], plants.coords[i,0], plants.coords[i,1], plants.sp[i]]) + '\n')
#
#f.close()

####################

def create_animals(env_dim = 10000, animal_sp = 5, individuals = 30, stage = 0, file_name = "", 
                   inter_mov = 7.1e-5, sl_mov = -0.002):
    '''
    This function populates a continuous landscape of dimension env_dim (m) with "individuals" 
    animals of "animal_sp" different species. 
    '''
    
    # Animal agents
    class Animals:
        idd = np.arange(1, individuals+1)
        coords = np.array([[], []], dtype = float)
        sp = np.zeros((individuals,), dtype = int)
        seed_n = []
        seed_t = []
        travel = np.zeros((individuals,), dtype = float)
        perch = np.zeros((individuals,), dtype = float)
        p_interact = np.repeat(-2, individuals) # aqui esse array nao e usado, mas ja esta pronto para considerar a ultima planta visitada por cada animal
        
    animals = Animals()
    an_sp = animal_sp
    mov = np.repeat(Parms.mov_dist, individuals)
    time = np.repeat(Parms.daily_time, individuals)
    
    # Animal location
    x = env_dim*np.random.random(individuals)
    y = env_dim*np.random.random(individuals)
    animals.coords = np.append(x, y).reshape(2, individuals).transpose()
    
    # Defining plant species
    if Parms.a_abund_field:
        names, prop, mat = read_abund(file_name, stage)
        order = np.argsort(prop)[::-1]
        if Parms.frug:
            temp = mat[:,1]/100
            temp = temp[order]
        if not Parms.mov_dist:
            mass = mat[:,0]
            mov_sp = inter_mov*np.exp(sl_mov*mass)
            mov_sp = mov_sp[order]
        names = names[order]
        prop = np.sort(prop)[::-1]
        an_sp = len(prop)
    else: 
        names = np.arange(1, an_sp+1).astype('S10')
        spid = np.arange(1, an_sp+1) # species id
        s = 5 # parameter of the lognormal distribution
        prop = lognorm.pdf(spid, s, loc=0, scale=1) # proportions come from a discretized log normal distribution, with parameters (mulog = 0, sdlog = 2)
    
    prop = prop/sum(prop) # proportion of individuals of each species
    nind = np.rint(prop*individuals).astype(int) # number of individuals of each species
    if Parms.a_remove_absent:    
        names = names[np.where(nind != 0)]  
        if Parms.frug:        
            temp = temp[np.where(nind != 0)]
        if not Parms.mov_dist:
            mov = mov[np.where(nind != 0)]
        nind = nind[np.where(nind != 0)]
        an_sp = len(nind)
    while nind.sum() < individuals:
        x = np.random.choice(np.arange(an_sp))
        nind[x] = nind[x] + 1
    while nind.sum() > individuals:
        x = np.random.choice(np.arange(an_sp))
        if nind[x] > 1:   # I am not removing any species... meybe we have to think about that!
            nind[x] = nind[x] - 1
    prop = nind.astype(float)/nind.sum() # proportion of individuals of each species
    #propacum = prop.cumsum(0) # Cumulative probability for each species
        
    #sp = animal_sp-1
    for ind in np.arange(individuals):
        #if nind[sp] <= 0:
        #    sp = sp-1
        
        #animals.sp[ind] = sp+1
        #nind[sp] = nind[sp] - 1
        animals.seed_n.append([])
        animals.seed_t.append([])

    i = 0
    sp = 0     
    while sp < an_sp:        
        if(nind[sp] > 0):         
            animals.sp[i:(i+nind[sp])] = sp+1
            if Parms.frug:            
                time[i:(i+nind[sp])] *= temp[sp]
            if not Parms.mov_dist:
                mov[i:(i+nind[sp])] = mov_sp[sp]
            i = i+nind[sp]
        sp = sp+1
        
    return animals.idd, animals.sp, animals.coords, time, animals.seed_n, animals.seed_t, animals.travel, animals.perch, animals.p_interact, mov, names
    
# Testing environment + animals
#p_idd, p_sp, p_fruits, p_coords, p_names = create_environment(plant_sp = 30, individuals = 500, rho = 0.0, file_name = "plantas_abund_rel3.csv")
#a_idd, a_sp, a_coords, an_time, seed_n, seed_t, an_trav, perch_t, p_interact, a_mov, a_names = create_animals(individuals = 220, file_name = "aves_abund_rel3.csv")
#
#fig = plt.figure()
#fig1 = fig.add_subplot(111)
#  
#aa = fig1.scatter(p_coords[:,0], p_coords[:,1], c = p_sp, s = p_fruits, alpha = 0.5)
#aa2 = fig1.scatter(a_coords[:,0], a_coords[:,1], c = a_sp, s = np.repeat(40, 30), marker = "s")
#plt.setp(aa2)

################

def initialize(stage = 0):
    '''
    This function initializes the landscape (put plant individuals in space) and provides initial
    positions to animals
    '''
    global p_idd, p_sp, p_coords, a_idd, a_sp, a_mov
    
    p_idd, p_sp, p_fruits, p_coords, p_names = create_environment(plant_sp = Parms.plant_sp, individuals = Parms.plants, rho = Parms.p_aggregation, stage = Parms.stage, file_name = Parms.p_file)
    a_idd, a_sp, a_coords, a_time, seed_n, seed_t, a_trav, perch_t, p_interact, a_mov, a_names = create_animals(animal_sp = Parms.animal_sp, individuals = Parms.animals, stage = Parms.stage, file_name = Parms.a_file)
    
    return p_fruits, a_coords, a_time, seed_n, seed_t, a_trav, perch_t, p_interact
    
def actualize(p_fruits):
    '''
    This function reinitializes seeds and seed times, travel and perch times and animal activity time
    after a day ends. Also, it calls the regrowth model for plants' fruits
    
    Depois tem que incorporar aqui 
    '''
    
    for i in np.arange(Parms.plants):
        p_fruits[i] = regrowth_model(p_fruits[i], max_fruits = 100, rate = 20)
        
    time = np.repeat(Parms.daily_time, Parms.animals)
    travel = np.zeros((Parms.animals,), dtype = float)
    perch = np.zeros((Parms.animals,), dtype = float)
    
    seed_n = []
    seed_t = []
    for ind in np.arange(Parms.animals):
        seed_n.append([])
        seed_t.append([])
    
    return p_fruits, time, seed_n, seed_t, travel, perch
    
def plot_landscape(p_fruits):
    '''
    This function plots the landscape (plant individuals)
    '''
    plt.scatter(p_coords[:,0], p_coords[:,1], c = p_sp, s = p_fruits)
    
# Ver a possibilidade de colocar uma funcao que exporta o ambiente para um arquivo
    
def plot_landscape_animals(p_fruits, a_coords):
    '''
    This function plots plant and animal individuals, in their current position
    '''
    fig = plt.figure()
    fig1 = fig.add_subplot(111)
      
    aa = fig1.scatter(p_coords[:,0], p_coords[:,1], c = p_sp, s = p_fruits, alpha = 0.5)
    aa2 = fig1.scatter(a_coords[:,0], a_coords[:,1], c = a_sp, s = np.repeat(40, 30), marker = "s")
    plt.setp(aa2)
    
def save_dynamics_step(a_coords, pl_fruits, day, time, pshow = True, ashow = True, intshow = True, 
                       intrareshow = True, intbinshow = True, intbinrareshow = True, name = 'default'):
    '''
    This function saves in a txt file the positions, identities and properties of plants and animals
    '''
        
    if pshow:
        if day == 0 and time == 0.0:
            file_n = name+'_plants_'+str(day)+'_'+str(time)+'.txt'
            
            f = open(file_n, 'w')
            
            f.write(';'.join(['id', 'x', 'y', 'sp', 'fruits']) + '\n')
            for i in np.arange(Parms.plants):
                f.write(';'.join(str(x) for x in [p_idd[i], p_coords[i,0], p_coords[i,1], p_sp[i], pl_fruits[i]]) + '\n')
            f.close() 
        elif day == 0 and time == 1.0:
            file_n = name+'_plants_fruits.txt'
            
            f = open(file_n, 'w')
            
            #f.write(';'.join(str(x) for x in [day, time]))
            f.write(str(day)+'_'+str(time)+';')
            for i in np.arange(Parms.plants):
                if i < (Parms.plants-1):
                    f.write(str(pl_fruits[i])+';')
                else:
                    f.write(str(pl_fruits[i])+'\n')
            f.close()
        else:
            file_n = name+'_plants_fruits.txt'
            
            f = open(file_n, 'a')
            
            f.write(str(day)+'_'+str(time)+';')
            for i in np.arange(Parms.plants):
                if i < (Parms.plants-1):
                    f.write(str(pl_fruits[i])+';')
                else:
                    f.write(str(pl_fruits[i])+'\n')
            f.write('\n')
            f.close()
            
    if ashow:
        file_n = name+'_animals.txt'
        if day == 0 and time == 0.0:
            f = open(file_n, 'w')
            
            f.write(';'.join(['day', 'time', 'id', 'x', 'y', 'sp']) + '\n')
            for i in np.arange(Parms.animals):
                f.write(';'.join(str(x) for x in [day, time, a_idd[i], a_coords[i,0], a_coords[i,1], a_sp[i]]) + '\n')
            f.close() 
        else:
            f = open(file_n, 'a')
            
            #f.write(';'.join(['day', 'time', 'id', 'x', 'y', 'sp']) + '\n')
            for i in np.arange(Parms.animals):
                f.write(';'.join(str(x) for x in [day, time, a_idd[i], a_coords[i,0], a_coords[i,1], a_sp[i]]) + '\n')
            f.close()
            
    if intshow:
        file_n = name+'_interactions.txt'
        if day == 0 and time == 0.0:
            f = open(file_n, 'w')
            
            f.write('day;time;animal_sp;')
            f.write(';'.join(str(x) for x in (np.arange(Parms.plant_sp)+1)) + '\n')
            for i in np.arange(Parms.animal_sp):
                f.write(str(day)+';'+str(time)+';'+str(i+1)+';')
                for j in np.arange(Parms.plant_sp):
                    if j < (Parms.plant_sp-1):
                        f.write(str(interaction_matrix[i,j])+';')
                    else:
                        f.write(str(interaction_matrix[i,j])+'\n')
            f.close() 
        else:
            f = open(file_n, 'a')
            
            #f.write('day;time;animal_sp;')
            #f.write(';'.join(str(x) for x in (np.arange(Parms.plant_sp)+1)) + '\n')
            for i in np.arange(Parms.animal_sp):
                f.write(str(day)+';'+str(time)+';'+str(i+1)+';')
                for j in np.arange(Parms.plant_sp):
                    if j < (Parms.plant_sp-1):
                        f.write(str(interaction_matrix[i,j])+';')
                    else:
                        f.write(str(interaction_matrix[i,j])+'\n')
            f.close()

    if intbinshow:
        file_n = name+'_int_bin.txt'
        if day == 0 and time == 0.0:
            f = open(file_n, 'w')
            
            f.write('day;time;animal_sp;')
            f.write(';'.join(str(x) for x in (np.arange(Parms.plant_sp)+1)) + '\n')
            for i in np.arange(Parms.animal_sp):
                f.write(str(day)+';'+str(time)+';'+str(i+1)+';')
                for j in np.arange(Parms.plant_sp):
                    if j < (Parms.plant_sp-1):
                        f.write(str(binary_matrix[i,j])+';')
                    else:
                        f.write(str(binary_matrix[i,j])+'\n')
            f.close() 
        else:
            f = open(file_n, 'a')
            
            #f.write('day;time;animal_sp;')
            #f.write(';'.join(str(x) for x in (np.arange(Parms.plant_sp)+1)) + '\n')
            for i in np.arange(Parms.animal_sp):
                f.write(str(day)+';'+str(time)+';'+str(i+1)+';')
                for j in np.arange(Parms.plant_sp):
                    if j < (Parms.plant_sp-1):
                        f.write(str(binary_matrix[i,j])+';')
                    else:
                        f.write(str(binary_matrix[i,j])+'\n')
            f.close()
            
    if intrareshow:
        file_n = name+'_interactions_raref.txt'
        if day == 0 and time == 0.0:
            f = open(file_n, 'w')
            
            f.write('day;time;animal_sp;')
            f.write(';'.join(str(x) for x in (np.arange(Parms.plant_sp)+1)) + '\n')
            for i in np.arange(Parms.animal_sp):
                f.write(str(day)+';'+str(time)+';'+str(i+1)+';')
                for j in np.arange(Parms.plant_sp):
                    if j < (Parms.plant_sp-1):
                        f.write(str(int_matrix_rarefied[i,j])+';')
                    else:
                        f.write(str(int_matrix_rarefied[i,j])+'\n')
            f.close() 
        else:
            f = open(file_n, 'a')
            
            #f.write('day;time;animal_sp;')
            #f.write(';'.join(str(x) for x in (np.arange(Parms.plant_sp)+1)) + '\n')
            for i in np.arange(Parms.animal_sp):
                f.write(str(day)+';'+str(time)+';'+str(i+1)+';')
                for j in np.arange(Parms.plant_sp):
                    if j < (Parms.plant_sp-1):
                        f.write(str(int_matrix_rarefied[i,j])+';')
                    else:
                        f.write(str(int_matrix_rarefied[i,j])+'\n')
            f.close()
            
    if intbinrareshow:
        file_n = name+'_int_bin_raref.txt'
        if day == 0 and time == 0.0:
            f = open(file_n, 'w')
            
            f.write('day;time;animal_sp;')
            f.write(';'.join(str(x) for x in (np.arange(Parms.plant_sp)+1)) + '\n')
            for i in np.arange(Parms.animal_sp):
                f.write(str(day)+';'+str(time)+';'+str(i+1)+';')
                for j in np.arange(Parms.plant_sp):
                    if j < (Parms.plant_sp-1):
                        f.write(str(bin_matrix_rarefied[i,j])+';')
                    else:
                        f.write(str(bin_matrix_rarefied[i,j])+'\n')
            f.close() 
        else:
            f = open(file_n, 'a')
            
            #f.write('day;time;animal_sp;')
            #f.write(';'.join(str(x) for x in (np.arange(Parms.plant_sp)+1)) + '\n')
            for i in np.arange(Parms.animal_sp):
                f.write(str(day)+';'+str(time)+';'+str(i+1)+';')
                for j in np.arange(Parms.plant_sp):
                    if j < (Parms.plant_sp-1):
                        f.write(str(bin_matrix_rarefied[i,j])+';')
                    else:
                        f.write(str(bin_matrix_rarefied[i,j])+'\n')
            f.close()
    

    #    f = open('animals_'+name+'.txt', 'w')
    #    
    #    f.write(';'.join(['id', 'x', 'y', 'sp']) + '\n')
    #    for i in np.arange(Parms.individuals):
    #        f.write(';'.join(str(x) for x in [plants.idd[i], plants.coords[i,0], plants.coords[i,1], plants.sp[i]]) + '\n')
    #    f.close()

# Using
#pl_fruits, an_coords, an_time, seed_n, seed_t, an_trav, perch_t, p_interact = initialize()
#
#plot_landscape(pl_fruits)
#plot_landscape_animals(pl_fruits, an_coords)
 
################

def movement_step(a_coords, p_fruits, fruit_c = 0.001, dist_c = 0.00005):
    '''
    O que fazer se o animal so fica nas mesmas plantas?
    '''
    attr_fruit = np.tanh(fruit_c * p_fruits ** 2)
    dists = np.array(map(distance, np.repeat(a_coords.reshape((1,2)), len(p_idd), axis=0), p_coords))
    np.putmask(dists, dists == 0., 100000)
    attr_dist = np.tanh(-dist_c * dists ** 2) + 1
    
    attr = attr_fruit * attr_dist
    if np.all(attr == 0.0):
        index = np.where( np.amin(dists) == dists )
    else:
        attr_prob = attr/sum(attr)
        attr_cumprob = attr_prob.cumsum(0)
    
        nrand = np.random.uniform()
        index = np.where( nrand < attr_cumprob )
        
    plant_index = np.amin(index)
    
    x = p_coords[plant_index][0]
    y = p_coords[plant_index][1]
    
    dist = distance(a_coords, np.array([x, y]))
    
    return np.array([x, y]), plant_index, dist
    
# Using
#ind = 0
#p_interact = None
##
#an_coords[ind], p_interact = movement_step(an_coords[ind], pl_fruits, last = p_interact)
#an_coords[ind], p_interact, dis = movement_step(an_coords[ind], pl_fruits)
#print an_coords[ind], p_interact
    
def perch_time(nind, time_left, shape = 4, scale = 1.25):
    '''
    
    '''
    #left = time_left    
    
    time = gamma.rvs(shape, scale = scale, size = nind)
    
    if time > time_left:
        time = np.array([time_left])
    
    #for i in np.arange(nind):   
    #    if time[i] > time_left[i]:
    #        time[i] = time_left[i]
    return time
    
def consume(fruits, filled, gutsize = 15, alfa = 10, delta = 9):
    '''
    Como sera a estrutura de sementes ingeridas e tempo de passagem? Pensar depois
    '''
    fruits = fruits.item()
    #filled = filled.item()
    
    hip = (alfa * fruits) / (delta + fruits)
    frac = (1 - filled/gutsize) * gutsize
    amount = int(min(hip, frac))
    
    time = gamma.rvs(3, scale = 9, size = 1).item()
    
    return amount, time

        
def interact(plants, p_fruits, ind, seed_n, seed_t, perch, day, time):
    '''
    Aqui preciso colocar - consumo, espera, defecar
    Defecar precisa ser no meio do caminho, mas por enquanto vamos fazer so no poleiro
    cada individuo
    
    PRECISAMOS COLOCAR PARA PODER DEFECAR NO MEIO DO CAMINHO TAMBEM
    DA PRA COLOCAR UM IF - SE ESTA NO MEIO DO CAMINHO, SO DEFECA; SE ESTA PERCHED, PODE DEFECAR E/OU COMER
    '''    
    # Defecate
    if len(seed_t) > 0:
        for j in np.arange(len(seed_n))[::-1]: 
            seed_t[j] -= perch
            if seed_t[j] <= 0:
                del seed_t[j]
                del seed_n[j]
                
    # Eat
    amount, s_time = consume(p_fruits[plants], sum(seed_n))
    if amount > 0:            
        seed_n.append(amount)
        seed_t.append(s_time)
        p_fruits[plants] -= amount
        interaction_matrix[(a_sp[ind]-1), (p_sp[plants]-1)] += amount # mudar aqui se for so a interacao
        binary_matrix[(a_sp[ind]-1), (p_sp[plants]-1)] += 1
        if day >= 4 and time < 30:
            int_matrix_rarefied[(a_sp[ind]-1), (p_sp[plants]-1)] += amount
            bin_matrix_rarefied[(a_sp[ind]-1), (p_sp[plants]-1)] += 1
        #print a_sp[ind], p_sp[plants], plants
    
    return p_fruits, seed_n, seed_t   
                            
# Using
#global interaction_matrix
#interaction_matrix = np.zeros((Parms.animal_sp, Parms.plant_sp), dtype=np.int)
#pl_fruits, an_coords, an_time, seed_n, seed_t, an_trav, perch_t, p_interact = initialize()
#
#ind = 0
#
#an_coords[ind], p_interact, dis = movement_step(an_coords[ind], pl_fruits) # coferir movimento
#print an_coords[ind], p_interact, pl_fruits[p_interact], seed_n[ind], seed_t[ind]
#
#pl_fruits, seed_n[ind], seed_t[ind] = interact(p_interact, pl_fruits, ind, seed_n[ind], seed_t[ind], 5.0) # qual o prob de intract?

def move_and_interact(pl_fruits, an_coords, an_time, seed_n, seed_t, travel, perch_t, p_interact, day, minute):
    '''
    colocar um move e uma interacao para cada individuo
    
    GARANTIR QUE USE_SE O TEMPO ATE O FIM COMO PARAMETRO PARA ANDAR E INTERAGIR OU NAO
    '''
    
    for ind in np.arange(Parms.animals):
        if an_time[ind] > 0:
            
            elapsed_t = 0.0
            while elapsed_t < Parms.time_step:
                
                if travel[ind] > (Parms.time_step - elapsed_t): # continue traveling, if already traveling
                    travel[ind] -= (Parms.time_step - elapsed_t)
                    elapsed_t = Parms.time_step
                    perch_t[ind] = 0.0
                elif perch_t[ind] > (Parms.time_step - elapsed_t): # continue perched, if already perched
                    perch_t[ind] -= (Parms.time_step - elapsed_t)
                    elapsed_t = Parms.time_step
                    travel[ind] = 0.0            
                elif ((perch_t[ind] <= (Parms.time_step - elapsed_t)) and (travel[ind] == 0.0)): # new move if it is time to move (if not stopped)
                    
                    an_coords[ind], p_interact[ind], dist = movement_step(an_coords[ind], pl_fruits, dist_c = a_mov[ind])
                
                    travel[ind] = dist/Parms.speed
                    # if there is still some time for the animal to leave the plant
                    if perch_t[ind] > 0:
                        travel[ind] -= perch_t[ind]
                    if travel[ind] > (Parms.time_step - elapsed_t):
                        travel[ind] -= (Parms.time_step - elapsed_t)
                        elapsed_t = Parms.time_step
                    else:
                        elapsed_t += (perch_t[ind] + travel[ind])
                    perch_t[ind] = 0.0 
                elif ((perch_t[ind] == 0.0) and (travel[ind] < (Parms.time_step - elapsed_t))): # continue moving and stoping
                    elapsed_t += (travel[ind])
                    
                for i in range(len(seed_t[ind])):
                    seed_t[ind][i] -= elapsed_t
                    
                # new interaction if the animal gets to the new plant in this step
                if( ((perch_t[ind] + travel[ind]) <= Parms.time_step) and (elapsed_t < Parms.time_step) ):
                    
                    travel[ind] = 0.0
                    new_perch_t = perch_time(1, an_time[ind])[0]
                    if(new_perch_t < (Parms.time_step - elapsed_t)):                        
                        pl_fruits, seed_n[ind], seed_t[ind] = interact(p_interact[ind], pl_fruits, ind, seed_n[ind], seed_t[ind], new_perch_t, day, minute)
                        elapsed_t += new_perch_t
                        perch_t[ind] = 0.0
                    else:
                        pl_fruits, seed_n[ind], seed_t[ind] = interact(p_interact[ind], pl_fruits, ind, seed_n[ind], seed_t[ind], (Parms.time_step - elapsed_t), day, minute)
                        perch_t[ind] = new_perch_t - (Parms.time_step - elapsed_t)
                        elapsed_t = Parms.time_step
                
            #plants.append(p_interact)
    
    for i in np.arange(len(an_time)):
        if(an_time[i] > 0.0):  
            an_time[i] = an_time[i] - Parms.time_step
    
    return pl_fruits, an_coords, an_time, seed_n, seed_t, travel, perch_t, p_interact

#    while( np.any(an_time > 0) ):    
#        plants = []
#        # move
#        for ind in np.arange(Parms.animals):
#            if(an_time[ind] > 0):
#                an_coords[ind], p_interact, dist = movement_step(an_coords[ind], pl_fruits)
#                
#                time_trav = dist/Parms.speed
#                an_time[ind] = an_time[ind] - time_trav if (an_time[ind] - time_trav) > 0 else 0
#                seed_t[ind] -= time_trav; seed_t[ind] = seed_t[ind].tolist()
#                plants.append(p_interact)
#         
#        # interact
#        perch = perch_time(Parms.animals, an_time)
#        
#        an_time, pl_fruits, seed_n, seed_t = interact(plants, an_time, pl_fruits, seed_n, seed_t, perch)

class Parms:
    # animal traits
    animals = 220
    animal_sp = 60
    speed = 360 # speed = 6m/s = 360m/min
    a_abund_field = 1 # 0 for theorical abundances, 1 for field abundances (need input files)
    a_remove_absent = 0 # 1 = remove absent species; 0 = maintain all species    
    mov_dist = 0.0
    frug = 1    
    
    # plants
    plants = 500 #5000
    plant_sp = 50
    p_aggregation = 0.7 # incorporar
    p_abund_field = 1 # 0 for theorical abundances, 1 for field abundances (need input files)
    p_remove_absent = 1 # 1 = remove absent species; 0 = maintain all species    
    
    # simulation
    time_step = 1.0 # 1 minute
    daily_time = 300.0 # 5h of daily activity
    time_extension = 1 #15 # 15 dias
    
    stage = 0
    p_file = "plantas_abund_rel3.csv"
    a_file = "aves_abund_rel3.csv"

def simulation(sim = ""):
    '''
    colocar 1 simulacao, por 1 mes, aqui

    '''
    
    os.chdir("/media/windows/Users/ber/Documents/Documentos_ber/Atividades 2014/Colaboracao_Ricardo/simulacoes_interacao")
    #os.chdir("C:\Users\ber\Documents\Documentos_ber\Atividades 2014\Colaboracao_Ricardo\simulacoes_interacao\results")
    
    start_time = time.time()
    
    pl_fruits, an_coords, an_time, seed_n, seed_t, an_trav, perch_t, p_interact = initialize()
    
    # redefining number of species, if there is a different number in the field data
    if Parms.p_abund_field:
        Parms.plant_sp = len(np.unique(p_sp))
    if Parms.a_abund_field:
        Parms.animal_sp = len(np.unique(a_sp))
        
    if not Parms.mov_dist:
        mov_str = 'movvar'
    else:
        mov_str = str(Parms.mov_dist)
    
    global interaction_matrix, int_matrix_rarefied, binary_matrix, bin_matrix_rarefied
    interaction_matrix = np.zeros((Parms.animal_sp, Parms.plant_sp), dtype=np.int)
    int_matrix_rarefied = np.zeros((Parms.animal_sp, Parms.plant_sp), dtype=np.int)
    binary_matrix = np.zeros((Parms.animal_sp, Parms.plant_sp), dtype=np.int)
    bin_matrix_rarefied = np.zeros((Parms.animal_sp, Parms.plant_sp), dtype=np.int)
    
    print("--- %s seconds ---" % (time.time() - start_time))    
    
    #plot_landscape_animals(pl_fruits, an_coords)        
    
    # fazer um loop ate completar o tempo diario de atividade    
    day = 0
    while day < Parms.time_extension:
        
        start_time2 = time.time()
        minutes = 0.0
        while minutes < Parms.daily_time:
            
            #save_dynamics_step(an_coords, pl_fruits, day, minutes, name = 'teste')
            pl_fruits, an_coords, an_time, seed_n, seed_t, an_trav, perch_t, p_interact = move_and_interact(pl_fruits, an_coords, an_time, seed_n, seed_t, an_trav, perch_t, p_interact, day, minutes)
            minutes += Parms.time_step
            print minutes
            
        pl_fruits, an_time, seed_n, seed_t, an_trav, perch_t = actualize(pl_fruits)
        print day
        print("--- %s seconds ---" % (time.time() - start_time2)) 
        day += 1
   
    os.chdir("/media/windows/Users/ber/Documents/Documentos_ber/Atividades 2014/Colaboracao_Ricardo/simulacoes_interacao/results")
    namearq = 'sim_'+sim+'_stg_'+str(Parms.stage)+'_'+mov_str+'_frug_'+str(Parms.frug)+'_an_'+str(Parms.animal_sp)+'_'+str(Parms.animals)+'_pl_'+str(Parms.plant_sp)+'_'+str(Parms.plants)+'_'+str(Parms.p_aggregation)
    save_dynamics_step(an_coords, pl_fruits, day, minutes, pshow = False, ashow = False, name = namearq)
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
def main():
    
    initial_t = time.time()
    
    sim = 1
    while sim <= 10:
        simulation(sim = str(sim))
        sim += 1
        
    final_t = time.time()
    
    print("--- The simulation took %s seconds ---" % (final_t - initial_t))


    

# Depois
# criar uma classe, seeds, para colocar, para cada semente:# coords iniciais, coords finais, sp, t, sp animal que dispersou, outras infos

create_animals(10000, 60, 200)


read_abund("plantas_abund_rel2.csv", 2)
