# -*- coding: utf-8 -*-
#
# Created by Pádraig Mac Carron

########################
#Import Libraries
import networkx as nx 
import matplotlib.pyplot as plt
import random as rd
import pylab
import numpy as np
import matplotlib.animation as animation
########################


'''
Axelrod -- Agreement Threshold
==============================

30-August-2019

This is an adjustment of the original agreement threshold
(named _bc_traits) suggested by Guillaume Deffaunt

- Pick a random couple of neighbouring agents:

- Compute proba of interaction = (nb of features with difference
    <= threshold) / total nb of features

- With this proba, choose randomly one feature where the
    agents are different with a difference below threshold.
    Then one copies the feature of the other.
 
 

'''


#########################
#Choose initial variables

#Choose grid size
Grid_x = 10
Grid_y = 10

#Choose number of features
n_features = 6

#Choose number of traits
n_traits = 3

#Agreement threshold
at = 1

#Choose total number of events
t_max = 100000

#Thickness of biggest line
thickness = 10

#Time-step multiple to display new figure at
display = 500

#line_style -- need to make this longer than n_features
line_style = ['-.',':','--','-','-']
if n_features > 4:
    line_style = ['-.',':',':','--','--','-','-','-','-','-','-']


#############
# for saving
Writer = animation.writers['ffmpeg']
writer = Writer(fps=40, metadata=dict(artist='Me'), bitrate=600)


#########################
#draw initial nodes
G = nx.grid_2d_graph(Grid_x,Grid_y)#,periodic=True)

pylab.ion()
fig = plt.figure(figsize=(8,8.5),facecolor='0.15')

pos = {i:i for i in G.nodes()}
nx.draw_networkx_nodes(G,pos,node_size=80,node_color='#348ABD')
plt.axis('off')



#############
#initial edges

traits = range(n_traits)
trait_min, trait_max = min(traits), max(traits)

features = {u:[rd.randint(trait_min,trait_max) for
               _ in range(n_features)] for u in G.nodes}


ims = []
#This loop gets the edge width for each interaction
for u,v in G.edges():
    f_u, f_v = features[u], features[v]
    common = len([i for i, j in zip(f_u, f_v) if i == j])
    #1- to get nothing for 100% match, full line for 0% match
    G.get_edge_data(u,v)['width'] = 1-(common/n_features)



#This draws the initial state for 10 runs
ims += [[plt.text(0,Grid_y+.05,'time=0',fontsize=16, color='0.15'),
         nx.draw_networkx_edges(G,pos,edge_color='0.85',
                                width=[thickness*d['width'] for u,v,d in G.edges(data=True)],
                                style=[line_style[int(np.ceil(n_features*d['width']))] for u,v,d in G.edges(data=True)],)]
        for _ in range(50)]


############
#Model


for event in range(1,t_max+1):

    #pick a random node and a random neighbour
    u = rd.choice(list(G.nodes()))
    v = rd.choice(list(nx.neighbors(G,u)))

    #get the number of common features and the index of non-common features
    f_u, f_v = features[u], features[v]
    
    if f_u != f_v:
        #non_common_features_index = [i for i,j in enumerate(zip(f_u, f_v)) if j[0] != j[1]]
        #non_common_features = [i for i, j in zip(f_u, f_v) if i != j]
        #common = (n_features - len(non_common_features_index))/n_features
        diff = np.absolute(np.array(f_u)-np.array(f_v))
        common = diff[diff<=at].size/n_features
        within = np.where((diff<=at) & (diff>0))[0]
        
        #with probability of common features change one of the non-common features of u
        prob = rd.random()
        if prob <= common and within.size>0:     
            feature_to_change = rd.choice(within)
            f_c = features[u][feature_to_change]
            features[u][feature_to_change] = f_v[feature_to_change]
                
            #Need to check u's neighbours and recalculate width
            for n in G.neighbors(u):
                f_u, f_n = features[u], features[n]
                common = len([i for i, j in zip(f_u, f_n) if i == j])
                G.get_edge_data(u,n)['width'] = 1-(common/n_features)

    if event % display == 0:
        ims += [[plt.text(0,Grid_y+.05,'time='+str(event),fontsize=16,color='w'),
                 nx.draw_networkx_edges(G,pos,edge_color='0.85',
                                        width=[thickness*d['width'] for u,v,d in G.edges(data=True)],
                                        style=[line_style[int(np.ceil(n_features*d['width']))] for u,v,d in G.edges(data=True)],)]]
    

im_ani = animation.ArtistAnimation(fig, ims, interval=1, repeat=False, )#blit=False)
#im_ani.save('Axelrod_change_by_one.mp4',writer=writer)
plt.show()

#plt.close()
