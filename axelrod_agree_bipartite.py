# -*- coding: utf-8 -*-
#
# Created by Pádraig Mac Carron

########################
#Import Libraries
import networkx as nx 
import matplotlib.pyplot as plt
import random as rd
import numpy as np
import pylab
from networkx.drawing.nx_agraph import graphviz_layout
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
 

 

This version shows it on a bipartite graph.
'''

#########################
#Choose initial variables

#Choose grid size
Grid_x = 10
Grid_y = 10

#Choose number of features
n_features = 7

#Choose number of traits
n_traits = 3
'''For n_traits>3 need to change the weight threshold in the 
    people space around line 256 below
    (and in the edges, subtract higher number in inequality)
'''

#Bounded-confidence threshold
bc = 1
at = bc

#Choose total number of events
t_max = 50000

#edge thickness
thickness = 1.6

#This is how much to incriment the x for each trait box (reduce for lots of traits)
xi = 0.45

#Time-step multiple to display new figure at
display = 500

############
# for saving
Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='Me'), bitrate=2400)


#########################
#draw initial nodes

#creates lattice-like graph
G = nx.grid_2d_graph(Grid_x,Grid_y,periodic=True)


#turns interactive on
pylab.ion()
fig = plt.figure(figsize=(15,9.5),facecolor='0.15')

#pos = {i:i for i in G.nodes()}
#pos = nx.shell_layout(G)
pos = nx.circular_layout(G, scale=4)
nx.draw_networkx_nodes(G,pos,node_size=40,node_color='#348ABD')
nx.draw_networkx_edges(G,pos,width=0.6,edge_color='g',style='dotted',alpha=0.8)
plt.axis('off')


##For the bipartite graph I have to create a new graph as the original 
#  has edges already that we need

H = nx.Graph()
H.add_nodes_from(G.nodes(),ntype='person')

#triangles = [(-1,-1),(1,1),(1,-1),(-1,1),(0,-1),(-1,0),(1,0),(0,1)]
###orig = [-0.07,-0.07]
colours = ['#FBC15E','#E24A33','#8EBA42','#988ED5','c','brown','pink','darkgreen']
#triangles = [(-0.2,-0.2),(0.2,0.2),(0.2,-0.2),(-0.2,0.2),(0,-0.2),(-0.,0),(0.2,0),(0,0.2)]
### Squares in straight layout
triangles = [(-3,1),#(1,1),
             (-3,0),#(1,0),
             (-3,-1),(0.8,-1),(0.8,0),(0.8,1),
             (-3,-1.7),(0.8,-1.7)]

pos_h = pos.copy()
pos_label = {}
new_nodes = []
c = 0
for f in range(n_features):
    add = 0
    name = 'f' + str(f+1) + '_q'
    for t in range(n_traits):
        node = name + str(t+1)
        H.add_node(node,ntype='att')
        new_nodes.append(node)

        tri = list(triangles[f])
        if t > 0:
            add += xi
        tri = list(triangles[f])
        tri[0] = tri[0] + add            
 
        pos_h[node] = tri
        pos_label[node] = [pos_h[node][0],pos_h[node][1]+0.2]
            
        c += 1
##        nx.draw_networkx_nodes(H, pos_h, nodelist=[node], nodesize=50,
##                               node_color=colours[f],node_shape='s',)

###Draw H's nodes and labels
##nx.draw_networkx_nodes(H,pos_h,nodelist=new_nodes,nodesize=50, node_color='y',
##                       node_shape='^',#edgecolors=,linewidths=,
##                       )

nx.draw_networkx_labels(H,pos_label, labels={u:u for u in new_nodes},size=18,font_color='w')

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
    G.get_edge_data(u,v)['width'] = thickness*(common/n_features)


#add edges to H, bipartite
for u in G:
    f_u = features[u]
    for i in range(n_features):
        ft = 'f' + str(i+1) + '_q' + str(f_u[i]+1)
        H.add_edge(u,ft,col=colours[i])
ecol = [d['col'] for u,v,d in H.edges(data=True)]

#This draws the initial state for 10 runs
ims += [[plt.text(-4,-4,'time=0',fontsize=16,color='w'),
         plt.text(-4,3,'a='+str(bc),fontsize=16,color='w'),
         #nx.draw_networkx_edges(G,pos,width=[d['width'] for u,v,d in G.edges(data=True)],edge_color='0.85'),
         nx.draw_networkx_edges(H, pos_h, width=1., edge_color=ecol,
                                style='solid',alpha=0.5)]
        for _ in range(20)]


############
#Model


for event in range(1,t_max+1):

    #pick a random node and a random neighbour
    u = rd.choice(list(G.nodes()))
    v = rd.choice(list(nx.neighbors(G,u)))

    #get the number of common features and the index of non-common features
    f_u, f_v = features[u], features[v]
    if f_u != f_v:
        #non_common_features_index = [i for i, j in enumerate(zip(f_u, f_v)) if j[0] != j[1]]
        #non_common_features = [i for i, j in zip(f_u, f_v) if i != j]
        #common = (n_features - len(non_common_features_index))/n_features
        diff = np.absolute(np.array(f_u)-np.array(f_v))
        common = diff[diff<=at].size/n_features
        within = np.where((diff<=at) & (diff>0))[0]
        
        #with probability of common features change one of the non-common features of u
        prob = rd.random()
        if prob <= common and within.size>0:
            feature_to_change = rd.choice(within)
            ft = 'f' + str(feature_to_change+1) + '_q' + str(f_u[feature_to_change]+1)
            f_c = features[u][feature_to_change]
            if np.absolute(f_c - f_v[feature_to_change]) <= bc:
                features[u][feature_to_change] = f_v[feature_to_change]

               #Need to check u's neighbours and recalculate width
                for n in G.neighbors(u):
                    f_u, f_n = features[u], features[n]
                    common = len([i for i, j in zip(f_u, f_n) if i == j])
                    G.get_edge_data(u,n)['width'] = thickness*(common/n_features)

                #H's edges need to be changed
                H.remove_edge(u,ft)
                ft = 'f' + str(feature_to_change+1) + '_q' + str(f_v[feature_to_change]+1)
                H.add_edge(u,ft,col=colours[feature_to_change])
            

    if event % display == 0:
        ecol = [d['col'] for u,v,d in H.edges(data=True)]
        ims += [[plt.text(-4,-4,'time='+str(event),fontsize=16,color='w'),
                 plt.text(-4,3,'a='+str(bc),fontsize=16,color='w'),
                 #nx.draw_networkx_edges(G,pos,width=[d['width'] for u,v,d in G.edges(data=True)]
                 #                       ,edge_color='0.85'),
                 nx.draw_networkx_edges(H, pos_h, width=1., edge_color=ecol,
                                        style='solid',alpha=0.5)
                 ]]
    

ecol = [d['col'] for u,v,d in H.edges(data=True)]
#This draws the final state for 50 runs
ims += [[plt.text(-4,-4,'time='+str(event),fontsize=16,color='w'),
         plt.text(-4,3,'a='+str(bc),fontsize=16,color='w'),
         #nx.draw_networkx_edges(G,pos,width=[d['width'] for u,v,d in G.edges(data=True)],edge_color='0.85'),
         nx.draw_networkx_edges(H, pos_h, width=1., edge_color=ecol,
                                style='solid',alpha=0.5)]
        for _ in range(50)]

#############
#Produce animation

im_ani = animation.ArtistAnimation(fig, ims, interval=1, repeat=False, )#blit=False)
im_ani.save('Axelrod_agreethres_fq_' + str(n_features) +str(n_traits) +'_a'+str(bc)+'.mp4',writer=writer, savefig_kwargs={'facecolor':'0.2'})

#cols = [ for i,u in enumerate(G)


#############
#Final State

plt.figure(figsize=(15,9.5),facecolor='0.15')
nx.draw_networkx_nodes(G,pos,node_size=40,node_color='#348ABD')
plt.axis('off')

nx.draw_networkx_labels(H,pos_label, labels={u:u for u in new_nodes},size=28,font_color='w')
nx.draw_networkx_edges(H, pos_h, width=1., edge_color=ecol,style='solid',alpha=0.5)
#plt.savefig('figs/trait_bounded_confidence_final_state' + str(n_features) +str(n_traits) +'.png',bbox_inches='tight',facecolor='0.15')

#############
#People space
weight_thres = n_features-1
g = nx.bipartite.weighted_projected_graph(H, G.nodes())
h = nx.Graph()
h.add_weighted_edges_from([(u,v,d['weight']) for u,v,d in g.edges(data=True)
                               if d['weight'] >= weight_thres])
plt.figure(figsize=(15,9),facecolor='0.15')
plt.axis('off')
pos_p = nx.nx_agraph.graphviz_layout(h)

nx.draw_networkx_nodes(h, pos_p, node_size=50, node_color='#FBC15E')
weak = [(u,v) for u,v,d in h.edges(data=True) if d['weight'] < n_features-2]# and h.has_edge(u,v)]
nx.draw_networkx_edges(h,pos_p,edgelist=weak,width=0.6,style='dotted',
                       edge_color='#30B838',alpha=0.95)
med = [(u,v) for u,v,d in h.edges(data=True) if d['weight'] == n_features-1]
nx.draw_networkx_edges(h,pos_p,edgelist=med,width=0.6,style='dashed',
                       edge_color='#348ABD',alpha=0.8)


med2 = [(u,v) for u,v,d in g.edges(data=True) if d['weight'] == n_features-2]
nx.draw_networkx_edges(h,pos_p,edgelist=med2,width=0.4,style='dotted',
                       edge_color='#E24A33',alpha=0.8)


strong = [(u,v) for u,v,d in g.edges(data=True) if d['weight'] == n_features]
nx.draw_networkx_edges(h,pos_p,edgelist=strong,width=0.8
                       ,edge_color='0.85',alpha=0.9,facecolor='0.15')

plt.savefig('trait_bc_people' + str(n_features) +str(n_traits) +'.svg',bbox_inches='tight',facecolor='0.15')
################
#Attitude space

att = nx.bipartite.weighted_projected_graph(H, new_nodes)
att = max(nx.connected_component_subgraphs(att), key=len)
weight = [d['weight'] for u,v,d in att.edges(data=True)]

pos_a = nx.spring_layout(att)
pos_lab = {u:[pos_a[u][0],pos_a[u][1]+0.03] for u in pos_a}
plt.figure(figsize=(12,9),facecolor='0.15')
plt.axis('off')
nx.draw_networkx_nodes(att, pos_a, node_size=80, node_color='#FBC15E')
nx.draw_networkx_labels(att,pos_lab,font_color='w')

low = max(weight)/3#45
high = 2*max(weight)/3#89

e_weak = [(u,v) for u,v,d in att.edges(data=True) if d['weight'] < low]
w_weak = [d['weight'] for u,v,d in att.edges(data=True) if d['weight'] < low]

e_med = [(u,v) for u,v,d in att.edges(data=True)
   if d['weight'] >=low and d['weight'] < high]
w_med = [d['weight'] for u,v,d in att.edges(data=True)
      if d['weight'] >=low and d['weight'] < high]

e_strong =[(u,v) for u,v,d in att.edges(data=True) if d['weight'] >= high]
w_strong = [d['weight'] for u,v,d in att.edges(data=True) if d['weight'] >= high]

nx.draw_networkx_edges(att, pos_a, edgelist=e_weak, edge_color='#FFB5B8',
		   style='dotted',width=np.array(w_weak)/40,alpha=0.99)
nx.draw_networkx_edges(att, pos_a, edgelist=e_med, edge_color='#8EBA42',
		   style='dashed',width=np.array(w_med)/60,alpha=0.9)
nx.draw_networkx_edges(att, pos_a, edgelist=e_strong, edge_color='0.9',
		   style='solid',width=np.array(w_strong)/70)

#plt.savefig('figs/trait_bounded_confidence_attitudes' + str(n_features) +str(n_traits) +'.png',bbox_inches='tight',facecolor='0.15')

################

plt.show()
#plt.close()

'''
--Use stochastic block model to get the number of clusters, then compare with
standard Axelrod for standard block model

'''
