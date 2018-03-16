import sys
import networkx as nx
import matplotlib.pyplot as plt
import random

name_file = sys.argv[1]

G = nx.DiGraph()

with open(name_file) as file:
	line = file.read()

edges = line.replace("{","").split("}")
edges = edges[:-1] #take out last line

nodes = []
for i in range(len(edges)):
	edges[i] = edges[i].split(",")	
	node1 = int(edges[i][0])
	node2 = int(edges[i][1])
	if not node1 in nodes:
		nodes.append(node1)
	if not node2 in nodes:
		nodes.append(node2)

#Apparently I should not add the nodes
for i in range(len(nodes)):
	print(nodes[i])
	G.add_node(nodes[i], pos=(nodes[i],0))

for i in range(len(edges)):
	#print(edges[i][0], edges[i][1])
	if len(edges[i]) == 3:
		w = edges[i][2]
	else:
		w = round(random.random(),2)
	G.add_edge(int(edges[i][0]), int(edges[i][1]), weight=w)

print("drawing")

#pos=nx.get_node_attributes(G,'pos')
nx.draw(G,with_labels=True)
#labels = nx.get_edge_attributes(G,'weight')
#nx.draw_networkx_edge_labels(G,pos,edge_labels=labels)
plt.show()
#plt.savefig("fig.png")
