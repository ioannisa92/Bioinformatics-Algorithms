#!/usr/bin/env python3
# Name: Ioannis Anastopoulos
# Date: 11/01/2017

"""Program finds the best path in a DAG.
Program accepts DAG as an input and outputs the best path."""


import sys
from collections import defaultdict
from random import choice
import math
from copy import deepcopy

class UsageError(Exception):
	def __init__(self, message):
		self.message = message

class longestPath():
	def __init__(self):
		self.source = None
		self.sink = None

	def Graph(self):
		'''Method returns adjecency list dict and backpointer dict'''
		adjacency_dict = defaultdict(dict)
		backpointer_dict = defaultdict(dict)
		with sys.stdin as fn:
			lines = fn.readlines()
			self.source = int(lines[0].rstrip())
			self.sink = int(lines[1].rstrip())
			
			weight = None
			for line in lines[2:]:
				nodes_weights = line.rstrip().split('->')
				nodes_in = int(nodes_weights[0])
				nodes_out = nodes_weights[1].split(':')
				'''Storing highest edge weight, in thecase where there are multiple weights for the same edge'''
				if nodes_in in adjacency_dict and int(nodes_out[0]) in adjacency_dict[nodes_in]:
					current_weight = int(nodes_out[1])
					if current_weight> adjacency_dict[nodes_in][int(nodes_out[0])]:
						adjacency_dict[nodes_in][int(nodes_out[0])]=current_weight
				else:
					adjacency_dict[nodes_in][int(nodes_out[0])] = int(nodes_out[1])

			'''basically reversing adjacency_dict, retaining the weight of the edge'''
			for k,v in adjacency_dict.items():
				for i,j in v.items():
					backpointer_dict[i][k]=j
		return adjacency_dict,backpointer_dict
	
	def zeroInDegreeNode(self, adjacency_graph):
		'''method calculates in degress and out degrees, and returns nodes with zero indegrees'''
		zero_inDegree_nodes = []
		nodes = adjacency_graph.keys()
		for node in nodes:
			in_node_count = list()
			for k,v in adjacency_graph.items():
				if node in v:
					in_node_count.append(1)
			if sum(in_node_count) == 0:
				zero_inDegree_nodes.append(node)
		return zero_inDegree_nodes

	def topologicalOrder(self, zero_in_nodes, adjacency_graph, backpointer_dict):
		'''Method orders the graph in topological order'''
		topological_list = list()
		candidates = zero_in_nodes #list with zero_indegree nodes
		while len(candidates) !=0:
			candidate_node = choice(candidates)
			topological_list.append(candidate_node)
			candidates.remove(candidate_node)
			
			outgoing_nodes = adjacency_graph[candidate_node].keys() #nodes connected to the current candidate node

			for nodes in outgoing_nodes:
				del backpointer_dict[nodes][candidate_node] #deleting edge
				if len(backpointer_dict[nodes]) == 0: #checking if outnode still has incoming edges
					candidates.append(nodes)
		for k in backpointer_dict:
			if len(backpointer_dict[k]) >0:
				raise UsageError("The input grah is not a DAG")
			else:
				return topological_list

	def LongestPath(self, topological_list, adjacency_graph,backpointer_dict):
		'''Method returns max score at each node, and the node used to calculate that score'''
		score_dict = dict() #saving computations
		for node in adjacency_graph:
			score_dict[node] = -math.inf

		maxScore_predecessor_node_list = set() #list of tuples of max score, and the node used to calculate the max score
		for nodes in topological_list:
			scores_list = []
		
			predecessor_and_weights = backpointer_dict[nodes] #nodes, current node in topological order is connected to, and the weight of each edge	
			max_score = None
			for k,v in predecessor_and_weights.items():
				if k == self.source:
					score_dict[k] = 0 #setting source to 0
				score = score_dict[k]+v
				if max_score is None or score>max_score:
					max_score = score
					max_score_node = k
				scores_list.append(score_dict[k]+v) 
				maxScore_predecessor_node = (max_score, max_score_node) #(maxscore, node): tuple of max score at node, and predecessor node that was used to compute score at node
				maxScore_predecessor_node_list.add(maxScore_predecessor_node)
			score_dict[nodes] = max(scores_list) if scores_list else score_dict[nodes]

		if self.sink not in score_dict or score_dict[self.sink]==-math.inf: #taking care of case where sink is not in the adjency list
			raise UsageError('No path exists')
		return(score_dict, maxScore_predecessor_node_list)

	def backTracking(self, longest_path_dict, maxScore_predecessorNode_list):
		'''Method backtracks from sink to source, returning the longest path'''
		path = list()
		path.append(self.sink)
		while path[-1] != self.source: #stop appending to path when encountering source.
			for scores_prenodes in maxScore_predecessorNode_list:
				if longest_path_dict[path[-1]] == scores_prenodes[0]:
					if longest_path_dict[scores_prenodes[1]] != -math.inf:
						path.append(scores_prenodes[1])
		return('->'.join(str(char) for char in path[::-1]))#print path list backwards with arrows
	
	def getPath(self):
		'''method prints length of longest path, and the longest path'''
		adjacency_graph, backpointer_dict = self.Graph()
		if self.source == self.sink:
			raise UsageError('No path exists')
		backpointer_dict_copy = deepcopy(backpointer_dict) #deepcopy made because backpointer dict is destroyed in topological ordering
		adjacency_graph_copy = deepcopy(adjacency_graph)
		zeroInNodes = (self.zeroInDegreeNode(adjacency_graph))
		topological_list = (self.topologicalOrder(zeroInNodes, adjacency_graph, backpointer_dict))
		longest_path_dict, maxScore_predecessore_node_list = (self.LongestPath(topological_list, adjacency_graph_copy,backpointer_dict_copy))
		return(longest_path_dict[self.sink],self.backTracking(longest_path_dict, maxScore_predecessore_node_list))
			
def main():
	res = longestPath()
	sink_score , path = res.getPath()
	print(sink_score)
	print(path)
if __name__ == "__main__":
	main()