#!/usr/bin/env python3
# Name: Ioannis Anastopoulos
# Date: 11/09/2017
"""Viterbi Algorithm"""

'''Program calculates path that maximizes the (unconditional) probability Pr(x, π) over all possible paths π. Class mostLikelyHiddenPath
calcualtes max score for each state, and backtracks the graph from sink to source, picking the node that gave rise to max score, thereby 
returning the most likely path'''

import sys
from collections import defaultdict
import numpy as np
import math

class UsageError(Exception):
	def __init__(self, message):
		self.message = message

class mostLikelyHiddenPath():
	'''Class calculates most likley path over all possibile paths in a viterbi graph.Method LongestPath() calculates max score and each 
	node(state) and the node the max score came from. Method backTracking() backtracks from sink to source returning most likley path'''
	def __init__(self, emitted_sequence,emitted_alphabet,transition_alphabet, transition_matrix, emission_matrix):
		self.string_x = emitted_sequence
		self.emitted_alphabet=emitted_alphabet
		self.transition_alphabet=transition_alphabet
		self.transition_matrix=transition_matrix
		self.emission_matrix=emission_matrix
	
	def LongestPath(self):
		'''Method returns max score at each node, and the node used to calculate that score'''
		score_dict = defaultdict(dict) #saving score computations
		maxScore_predecessor_node_dict = dict() #dict of max_score: node used to calculate the max score

		for position in range(len(self.string_x)):
			for state in self.transition_matrix:
				score_dict[position][state] = -math.inf
				if position ==0: #scoring states at position 0, based on source_score=1
					source_to_state_prob = 1/len(self.transition_matrix.keys()) #prob of transitioning to either state is equal
					score = np.float128(1)*np.float128(source_to_state_prob)*np.float128(self.emission_matrix[state][self.string_x[position]]) #score at source is 1, and prob of source going to either state in position 0 is equal
					score_dict[position][state] = score

		for position in range(1,len(self.string_x)):
			for states in self.transition_matrix.keys():
				scores_list = []
				predecessor_and_prob = self.transition_matrix[states] #state current state is connected to, and the transition prob of the edge

				max_score = None
				for predecessor, prob in predecessor_and_prob.items():
					score = np.float128(score_dict[position-1][predecessor])*np.float128(self.transition_matrix[predecessor][states])*np.float128(self.emission_matrix[states][self.string_x[position]])
					scores_list.append((score))
					'''storing max score, and predecessor node used to calculate max score'''
					if max_score is None or score>max_score:
						max_score = score
						max_score_node = predecessor
					maxScore_predecessor_node_dict[max_score] = max_score_node #finding the max score from all states and edges connected to current state

				score_dict[position][states] = max(scores_list)
		return(score_dict, maxScore_predecessor_node_dict)

	def backTracking(self,longest_path_dict, maxScore_predecessor_node_dict):
		'''Method backtracks from sink to source, returning the longest path'''
		path = []			
		max_score =None
		for states_scores in longest_path_dict[len(self.string_x)-1].items(): #finding max of sink at last position of string_x to initialize backtracking
			score = states_scores[1]
			if max_score is None or score>max_score:
				max_score = score #max score
				max_score_node = states_scores[0] #state that max score at sink came from			
		path.append(max_score_node) #initializing path from sink
	
		i = sorted(range(len(self.string_x)), reverse=True) #transversing graph backwards
		for n in i:
			if longest_path_dict[n][path[-1]] in maxScore_predecessor_node_dict.keys():	
				path.append(maxScore_predecessor_node_dict[longest_path_dict[n][path[-1]]]) #appending node that gave rise to max score.
		return(''.join(str(char) for char in path[::-1]))
	
	def printPath(self):
		for key in self.transition_matrix.keys():
			if sum(self.transition_matrix[key].values())<0.998: #allows for sum of probabilities of 0.999
				raise UsageError('Transition Probability ERROR: not summing to 1')
		for key in self.emission_matrix.keys():
			if sum(self.emission_matrix[key].values())<0.998:#allows for sum of probabilities of 0.999
				raise UsageError('Emission Probability ERROR: not summing to 1')		
		longest_path_dict, maxScore_predecessor_node_dict = self.LongestPath()
		return(self.backTracking(longest_path_dict, maxScore_predecessor_node_dict))

def file_parse():
	'''Function parses file returning initial emissions and transiiton matrices
	emission and transition alphabets, and emitted sequence'''
	transition_matrix = defaultdict(dict)
	emission_matrix = defaultdict(dict)
	with sys.stdin as fn:
		lines = fn.readlines()
		emitted_sequence = lines[0].rstrip()
		emitted_alphabet=lines[2].split()
		transition_alphabet=lines[4].split()
		transition_prob = lines[7:7+len(transition_alphabet)]
		trans_probabilities=[]
		for prob in transition_prob:
			values = prob.split()
			del values[0]#deleting state from input, not assuming that row labels are in order
			trans_probabilities.append(values)
			
		for i in range(len(transition_alphabet)):
			for j in range(len(transition_alphabet)):
				transition_matrix[transition_alphabet[i]][transition_alphabet[j]]=float(trans_probabilities[i][j])

		emission_prob = lines[len(lines)-len(transition_alphabet):]

		emssn_probabilities=[]
		for prob in emission_prob:
			values=prob.split()
			del values[0]#deleting state from input, not assuming that row labels are in order
			emssn_probabilities.append(values)
			
		for i in range(len(transition_alphabet)):
			for j in range(len(emitted_alphabet)):
				emission_matrix[transition_alphabet[i]][emitted_alphabet[j]]=float(emssn_probabilities[i][j])
	return emitted_sequence,emitted_alphabet,transition_alphabet, transition_matrix, emission_matrix				
def main():
	try:
		emitted_sequence,emitted_alphabet,transition_alphabet, transition_matrix, emission_matrix=file_parse()
		res = mostLikelyHiddenPath(emitted_sequence,emitted_alphabet,transition_alphabet, transition_matrix, emission_matrix)
		print(res.printPath())
	except UsageError as error:
		print(error.message)

if __name__ == "__main__":
	main()