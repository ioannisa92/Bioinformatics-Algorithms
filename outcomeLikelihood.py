#!/usr/bin/env python3
# Name: Ioannis Anastopoulos
# Date: 11/09/2017

'''Program calculates Pr(x) over all possible paths π. Class mostLikelyHiddenPath
calculates the sum of scores from all predecessors of each state. The sum of the predecessors of the sink
returns the Pr(x) that the HMM emits the given string'''
import sys
from collections import defaultdict
import numpy as np
import math

class UsageError(Exception):
	def __init__(self, message):
		self.message = message

class mostLikelyHiddenPath():
	'''class calculates Pr(x) over all possibile paths in a viterbi graph.Method LongestPath() calculates max score and each 
	node(state) and the node the max score came from. Method backTracking() backtracks from sink to source returning most likley path'''
	def __init__(self, string_x, states, symbols, transition_matrix, emission_matrix):
		self.string_x = string_x
		self.states=states
		self.symbols=symbols
		self.transition_matrix=transition_matrix
		self.emission_matrix=emission_matrix

	def LongestPath(self):
		'''Method returns dict of sum of scores for each state at each position on string_x'''
		score_dict = defaultdict(dict) #saving computations

		for position in range(len(self.string_x)):
			for state in self.transition_matrix:
				score_dict[position][state] = -math.inf
				if position ==0:#scoring states at position 0, based on source_score=1
					source_to_state_prob = 1/len(self.transition_matrix.keys()) #prob of transitioning to either state is equal
					score = (1)*source_to_state_prob*self.emission_matrix[state][self.string_x[position]] #score at source is 1, and prob of source going to either state in position 0 is equal
					score_dict[position][state] = score

		for position in range(1,len(self.string_x)):
			for states in self.transition_matrix.keys():
				scores_list = []
				predecessor_and_prob = self.transition_matrix[states] #state current state is connected to, and the transition prob of the edge
				for predecessor, prob in predecessor_and_prob.items():
					score = score_dict[position-1][predecessor]*(self.transition_matrix[predecessor][states]*self.emission_matrix[states][self.string_x[position]])
					scores_list.append((score))

				score_dict[position][states] = sum(scores_list) #summing scores from all states and edged connected to current state
		return(score_dict)

	def probabilityOfEmission(self, longest_path_dict):
		'''Method calculates Pr(x)'''
		sink_score = sum(longest_path_dict[len(self.string_x)-1].values())
		return(sink_score) #sum of state scores that lead to sink, is Pr(x)

	def printProbability(self):
		'''Method print probability result'''
		for key in self.transition_matrix.keys():
			if sum(self.transition_matrix[key].values())<0.998:#allows for sum of probabilities of 0.999
				raise UsageError('Transition Probability ERROR: not summing to 1')
		for key in self.emission_matrix.keys():
			if sum(self.emission_matrix[key].values())<0.998:#allows for sum of probabilities of 0.999
				raise UsageError('Emission Probability ERROR: not summing to 1')
		longest_path_dict = self.LongestPath()
		return(self.probabilityOfEmission(longest_path_dict))

def file_parse():
	transition_matrix = defaultdict(dict)
	emission_matrix = defaultdict(dict)
	with sys.stdin as fn:
		lines = fn.readlines()
		string_x = lines[0].rstrip()

		states = lines[4].split()
		transition_prob = lines[7:7+len(states)]

		for prob in transition_prob:
			values = prob.split()
			for i in range(len(states)):
				transition_matrix[values[0]][states[i]] = float(values[i+1])

		symbols = lines[2].split()
		emission_prob = lines[len(lines)-len(states):]

		for emissions in emission_prob:
			prob = emissions.rstrip().split()
			for i in range(len(symbols)):
				emission_matrix[prob[0]][symbols[i]] = float(prob[i+1])
	return string_x, states, symbols, transition_matrix, emission_matrix

def main():
	try:
		string_x, states, symbols, transition_matrix, emission_matrix=file_parse()
		res = mostLikelyHiddenPath(string_x, states, symbols, transition_matrix, emission_matrix)
		print(res.printProbability())
	except UsageError as error:
		print(error.message)


if __name__ == "__main__":
	main()