#!/usr/bin/env python3
# Name: Ioannis Anastopoulos
# Date: 12/01/2017


'''Program calculates forward and backward scores of a given transition and emission matrix in clsas SoftDecoding. Method printProbability in this class, 
calculates Pr(πi = k|x) that the HMM was in state k at step i (for each state k and each step i) and returns a numpy matrix. Matrix is printed out in main.
File is parsed in function file_parse'''


import sys
from collections import defaultdict
import numpy as np
import math

class SoftDecoding():
	'''Class calculates forward and backwards scores, returing node and edge responsibility matrixes'''
	def __init__(self, emitted_sequence,emitted_alphabet, transition_alphabet,transition_matrix, emission_matrix):
		self.string_x = emitted_sequence
		self.transition_matrix=transition_matrix
		self.emission_matrix=emission_matrix
		self.transition_alphabet=transition_alphabet
		self.emission_alphabet=emitted_alphabet
	
	def Forward(self, transition_matrix, emission_matrix):
		'''Method calculates forward scores in dict of each position in the emitted sequence'''
		forward_dict = defaultdict(dict) #saving computations in dict of dicts {string_index:state: forward score}
		for i in range(len(self.string_x)):
			for states in self.transition_alphabet:
				forward_dict[i][states]=0#initialixing forward dict with zero counts
				if i==0:#initializing forward dictionary: transition prob from source to states if first position(index:0) is equal
					forward_dict[i][states]=1/len(self.transition_alphabet)*emission_matrix[states][self.string_x[i]]
				else:
					for prev_state,score in forward_dict[i-1].items():#predecessor states current state is connected to, and the score of these predecessors
						forward_dict[i][states]+=forward_dict[i-1][prev_state]*transition_matrix[prev_state][states]*emission_matrix[states][self.string_x[i]]#calculating backward scores, based on scores from predecesors
		forward_sink=sum(forward_dict[len(self.string_x)-1].values())
		return(forward_dict, forward_sink)

	def Backward(self, transition_matrix, emission_matrix):
		'''Method returns dict of sum of scores for each state at each position on string_x'''
		backward_dict = defaultdict(dict) #saving computations in dict of dicts {string_index:state: backward score}

		for state in self.transition_alphabet:
			backward_dict[max(list(range(len(self.string_x))))][state] = 1
		
		for position in sorted(list(range(len(self.string_x))), reverse=True):
			if position!=max(list(range(len(self.string_x)))):
				for states in self.transition_alphabet:#current state
					scores_list = []#list of backward calculations to be summed later for each of the current states
					for predecessor,p_score in backward_dict[position+1].items(): #predecessor states current state is connected to, and the score of these predecessors
						score = p_score*(transition_matrix[states][predecessor]*emission_matrix[predecessor][self.string_x[position+1]])
						scores_list.append((score))
					backward_dict[position][states] = sum(scores_list) #summing scores from all states and edges connected to current state
		return(backward_dict)
	
	def printProbability(self):
		'''Method print probability Pr(πi = k|x) that the HMM was in state k at step i'''
		forward, forward_sink = self.Forward(self.transition_matrix, self.emission_matrix)
		backward = self.Backward(self.transition_matrix, self.emission_matrix)
		output = np.zeros(shape=(len(self.string_x), len(self.transition_alphabet)))
		for i in range(len(self.string_x)):
			for state in self.transition_alphabet:
				output[i,self.transition_alphabet.index(state)]=forward[i][state]*backward[i][state]/forward_sink
		return(output)

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
	emitted_sequence,emitted_alphabet,transition_alphabet, transition_matrix, emission_matrix=file_parse()
	res = SoftDecoding(emitted_sequence,emitted_alphabet,transition_alphabet, transition_matrix, emission_matrix)
	### -------------Printing matrix ------------ ####
	output = (res.printProbability())
	print("\t".join(char for char in res.transition_alphabet))
	for i,e in enumerate(res.string_x):
		for index in range(len(output)):
			if index==i:
				print('\t'.join(format(num,".4f") for num in output[index]))
	### -------------Printing matrix ------------ ####
if __name__ == "__main__":
	main()