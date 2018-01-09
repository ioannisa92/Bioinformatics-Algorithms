#!/usr/bin/env python3
# Name: Ioannis Anastopoulos
# Date: 12/01/2017

"""CAUTION: python3.6 is recommended with this program - there is an issue with the sorting of dictionaries
in earlier versions of python, which python3.6 takes care of"""

'''Program conducts Baum-Welch learning. Class BaumWelchLearning accepts initial transition and emission matrices
and passes them to the SoftDecoding class, which calculates Pr(πi = k|x) Pr(πi = l, πi+1=k|x), in other words Pi* and Pi**. 
Class BaumWelchLearning converts Pi* and PI** to dict of dicts and calculates new transition and emission matrices for any
number of iterations the file specifies''' 

import sys
from collections import defaultdict
import numpy as np
import itertools

class SoftDecoding():
	'''Class calculates forward and backwards scores, returing node and edge responsibility matrixes'''
	def __init__(self, emitted_sequence, transition_alphabet,transition_matrix, emission_matrix):
		self.string_x = emitted_sequence
		self.transition_matrix=transition_matrix
		self.emission_matrix=emission_matrix
		self.transition_alphabet=transition_alphabet

	def Forward(self):
		'''Method calculates forward scores in dict of each position in the emitted sequence'''
		forward_dict = defaultdict(dict) #saving computations in dict of dicts {string_index:state: forward score}

		for i in range(len(self.string_x)):
			for states in self.transition_alphabet:
				forward_dict[i][states]=0 #initialixing forward dict with zero counts
				if i==0: #initializing forward dictionary: transition prob from source to states if first position(index:0) is equal
					forward_dict[i][states]=1/len(self.transition_alphabet)*self.emission_matrix[states][self.string_x[i]]
				else:
					for prev_state,score in forward_dict[i-1].items():#predecessor states current state is connected to, and the score of these predecessors
						forward_dict[i][states]+=forward_dict[i-1][prev_state]*self.transition_matrix[prev_state][states]*self.emission_matrix[states][self.string_x[i]] #calculating backward scores, based on scores from predecesors
		forward_sink=sum(forward_dict[len(self.string_x)-1].values())
		return(forward_dict, forward_sink)

	def Backward(self):
		'''Method returns dict of sum of scores for each state at each position on string_x'''
		backward_dict = defaultdict(dict) #saving computations in dict of dicts {string_index:state: backward score}

		for state in self.transition_matrix:
			backward_dict[max(list(range(len(self.string_x))))][state] = 1 #initializing backward dictionary with a score of 1 for the states in the last position
		
		for position in sorted(list(range(len(self.string_x))), reverse=True):
			if position!=max(list(range(len(self.string_x)))):
				for states in self.transition_matrix: #current state 
					scores_list = [] #list of backward calculations to be summed later for each of the current states
					for predecessor,p_score in backward_dict[position+1].items(): #predecessor states current state is connected to, and the score of these predecessors
						score = p_score*(self.transition_matrix[states][predecessor]*self.emission_matrix[predecessor][self.string_x[position+1]])
						scores_list.append((score))
					backward_dict[position][states] = sum(scores_list) #summing scores from all states and edges connected to current state
		return(backward_dict)
	
	def responsibilityMatrices(self):
		'''Method retuns  probability Pr(πi = k|x) and Pr(πi = l, πi+1=k|x) that the HMM was in state k at step i'''
		forward, forward_sink = self.Forward()
		backward = self.Backward()
		
		# Pr(πi = k|x)
		pi_star = np.zeros(shape=(len(self.string_x), len(self.transition_alphabet)))
		for i in range(len(self.string_x)):
			for state in self.transition_alphabet:
				pi_star[i,self.transition_alphabet.index(state)]=forward[i][state]*backward[i][state]/forward_sink
		#creating list of all possible transitions
		all_state_transitions= [char for char in itertools.product(self.transition_alphabet, repeat=len(self.transition_alphabet))]
		
		#Pr(πi = l, πi+1=k|x)
		pi_star_star = np.zeros(shape=(len(all_state_transitions), len(self.string_x)-1))
		for i in range(len(self.string_x)-1):
			for state_trans in all_state_transitions:
				forward_prob=forward[i][state_trans[0]] #forward score at i
				backward_prob=backward[i+1][state_trans[1]] #backward score at i+1
				edge_weight=self.transition_matrix[state_trans[0]][state_trans[1]]*self.emission_matrix[state_trans[1]][self.string_x[i+1]] #edge weight from l,i, to k,i+1
				pi_star_star[all_state_transitions.index(state_trans), i]=forward_prob*backward_prob*edge_weight/forward_sink
		return(pi_star,pi_star_star, all_state_transitions)

class BaumWelchLearning():
	'''Class accepts initial transition and emission matrices from input file, passes them in SoftDecoding Class where
	Pr(πi = k|x) and Pr(πi = l, πi+1=k|x) are calculated. Class converts Pi* and Pi** to dict of dicts with matrixBuild method'''
	def __init__(self, iterations, emitted_sequence, emitted_alphabet, transition_alphabet, transition_matrix, emission_matrix):
		self.iterations=iterations
		self.emitted_sequence=emitted_sequence
		self.emitted_alphabet=emitted_alphabet
		self.transition_alphabet=transition_alphabet
		self.transition_matrix=transition_matrix
		self.emission_matrix=emission_matrix

	def matrixBuild(self, pi_star, pi_star_star, all_state_transitions):
		'''Method converts numpy matrices of pi_star, and pi_star_star to dict of dicts for transition and emission probabilities'''
		new_transition_matrix= defaultdict(dict)
		new_emission_matrix=defaultdict(dict)
	
		for states in self.transition_alphabet:
			for char in self.emitted_alphabet:#eliminating duplicate emission symbols in input
				new_emission_matrix[states][char]=0

		for i in range(len(pi_star)):
			for states in new_emission_matrix.keys():
				if self.emitted_sequence[i] in new_emission_matrix[states]:
					new_emission_matrix[states][self.emitted_sequence[i]]+=pi_star[i,self.transition_alphabet.index(states)] #summing prob of each char emission
		for key in new_emission_matrix:
			total=sum(new_emission_matrix[key].values())
			for subkey in new_emission_matrix[key]:
				new_emission_matrix[key][subkey]/=total #normalization of sums of each char probability

		for i,state_trans in enumerate(all_state_transitions):
			new_transition_matrix[state_trans[0]][state_trans[1]]=pi_star_star[i].sum()#summing prob of each state trasition probabilities
		for key in new_transition_matrix:
			total=sum(new_transition_matrix[key].values())
			for subkey in new_transition_matrix[key]:
				new_transition_matrix[key][subkey]=new_transition_matrix[key][subkey]/total#normalizing sum from above iteration.
		return(new_transition_matrix,new_emission_matrix )
	
	def Iterations(self):
		'''method iterates a given amount of times, creating new Pi* and Pi**, and new emission and transition matrices'''
		current_transition_matrix=self.transition_matrix #initial transition matrix
		current_emission_matrix=self.emission_matrix #initial emission matric
		for i in range(self.iterations):
			pi_star, pi_star_star, all_state_transitions = SoftDecoding(self.emitted_sequence, self.transition_alphabet,current_transition_matrix, current_emission_matrix).responsibilityMatrices()
			n_transition_matrix, n_emission_matrix=self.matrixBuild(pi_star, pi_star_star, all_state_transitions)
		
			current_transition_matrix=n_transition_matrix #new transition matrix after each iteration
			current_emission_matrix=n_emission_matrix #new emission matric after each iteration
		return(current_transition_matrix, current_emission_matrix)

def file_parse():
	'''Function parses file returning initial emissions and transiiton matrices
	emission and transition alphabets, numbers of iterations, and emitted sequence'''
	transition_matrix = defaultdict(dict)
	emission_matrix = defaultdict(dict)
	with sys.stdin as fn:
		lines = fn.readlines()
		iterations = int(lines[0].rstrip())
		emitted_sequence=lines[2].rstrip()
		emitted_alphabet=lines[4].split()
		transition_alphabet=lines[6].split()

		transition_prob = lines[9:9+len(transition_alphabet)]
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
	return iterations,emitted_sequence,emitted_alphabet,transition_alphabet, transition_matrix, emission_matrix

def main():
	iterations,emitted_sequence,emitted_alphabet,transition_alphabet, transition_matrix, emission_matrix=file_parse()
	res=BaumWelchLearning(iterations,emitted_sequence,emitted_alphabet,transition_alphabet, transition_matrix, emission_matrix)
	transition_matrix, emission_matrix= res.Iterations()
	print(transition_matrix.keys())
	### -------------Printing matrices ------------ ####
	print('\t'+'\t'.join(char for char in res.transition_alphabet))
	for i,e in enumerate(res.transition_alphabet):
		for index in range(len(transition_matrix)):
			if index==i:
				print(e+'\t'+'\t'.join(format(num,".3f") for num in transition_matrix[e].values()))	
	print('--------')
	print('\t'+'\t'.join(char for char in res.emitted_alphabet))
	for i,e in enumerate(res.transition_alphabet):
		for index in range(len(transition_matrix)):
			if index==i:
				print(e+'\t'+'\t'.join(format(num,".3f") for num in emission_matrix[e].values()))
	### -------------Printing matrices ------------ ####

if __name__ == "__main__":
	main()