#!/usr/bin/env python3
# Name: Ioannis Anastopoulos
# Date: 11/29/2017
import sys
import numpy as np
import math
from collections import defaultdict

class ViterbiAlgorithm():
	'''Class uses the viterbi algorithm and returns the best hidden state path.
	Method LongestPath() calculates max score of each state, and method backTracking(), backtracks through the 
	hidden states and returns the best hidden path'''
	def __init__(self, string_x, transition_matrix, emission_matrix):
		self.string_x=string_x
		self.transition_matrix=transition_matrix
		self.emission_matrix=emission_matrix
	def LongestPath(self):
		'''Method returns max score at each node, and the node used to calculate that score'''
		score_dict = defaultdict(dict) #saving computations in dict of dicts {string_index:state: forward score}
		maxScore_predecessor_node_dict = dict() #dict of {max_score: node} used to calculate the max score

		for position in range(len(self.string_x)):
			for state in self.transition_matrix:
				score_dict[position][state] = -math.inf
				if position ==0: #scoring states at position 0, based on source_score=1
					source_to_state_prob = 1/len(self.transition_matrix.keys()) #prob of transitioning to either state is equal
					score = 1*(source_to_state_prob*self.emission_matrix[state][self.string_x[position]]) #score at source is 1, and prob of source going to either state in position 0 is equal
					score_dict[position][state] = score

		for position in range(1,len(self.string_x)):
			for states in self.transition_matrix.keys():
				scores_list = [] #saving calculations in order to calculate max later
				predecessor_and_prob = self.transition_matrix[states] #state current state is connected to, and the transition prob of the edge

				max_score = None
				for predecessor, prob in predecessor_and_prob.items():
					score = score_dict[position-1][predecessor]*(self.transition_matrix[predecessor][states]*self.emission_matrix[states][self.string_x[position]])
					scores_list.append((score))
					'''storing max score, and predecessor node used to calculate max score'''
					if max_score is None or score>max_score:
						max_score = score 
						max_score_node = predecessor #precessor nore the max score came from
					maxScore_predecessor_node_dict[max_score] = max_score_node #finding the max score from all states and edges connected to current state

				score_dict[position][states] = max(scores_list)
		return(score_dict, maxScore_predecessor_node_dict)

	def backTracking(self,longest_path_dict, maxScore_predecessor_node_dict):
		'''Method backtracks from sink to source, returning the longest path'''
		path = []			
		max_score =None
		for states_scores in longest_path_dict[len(self.string_x)-1].items(): #finding max of sink at last position of string_x
			score = states_scores[1]
			if max_score is None or score>max_score:
				max_score = score
				max_score_node = states_scores[0] #state that max score at sink came from			
		path.append(max_score_node)
	
		i = sorted(range(len(self.string_x)), reverse=True) #transversing graph backwards
		for n in i:
			if longest_path_dict[n][path[-1]] in maxScore_predecessor_node_dict.keys():	
				path.append(maxScore_predecessor_node_dict[longest_path_dict[n][path[-1]]])
		return(''.join(str(char) for char in path[::-1]))
	
	def printPath(self):
		'''Method returns path from method backTracking()'''
		longest_path_dict, maxScore_predecessor_node_dict = self.LongestPath()
		return(self.backTracking(longest_path_dict, maxScore_predecessor_node_dict))	

class HMMParameterEstimation():
	'''Class calculates a matrix of transition probabilities and a matrix of emission probabilities 
	 that maximize Pr(x,π) over all possible matrices of transition and emission probabilities. Method transitionCounts()
	 counts number of transitions from the given path, and method emissionCounts() calculates emissions from the given string'''
	def __init__(self, emitted_sequence, emission_alphabet, new_path, transition_alphabet):
		self.emitted_sequence=emitted_sequence
		self.emission_alphabet=emission_alphabet
		self.path=new_path
		self.path_alphabet=transition_alphabet
	
	def transitionCounts(self):
		'''Transition counts from given path'''
		transition_counts=defaultdict(dict)
		for i in range(len(self.path)-1):
			for char in self.path_alphabet:
				if self.path[i]==char: #current state on the path
					if self.path[i] in transition_counts and self.path[i+1] in transition_counts[self.path[i]]:
						transition_counts[self.path[i]][self.path[i+1]]+=1 #state next to the currect state, meaning transition - adding 1 for that transition
					else:
						transition_counts[self.path[i]][self.path[i+1]]=1 #if state next to the currect state, hasnt been encountered before, equal the transition to 1
				if char not in transition_counts:
					transition_counts[char]={chars:0 for chars in self.path_alphabet} #zeros are added for states never encounterd
		for key in transition_counts:
			for char in self.path_alphabet:
				if char not in transition_counts[key]:
					transition_counts[key][char]=0 #if state has never been seen before equal it to 0
		return(transition_counts)
	def emissionCounts(self):
		'''Emission counts from given string'''
		emission_counts=defaultdict(dict)

		for i in range(len(self.emitted_sequence)):
			if self.path[i] not in emission_counts:#current state on the path
				emission_counts[self.path[i]]={} #dict of {state:emission:count}
			if self.emitted_sequence[i] not in emission_counts[self.path[i]]:
				emission_counts[self.path[i]][self.emitted_sequence[i]]=1 #if emission not encountered in current state, equal it to 1
			else:
				emission_counts[self.path[i]][self.emitted_sequence[i]]+=1 #if emission has been encountered before, add 1 to count
		for key in emission_counts:
			for symbol in self.emission_alphabet:
				if symbol not in emission_counts[key]:
					emission_counts[key][symbol]=0 #if emission has never been seen before equal it to 0
		for char in self.path_alphabet:
			for symbol in self.emission_alphabet:
				if char not in emission_counts:
					emission_counts[char]={chars:0 for chars in self.emission_alphabet}#if emission for current state has never been seen before equal it to 0
		return(emission_counts)
	
	def MatrixBuild(self,transition_counts,emission_counts):
		'''method accepts transition and emission counts, and returns normalized transition and emission dict of dicts'''
		tranitions_prob_dict=defaultdict(dict)
		emissions_prob_dict=defaultdict(dict)
		if len(transition_counts)==0:
			for i in range(len(self.path_alphabet)):
				for j in range(len(self.path_alphabet)):
					tranitions_prob_dict[self.path_alphabet[i]][self.path_alphabet[j]]=np.float64(1/len(set(self.path_alphabet)))
		for key in transition_counts:
			for subkey in transition_counts[key]:
				if sum(transition_counts[key].values())==0:
					tranitions_prob_dict[key][subkey]=np.float64(transition_counts[key][subkey])+1/len(transition_counts[key]) #equal transition probability if all coutns are 0
				elif key==subkey and len(transition_counts)==1:
					for i in range(len(self.path_alphabet)):
						for j in range(len(self.path_alphabet)):
							tranitions_prob_dict[self.path_alphabet[i]][self.path_alphabet[j]]=np.float64(transition_counts[key][subkey])/sum(transition_counts[key].values())
				else:
					tranitions_prob_dict[key][subkey]=np.float64(transition_counts[key][subkey])/sum(transition_counts[key].values()) #noramlization bu sum of coutns
		if len(emission_counts)==1 and len(self.path_alphabet)!=1:
			for i in range(len(self.path_alphabet)):
				for j in range(len(self.emission_alphabet)):
					emission_matrix[self.path_alphabet[i]][self.emission_alphabet[j]]=np.float64(emission_counts[self.path_alphabet[i]][self.emission_alphabet[j]])/sum(emission_counts[self.path_alphabet[i]].values())
		for key in emission_counts:
			for subkey in emission_counts[key]:
				if sum(emission_counts[key].values())==0:
					emissions_prob_dict[key][subkey]=emission_counts[key][subkey]+1/len(emission_counts[key]) #equal emission probability if all coutns are 0
				else:
					emissions_prob_dict[key][subkey]=emission_counts[key][subkey]/sum(emission_counts[key].values()) #normalization by sum of coutns
		return(tranitions_prob_dict,emissions_prob_dict)

	def wrap(self):
		'''Method wraps transitionCounts, emissionCounts, and MatrixBuild methods to return transition and emission dict of dicts'''
		transition_counts=self.transitionCounts()
		emission_counts=self.emissionCounts()
		return(self.MatrixBuild(transition_counts,emission_counts))

class ViterbiLearning():
	'''Class applies Viterbi Learning algorithm over all possible transition and emission matrices and over all hidden paths π.
	iteration() method uses ViterbiAlgorithm class to get a new path after every iteration, and HMMParameterEstimation class to get 
	a new set of transition and emission parameters after each iteration'''
	def __init__(self, iterations,emitted_sequence,emitted_alphabet,transition_alphabet, transition_matrix, emission_matrix):
		self.iterations=iterations
		self.emitted_sequence=emitted_sequence
		self.emitted_alphabet=emitted_alphabet
		self.transition_alphabet=transition_alphabet
		self.transition_matrix=transition_matrix
		self.emission_matrix=emission_matrix
		
	def iteration(self, transition_matrix, emission_matrix):
		current_transition_parameters=self.transition_matrix
		current_emission_parameters=self.emission_matrix

		for i in range(self.iterations):
			new_path=ViterbiAlgorithm(self.emitted_sequence,current_transition_parameters,current_emission_parameters).printPath()
			new_parameters=HMMParameterEstimation(self.emitted_sequence, self.emitted_alphabet, new_path, self.transition_alphabet).wrap()
			new_transition_parameters, new_emission_parameters= new_parameters
			
			current_path=new_path
			current_transition_parameters=new_transition_parameters
			current_emission_parameters=new_emission_parameters

		return(current_transition_parameters,current_emission_parameters)
	
	def MatrixBuild(self,current_transition_parameters,current_emission_parameters):
		'''Method accepts transition and emission matrices, converts them to numpy matrices, and normalizes them'''
		transitions_alphabet = {char:i for i,char in enumerate(self.transition_alphabet)} #{state:index}
		emissions_alphabet = {char:i for i,char in enumerate(self.emitted_alphabet)} #{char: index}
		### ------------------ Creating matrices with zeros ------------------###
		transition_matrix = np.zeros(shape=(len(self.transition_alphabet), len(self.transition_alphabet)))
		emission_matrix = np.zeros(shape=(len(self.transition_alphabet), len(self.emitted_alphabet)))
		
		for key in current_transition_parameters:
			for subkey in current_transition_parameters[key]:
				if sum(current_transition_parameters[key].values())==0:
					transition_matrix[transitions_alphabet[key],transitions_alphabet[subkey]]=1/len(current_transition_parameters[key])#equal transition probability if all coutns are 0
				else:
					transition_matrix[transitions_alphabet[key],transitions_alphabet[subkey]]=current_transition_parameters[key][subkey]/sum(current_transition_parameters[key].values())#normalization by sum of probabilities for each transition
		for key in current_emission_parameters:
			for subkey in current_emission_parameters[key]:
				if sum(current_emission_parameters[key].values())==0:
					emission_matrix[transitions_alphabet[key],emissions_alphabet[subkey]]=1/len(current_emission_parameters[key])#equal emission probability if all coutns are 0
				else:
					emission_matrix[transitions_alphabet[key],emissions_alphabet[subkey]]=current_emission_parameters[key][subkey]/sum(current_emission_parameters[key].values())#normalization by sum of probabilities for each emission
		return(transition_matrix,emission_matrix)	
	
	def wrap(self):
		'''Mehtod wraps iteration, and MatrixBuild to return new set of parameters, after the given number of iterations'''
		current_transition_parameters,current_emission_parameters = (self.iteration(self.transition_matrix,self.emission_matrix))
		return(self.MatrixBuild(current_transition_parameters,current_emission_parameters))

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
	res = ViterbiLearning(iterations,emitted_sequence,emitted_alphabet,transition_alphabet, transition_matrix, emission_matrix)
	transition_matrix, emission_matrix=res.wrap()

	### -------------Printing matrix ------------ ####
	print('\t'+'\t'.join(char for char in res.transition_alphabet))
	for i,e in enumerate(res.transition_alphabet):
		for index in range(len(transition_matrix)):
			if index==i:
				print(e+'\t'+'\t'.join(format(num,".3f") for num in transition_matrix[index]))			
	print('--------')
	print('\t'+'\t'.join(char for char in res.emitted_alphabet))
	for i,e in enumerate(res.transition_alphabet):
		for index in range(len(transition_matrix)):
			if index==i:
				print(e+'\t'+'\t'.join(format(num,".3f") for num in emission_matrix[index]))
	

if __name__ == "__main__":
	main()