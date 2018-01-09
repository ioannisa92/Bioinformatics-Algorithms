#!/usr/bin/env python3
# Name: Ioannis Anastopoulos
# Date: 11/29/2017

'''Program calculates matrix of transition probabilities 
and emission probabilities that maximize Pr(x,π) over all possible matrices of transition and emission probabilities.
Class HMMParameterEstimation() counts transitions and emissions with transitionCounts() and emissionCounts() methods, 
and converts these dict of dicts to normalized numpy arrays with MatrixBuild().'''

import sys
import numpy as np
from collections import defaultdict

class HMMParameterEstimation():
	def __init__(self, emitted_sequence,emission_alphabet,path,path_alphabet):
		self.emitted_sequence=emitted_sequence
		self.emission_alphabet=emission_alphabet
		self.path=path
		self.path_alphabet=path_alphabet
	
	def transitionCounts(self):
		'''Method calculates trasition counts by iterating through each char in the given path, and looking at the
		char next to that.'''
		transition_counts=defaultdict(dict)
		for i in range(len(self.path)-1):
			for char in self.path_alphabet:
				if self.path[i]==char:
					if self.path[i] in transition_counts and self.path[i+1] in transition_counts[self.path[i]]:
						transition_counts[self.path[i]][self.path[i+1]]+=1
					else:
						transition_counts[self.path[i]][self.path[i+1]]=1
				if char not in transition_counts:
					transition_counts[char]={chars:0 for chars in self.path_alphabet}
		for key in transition_counts:
			for char in self.path_alphabet:
				if char not in transition_counts[key]:
					transition_counts[key][char]=0
		return(transition_counts)
	def emissionCounts(self):
		'''Method calculates emission counts by iterating through indexes of both path and emission (same length)
		and counting the emission of each char in the path'''
		emission_counts=defaultdict(dict)

		for i in range(len(self.emitted_sequence)):
			if self.path[i] not in emission_counts:
				emission_counts[self.path[i]]={}
			if self.emitted_sequence[i] not in emission_counts[self.path[i]]:
				emission_counts[self.path[i]][self.emitted_sequence[i]]=1
			else:
				emission_counts[self.path[i]][self.emitted_sequence[i]]+=1
		for key in emission_counts:
			for symbol in self.emission_alphabet:
				if symbol not in emission_counts[key]:
					emission_counts[key][symbol]=0
		for char in self.path_alphabet:
			for symbol in self.emission_alphabet:
				if char not in emission_counts:
					emission_counts[char]={chars:0 for chars in self.emission_alphabet}
		return(emission_counts)
	
	def MatrixBuild(self,transition_counts,emission_counts):
		'''Method converts transition_counts and emission_counts from dict of dicts to numpy matrices
		   and normalizes them'''
		transitions_alphabet = {char:i for i,char in enumerate(self.path_alphabet)}
		emissions_alphabet = {char:i for i,char in enumerate(self.emission_alphabet)}
		### ------------------ Creating matrices with zeros ------------------###
		transition_matrix = np.zeros(shape=(len(self.path_alphabet), len(self.path_alphabet)))
		emission_matrix = np.zeros(shape=(len(self.path_alphabet), len(self.emission_alphabet)))
		### --- Filling in matrices and normalizing them at the same time --- ###
		if len(transition_counts)==0:
			for i in range(len(self.path_alphabet)):
				for j in range(len(self.path_alphabet)):
					transition_matrix[i,j]=1/len(set(self.path_alphabet))
		for key in transition_counts:
			for subkey in transition_counts[key]:
				if sum(transition_counts[key].values())==0:
					transition_matrix[transitions_alphabet[key],transitions_alphabet[subkey]]=1/len(transition_counts[key])
				elif key==subkey and len(transition_counts)==1:
					for i in range(len(self.path_alphabet)):
						for j in range(len(self.path_alphabet)):
							transition_matrix[i,j]=np.float64(transition_counts[key][subkey])/sum(transition_counts[key].values())
				else:
					transition_matrix[transitions_alphabet[key],transitions_alphabet[subkey]]=np.float64(transition_counts[key][subkey])/sum(transition_counts[key].values())

		if len(emission_counts)==1 and len(self.path_alphabet)!=1:
			for i in range(len(self.path_alphabet)):
				for j in range(len(self.emission_alphabet)):
					emission_matrix[i,j]=np.float64(emission_counts[self.path_alphabet[i]][self.emission_alphabet[j]])/sum(emission_counts[self.path_alphabet[i]].values())
		for key in emission_counts:
			for subkey in emission_counts[key]:
				if sum(emission_counts[key].values())==0:
					emission_matrix[transitions_alphabet[key],emissions_alphabet[subkey]]=1/len(set(self.emission_alphabet))
				else:
					emission_matrix[transitions_alphabet[key],emissions_alphabet[subkey]]=np.float64(emission_counts[key][subkey])/sum(emission_counts[key].values())
		### --- Filling in matrices and normalizing them at the same time --- ###
		return(transition_matrix,emission_matrix)

	def wrap(self):
		transition_counts=self.transitionCounts()
		emission_counts=self.emissionCounts()
		return(self.MatrixBuild(transition_counts,emission_counts))

def file_parser():
	with sys.stdin as fn:
		lines=fn.readlines()
		emitted_sequence=lines[0].rstrip()
		emission_alphabet=lines[2].split()
		path=lines[4].rstrip()
		path_alphabet=lines[6].split()
	return(emitted_sequence,emission_alphabet,path,path_alphabet)

def main():
	emitted_sequence,emission_alphabet,path,path_alphabet=file_parser()
	res=HMMParameterEstimation(emitted_sequence,emission_alphabet,path,path_alphabet)
	transition_matrix, emission_matrix=res.wrap()
	### -------------Printing matrix ------------ ####
	print('\t'+'\t'.join(char for char in res.path_alphabet))
	for i,e in enumerate(res.path_alphabet):
		for index in range(len(transition_matrix)):
			if index==i:
				print(e+'\t'+'\t'.join(format(num,".3f") for num in transition_matrix[index]))			
	print('--------')
	print('\t'+'\t'.join(char for char in res.emission_alphabet))
	for i,e in enumerate(res.path_alphabet):
		for index in range(len(transition_matrix)):
			if index==i:
				print(e+'\t'+'\t'.join(format(num,".3f") for num in emission_matrix[index]))
if __name__ == "__main__":
	main()