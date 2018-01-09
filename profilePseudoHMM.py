#!/usr/bin/env python3
# Name: Ioannis Anastopoulos
# Date: 11/20/2017

'''Program calculates transition and emission probabilities given a threshold theta, pseudocount, and a multiple alignment. Class
HMM() counts transition and emission ocurrences from multiple alignment, and build all non-forbidden transition states. Normalized
numpy arrays of transition and emission probabilities are printed '''

import sys
import numpy as np
from collections import defaultdict
class HMM():
	def __init__(self, theta,pseudo,alphabet,alignment_length,grey_alignment_columns,clear_alignment_columns,m_alignment):
		''' theta is the given threshold
			pseudo is the pseudocount given
		    alphabet is the chars that are emitted by the hidden states
		    alignment_length is the length of the alignment
		    grey_alignment_columns ar ethe insert columns that are above theta
		    clear_alignment_columns are the match/delete columns that are below theta'''
		self.theta = theta
		self.pseudo = pseudo
		self.alphabet = alphabet
		self.alignment_length = alignment_length
		self.grey_alignment_columns = grey_alignment_columns
		self.clear_alignment_columns=clear_alignment_columns
		self.m_alignment=m_alignment

	def transEmissionCounts(self):
		transition_counts=defaultdict(dict)
		transition_counts['S']={}#initializing transition dict
		transition_counts['S']['M1']=0
		transition_counts['S']['D1']=0
		transition_counts['S']['I0']=0
		emission_counts=defaultdict(dict)
		### --- Creating list of tuples(State,Emission) --- ###
		for item in self.m_alignment:#iterating through each sequence in alignment
			item_states=[]#list of tuples(state, emission) of each sequence in alignment
			counter=0
			for i, char in enumerate(item):
				if i in self.clear_alignment_columns and char!="-":
					counter+=1 #increment each time you are in clear column
					item_states.append(('M'+str(counter), char))
				if i in self.clear_alignment_columns and char=="-":
					counter+=1#increment each time you are in clear column
					item_states.append(('D'+str(counter), char))
				if i in self.grey_alignment_columns and char!="-":
					item_states.append(('I'+str(counter), char)) #counter not incremented if i not in clear
		### --- Creating list of tuples(State,Emission) --- ###
		
		### --- Counting Transitions and Emissions --- ###
			transition_counts['S'][item_states[0][0]]+=1

			if 'E' not in transition_counts[item_states[len(item_states)-1][0]]:
				transition_counts[item_states[len(item_states)-1][0]]['E']=0
			transition_counts[item_states[len(item_states)-1][0]]['E']+=1
			
			for i in range(len(item_states)-1):
				if item_states[i][0] in transition_counts and item_states[i+1][0] in transition_counts[item_states[i][0]]:
					transition_counts[item_states[i][0]][item_states[i+1][0]]+=1
				else:
					transition_counts[item_states[i][0]][item_states[i+1][0]]=1		

			for i in range(len(item_states)):
				if item_states[i][0] in emission_counts and item_states[i][1] in emission_counts[item_states[i][0]]:
					emission_counts[item_states[i][0]][item_states[i][1]]+=1
				else:
					emission_counts[item_states[i][0]][item_states[i][1]]=1
		### --- Counting Transitions and Emissions --- ###

		###---Filling in transition and emission matrices with allowed states and emissions---###
		if 'I'+str(self.alignment_length-len(self.grey_alignment_columns)) not in transition_counts:
			transition_counts['I'+str(self.alignment_length-len(self.grey_alignment_columns))]['E']=0
		if 'D'+str(self.alignment_length-len(self.grey_alignment_columns)) not in transition_counts:
			transition_counts['D'+str(self.alignment_length-len(self.grey_alignment_columns))]['E']=0
		if 'M'+str(self.alignment_length-len(self.grey_alignment_columns)) not in transition_counts:
			transition_counts['M'+str(self.alignment_length-len(self.grey_alignment_columns))]['E']=0
		for i in range(self.alignment_length-len(self.grey_alignment_columns)):
			if i==0:
				if 'I'+str(i) not in transition_counts:
					transition_counts['I'+str(i)]={}
					transition_counts['I'+str(i)]['I'+str(i)]=0
					transition_counts['I'+str(i)]['M'+str(i+1)]=0
					transition_counts['I'+str(i)]['D'+str(i+1)]=0
			else:
				if 'I'+str(i) not in transition_counts:
					transition_counts['I'+str(i)]={}
					transition_counts['I'+str(i)]['I'+str(i)]=0
					transition_counts['I'+str(i)]['M'+str(i+1)]=0
					transition_counts['I'+str(i)]['D'+str(i+1)]=0
				if 'M'+str(i) not in transition_counts:
					transition_counts['M'+str(i)]={}
					transition_counts['M'+str(i)]['I'+str(i)]=0
					transition_counts['M'+str(i)]['M'+str(i+1)]=0
					transition_counts['M'+str(i)]['D'+str(i+1)]=0	
				if 'D'+str(i) not in transition_counts:
					transition_counts['D'+str(i)]={}
					transition_counts['D'+str(i)]['I'+str(i)]=0
					transition_counts['D'+str(i)]['M'+str(i+1)]=0
					transition_counts['D'+str(i)]['D'+str(i+1)]=0
		for key in transition_counts:
			if key=='I0':
				if 'I'+str(int(key[1])) not in transition_counts[key]:
					transition_counts[key]['I'+str(int(key[1]))]=0
				if 'D'+str(int(key[1])+1) not in transition_counts[key]:
					transition_counts[key]['D'+str(int(key[1])+1)]=0
				if 'M'+str(int(key[1])+1) not in transition_counts[key]:
					transition_counts[key]['M'+str(int(key[1])+1)]=0
			if key!='S' and key!='I0':
				if key[1:]!=str(self.alignment_length-len(self.grey_alignment_columns)):
					if 'I'+str(int(key[1:])) not in transition_counts[key]:
						transition_counts[key]['I'+str(int(key[1:]))]=0
					if 'D'+str(int(key[1:])+1) not in transition_counts[key]:
						transition_counts[key]['D'+str(int(key[1:])+1)]=0
					if 'M'+str(int(key[1:])+1) not in transition_counts[key]:
						transition_counts[key]['M'+str(int(key[1:])+1)]=0
				elif key[1:]==str(self.alignment_length-len(self.grey_alignment_columns)):
					if 'I'+key[1:] not in transition_counts[key]:
						transition_counts[key]['I'+key[1:]]=0
					if 'E' not in transition_counts[key]:
						transition_counts[key]['E']=0

		for i in range(self.alignment_length-len(self.grey_alignment_columns)+1):
			if i==0:
				if 'I'+str(i) not in emission_counts:
					emission_counts['I'+str(i)]={}
			else:
				if 'M'+str(i) not in emission_counts:
					emission_counts['M'+str(i)]={}
			for char in self.alphabet:
				if i==0:
					if char not in emission_counts['I'+str(i)]:
						emission_counts['I'+str(i)][char]=0
				else:
					if char not in emission_counts['I'+str(i)]:
						emission_counts['I'+str(i)][char]=0
					if char not in emission_counts['M'+str(i)]:
						emission_counts['M'+str(i)][char]=0	
		for key, value in emission_counts.items():
			if '-' in value:
				del emission_counts[key]['-']

		###---Filling in transition and emission matrices with allowed states and emissions---###

		return(transition_counts, emission_counts)
	
	def MatrixBuild(self, transition_counts, emission_counts):
		'''Method converts transition_counts and emission_counts from dict of dicts to numpy matrices
		   and normalizes them'''
		transition_alphabet = [] #all non-forbiden states of the alignment
		for i in range(self.alignment_length-len(self.grey_alignment_columns)+1):
			for char in 'MD':
				if i!=0:
					transition_alphabet.append(char+str(i))
			transition_alphabet.append('I'+str(i))
		transition_alphabet.insert(0, 'S')
		transition_alphabet.append('E')

		transitions_alphabet = {char:i for i,char in enumerate(transition_alphabet)}
		emissions_alphabet = {char:i for i,char in enumerate(self.alphabet)}
		### ------------------ Creating matrices with zeros ------------------###
		transition_matrix = np.zeros(shape=(len(transition_alphabet), len(transition_alphabet)))
		emission_matrix = np.zeros(shape=(len(transition_alphabet), len(self.alphabet)))
		### --- Filling in matrices and normalizing them at the same time --- ###
		for key in transition_counts:
			for subkey in transition_counts[key]:
				if sum(transition_counts[key].values())==0 and len(transition_counts[key])==3:
					transition_matrix[transitions_alphabet[key],transitions_alphabet[subkey]]=1/3
				elif sum(transition_counts[key].values())==0 and len(transition_counts[key])==2:
					transition_matrix[transitions_alphabet[key],transitions_alphabet[subkey]]=1/2
				else:
					transition_matrix[transitions_alphabet[key],transitions_alphabet[subkey]]=np.float64(transition_counts[key][subkey]/sum(transition_counts[key].values())+self.pseudo)/(1+len(transition_counts[key])*self.pseudo)

		for key in emission_counts:
			for subkey in emission_counts[key]:
				if sum(emission_counts[key].values())==0:
					emission_matrix[transitions_alphabet[key],emissions_alphabet[subkey]]=1/len(emission_counts[key])
				else:
					emission_matrix[transitions_alphabet[key],emissions_alphabet[subkey]]=np.float64(emission_counts[key][subkey]/sum(emission_counts[key].values())+self.pseudo)/(1+len(emission_counts[key])*self.pseudo)

		### --- Filling in matrices and normalizing them at the same time --- ###
		
		return(transition_alphabet, transition_matrix, emission_matrix)
		
	def wrap(self):
		transition_counts, emission_counts =self.transEmissionCounts()
		return(self.MatrixBuild(transition_counts,emission_counts))

def file_parse():
	'''Function parses input file and returns threshold,emission alphabet,length of alignment,
	   list of column indexes of insert columns, list of column indexes of match/delete columns, and multiple alignment'''
	m_alignment = []
	grey_alignment_columns=[]
	clear_alignment_columns=[]
	with sys.stdin as fn:
		lines = fn.readlines()
		theta_pseudo= (lines[0]).split()
		theta=float(theta_pseudo[0].rstrip())
		pseudo=float(theta_pseudo[1].rstrip())
		alphabet = lines[2].split()
		for line in lines[4:]:
			m_alignment.append(line.rstrip())
			alignment_length = len(line)#length of alignment 
		space_occurence ={i:0 for i in range(alignment_length)} #creating dict of each position and count of spaces in each position

	for item in m_alignment:
		for i in range(len(item)):
			if item[i]=='-':
				space_occurence[i]+=1#counting spaces in each column

	for key in space_occurence:
		space_occurence[key] = space_occurence[key]/len(m_alignment)#occurence of space in each column
		if space_occurence[key] >= theta:
			grey_alignment_columns.append(key) #columns that exceed theta

	for item in list(range(alignment_length)):
		if item not in grey_alignment_columns:
			clear_alignment_columns.append(item)# columns that do not exceed theta

	return(theta,pseudo,alphabet,alignment_length,grey_alignment_columns,clear_alignment_columns,m_alignment)

def main():
	theta,pseudo,alphabet,alignment_length,grey_alignment_columns,clear_alignment_columns,m_alignment=file_parse()
	res=HMM(theta,pseudo,alphabet,alignment_length,grey_alignment_columns,clear_alignment_columns,m_alignment)
	transition_alphabet, transition_matrix, emission_matrix= res.wrap()
	### -------------Printing matrix ------------ ####
	print('\t'+'\t'.join(char for char in transition_alphabet))
	for i,e in enumerate(transition_alphabet):
		for index in range(len(transition_matrix)):
			if index==i:
				print(e+'\t'+'\t'.join(format(num,".3f") for num in transition_matrix[index]))			
	print('--------')
	print('\t'+'\t'.join(char for char in res.alphabet))
	for i,e in enumerate(transition_alphabet):
		for index in range(len(transition_matrix)):
			if index==i:
				print(e+'\t'+'\t'.join(format(num,".3f") for num in emission_matrix[index]))
	# print('\n')
if __name__ == "__main__":
	main()