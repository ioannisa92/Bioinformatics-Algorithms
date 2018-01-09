#!/usr/bin/env python3
# Name: Ioannis Anastopoulos
# Date: 12/10/2017

'''Program calculates transition and emission probabilities given a threshold theta, pseudocount, and a multiple alignment. Class
HMM() counts transition and emission ocurrences from multiple alignment, and build all non-forbidden transition states. Normalized
numpy arrays of transition and emission probabilities are printed. Class AlighnToHMMProfile() uses transition and emission matrices
from lass HMM() converts them to dict of dicts, and calculates the longest path in the Viterbi graph. AlighnToHMMProfile() used method backTracking()
to backtrack from sink to source in order to find the best (most likely path)'''
import sys
import numpy as np
from collections import defaultdict
import math
class HMM():
	def __init__(self, string_x,theta,pseudo,string_alphabet,alignment_length,grey_alignment_columns,clear_alignment_columns,m_alignment):
		''' string_x is the emitted string
			theta is the given threshold
			pseudo is the pseudocount given
		    alphabet is the chars that are emitted by the hidden states
		    alignment_length is the length of the alignment
		    grey_alignment_columns ar ethe insert columns that are above theta
		    clear_alignment_columns are the match/delete columns that are below theta'''
		self.string_x=string_x
		self.theta = theta
		self.pseudo = pseudo
		self.string_alphabet = string_alphabet
		self.alignment_length = alignment_length
		self.grey_alignment_columns = grey_alignment_columns
		self.clear_alignment_columns=clear_alignment_columns
		self.m_alignemt=m_alignment

	def transEmissionCounts(self):
		transition_counts=defaultdict(dict)
		transition_counts['S']={}#initializing transition dict
		transition_counts['S']['M1']=0
		transition_counts['S']['D1']=0
		transition_counts['S']['I0']=0
		emission_counts=defaultdict(dict)

		### --- Creating list of tuples(State,Emission) --- ###
		for item in self.m_alignemt:
			item_states=[]
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
		if 'I'+str(str(self.alignment_length-len(self.grey_alignment_columns))) not in transition_counts:
			transition_counts['I'+str(str(self.alignment_length-len(self.grey_alignment_columns)))]['E']=0
		if 'D'+str(str(self.alignment_length-len(self.grey_alignment_columns))) not in transition_counts:
			transition_counts['D'+str(str(self.alignment_length-len(self.grey_alignment_columns)))]['E']=0
		if 'M'+str(str(self.alignment_length-len(self.grey_alignment_columns))) not in transition_counts:
			transition_counts['M'+str(str(self.alignment_length-len(self.grey_alignment_columns)))]['E']=0
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
			for char in self.string_alphabet:
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
		'''method converts dict of dicts of transition and emission counts into numpy
		matrices and normalizes them'''
		transition_alphabet = [] #all non-forbiden states of the alignment
		for i in range(self.alignment_length-len(self.grey_alignment_columns)+1):
			for char in 'MD':
				if i!=0:
					transition_alphabet.append(char+str(i))
			transition_alphabet.append('I'+str(i))
		transition_alphabet.insert(0, 'S')
		transition_alphabet.append('E')

		transitions_alphabet = {char:i for i,char in enumerate(transition_alphabet)}
		emissions_alphabet = {char:i for i,char in enumerate(self.string_alphabet)}
		### ------------------ Creating matrices with zeros ------------------###
		transition_matrix = np.zeros(shape=(len(transition_alphabet), len(transition_alphabet)))
		emission_matrix = np.zeros(shape=(len(transition_alphabet), len(self.string_alphabet)))
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
		
		return(self.clear_alignment_columns,transition_alphabet, transition_matrix, emission_matrix)
		
	def wrap(self):
		transition_counts, emission_counts =self.transEmissionCounts()
		return(self.MatrixBuild(transition_counts,emission_counts))

class AlighnToHMMProfile():
	'''Class calculates most likley path over all possibile paths in a viterbi graph.Method Converter() converts numpy matrices
	of transition and emission probabilities from HMM() class. Method LongestPath() calculates max score and each 
	node(state) and the node the max score came from. Method backTracking() backtracks from sink to source returning most likley path'''
	def __init__(self,string_x,transition_alphabet, string_alphabet,transition_matrix, emission_matrix,clear_alignment_columns):
		self.string_x=string_x
		self.string_alphabet=string_alphabet
		self.transition_alphabet=transition_alphabet
		self.transition_matrix=transition_matrix
		self.emission_matrix=emission_matrix
		self.clear_alignment_columns=clear_alignment_columns

	def Converter(self):
		'''Method converts numpy matrices of transitions and emissions from class HMM() to dict of dicts'''
		transition_matrix_dict=defaultdict(dict)
		predecessor_states_dict=defaultdict(dict) #inverting transition_matrix_dict so that the dict is: {node:predecessors:transition prob from predecessor to node}
		emission_matrix_dict=defaultdict(dict)

		for i in range(len(self.transition_alphabet)-3): #current state
			for j in range(2,len(self.transition_alphabet)-2): #non forbidden transitions to states
				if i==0:
					transition_matrix_dict['S']['I'+str(i)]=self.transition_matrix[i,i+1]
					transition_matrix_dict['S']['M'+str(i+1)]=self.transition_matrix[i,i+2]
					transition_matrix_dict['S']['D'+str(i+1)]=self.transition_matrix[i,i+3]
				elif i==1:
					transition_matrix_dict['I'+str(i-1)]['I'+str(i-1)]=self.transition_matrix[i,i]
					transition_matrix_dict['I'+str(i-1)]['M'+str(i)]=self.transition_matrix[i,i+1]
					transition_matrix_dict['I'+str(i-1)]['D'+str(i)]=self.transition_matrix[i,i+2]
				else:
					if self.transition_matrix[i,j+2]!=0:#excluding forbidden transitions
						transition_matrix_dict[self.transition_alphabet[i]][self.transition_alphabet[j+2]]=self.transition_matrix[i,j+2]
		
		for k,v in transition_matrix_dict.items(): #current node
			for i,j in v.items(): #nodes current node is connected to and the transition proability
				predecessor_states_dict[i][k]=j #inverting transition matrix to get {node:predecessor:transition probability from predecessor to node}
		
		for i in range(1,len(self.transition_alphabet)): #current state
			for j in range(len(self.string_alphabet)): # emissions
				if 'D' in self.transition_alphabet[i]:
					emission_matrix_dict[self.transition_alphabet[i]][self.string_alphabet[j]]=np.float128(1) #since Delete states do not emit, probability is set to 1, to help scoring longest path later in LongestPath() method
				else:
					emission_matrix_dict[self.transition_alphabet[i]][self.string_alphabet[j]]=self.emission_matrix[i,j]
		return(transition_matrix_dict,predecessor_states_dict,emission_matrix_dict)

	def topologicalOrder(self):
		'''Method topologically orders viterbi graph'''
		topological_order_dict=defaultdict(list) #{key(index): list of states in each index}
		delete_states=[char for char in self.transition_alphabet if 'D' in char] #getting all delete states first
		topological_order=[states for states in self.transition_alphabet if states!='S' and states!='E']
		
		delete_states.insert(0,'S')#appending source in the beginning of the delete states
		for index in range(len(self.string_x)+1): #indexes are length fo string +1 to make up for the 0th index which does not emit
			if index == 0:
				topological_order_dict[index]=delete_states #source and delete states are in index 0
			else:
				topological_order_dict[index]=topological_order	
		return(topological_order_dict)

	def LongestPath(self, transition_matrix, predecessor_states_dict,emission_matrix, topological_order_dict):
		'''Method returns max score at each node, and the node used to calculate that score'''
		score_dict = defaultdict(dict) #saving score computations
		maxScore_predecessor_node_dict = dict() #dict of max_score: node used to calculate the max score
		for position in topological_order_dict:
			for state_index in range(len(topological_order_dict[position])):
				score_dict[position][topological_order_dict[position][state_index]] = -math.inf
				if position==0 and topological_order_dict[position][state_index]=='S':
					score_dict[position][topological_order_dict[position][state_index]]=np.longdouble(1) #initiating scoring at S=1
				elif position==0: #scoring all delete states at index=0, based on score of imediate previous delete state
					score_dict[position][topological_order_dict[position][state_index]]=np.longdouble(score_dict[position][topological_order_dict[position][state_index-1]])*(transition_matrix[topological_order_dict[position][state_index-1]][topological_order_dict[position][state_index]])
	
				else: #if delete state is NOT encountered, then the score of predecessor in PREVIOUS index in score_dict are used
					if 'D' not in topological_order_dict[position][state_index]:
						scores_list = []
						predecessor_and_prob = predecessor_states_dict[topological_order_dict[position][state_index]] #state current state is connected to, and the transition prob of the edge

						max_score = None
						for predecessor, prob in predecessor_and_prob.items():
							if predecessor in score_dict[position-1].keys():
			
								score = np.longdouble(score_dict[position-1][predecessor])*np.longdouble(prob)*np.longdouble(emission_matrix[topological_order_dict[position][state_index]][self.string_x[position-1]])
								scores_list.append((score))
								'''storing max score, and predecessor node used to calculate max score'''
								if max_score is None or score>=max_score:
									max_score = score
									max_score_node = predecessor
								maxScore_predecessor_node_dict[max_score] = max_score_node #finding the max score from all states and edges connected to current state
							score_dict[position][topological_order_dict[position][state_index]] = max(scores_list) if scores_list else score_dict[position][topological_order_dict[position][state_index]]
					if 'D' in topological_order_dict[position][state_index]: #if Delete state is encounteres, then scores of predecessores in the SAME index in score_dict are used
						scores_list = []
						predecessor_and_prob = predecessor_states_dict[topological_order_dict[position][state_index]] #state current state is connected to, and the transition prob of the edge

						max_score = None
						for predecessor, prob in predecessor_and_prob.items():
							if predecessor in score_dict[position].keys():
			
								score = np.longdouble(score_dict[position][predecessor])*np.longdouble(prob)*np.longdouble(emission_matrix[topological_order_dict[position][state_index]][self.string_x[position-1]])
								scores_list.append((score))
								'''storing max score, and predecessor node used to calculate max score'''
								if max_score is None or score>=max_score:
									max_score = score
									max_score_node = predecessor
								maxScore_predecessor_node_dict[max_score] = max_score_node #finding the max score from all states and edges connected to current state
							score_dict[position][topological_order_dict[position][state_index]] = max(scores_list) if scores_list else score_dict[position][topological_order_dict[position][state_index]] #if score_list is empty, then -inf is used
		return(score_dict, maxScore_predecessor_node_dict)

	def backTracking(self,longest_path_dict, maxScore_predecessor_node_dict):
		'''Method backtracks from sink to source, returning the longest path'''
		path = []			
		max_score =None
		for states_scores in longest_path_dict[len(self.string_x)].items(): #finding max of sink at last position of string_x to initialize backtracking
			if str(len(self.clear_alignment_columns))==states_scores[0][1:]: #selecting MDI states connected to sink (E)
				score = states_scores[1]
				if max_score is None or score>=max_score:
					max_score = score #max score
					max_score_node = states_scores[0] #state that max score at sink came from
		path.append(max_score_node) #initializing path from sink
		
		counter=len(self.string_x) #counting from last index to beginning in score_dict from LongestPath() method called longest_path_dict here

		while path[-1] !='S':
			if 'D' in path[-1]: #if delete state is appended to the path, we stay at the same index, meaning the same column in the Viterbi graph
				if longest_path_dict[counter][path[-1]] in maxScore_predecessor_node_dict.keys():	
					path.append(maxScore_predecessor_node_dict[longest_path_dict[counter][path[-1]]]) #appending node that gave rise to max score.
			else: #if ANY other state is appende to the path, we count down 1
				if longest_path_dict[counter][path[-1]] in maxScore_predecessor_node_dict.keys():	
					path.append(maxScore_predecessor_node_dict[longest_path_dict[counter][path[-1]]]) #appending node that gave rise to max score.
					counter-=1
		path.remove('S') #while loop append S as well, so it is removed
		return(' '.join(str(char) for char in path[::-1]))
	
	def wrap(self):
		topological_order_dict=self.topologicalOrder()
		transition_matrix_dict,predecessor_states_dict, emission_matrix_dict=self.Converter()
		longest_path_dict,maxScore_predecessor_node_dict=self.LongestPath(transition_matrix_dict,predecessor_states_dict, emission_matrix_dict,topological_order_dict)
		return(self.backTracking(longest_path_dict, maxScore_predecessor_node_dict))

def file_parse():
	'''Function parses input file and returns string, threshold,emission alphabet,length of alignment,
	list of column indexes of insert columns, list of column indexes of match/delete columns, and multiple alignment'''
	m_alignment = []
	grey_alignment_columns=[]
	clear_alignment_columns=[]
	with sys.stdin as fn:
		lines = fn.readlines()
		string_x=lines[0].rstrip()
		theta_pseudo= (lines[2]).split()
		theta=float(theta_pseudo[0].rstrip())
		pseudo=float(theta_pseudo[1].rstrip())
		string_alphabet = lines[4].split()
		for line in lines[6:]:
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

	return(string_x,theta,pseudo,string_alphabet,alignment_length,grey_alignment_columns,clear_alignment_columns,m_alignment)

def main():
	string_x,theta,pseudo,string_alphabet,alignment_length,grey_alignment_columns,clear_alignment_columns,m_alignment=file_parse() #getting inputs from file
	
	profile=HMM(string_x,theta,pseudo,string_alphabet,alignment_length,grey_alignment_columns,clear_alignment_columns,m_alignment) #passing inputs to HMM() class to get transition and emission matrices
	clear_alignment_columns,transition_alphabet, transition_matrix, emission_matrix= profile.wrap() #using wrap method of HMM() class, which returns transition and emission matrices as numpy matrices

	alignment =AlighnToHMMProfile(string_x,transition_alphabet,string_alphabet, transition_matrix, emission_matrix, clear_alignment_columns) #passing outputs from HMM() class to AlighnToHMMProfile() class
	transition_matrix_dict,predecessor_states_dict, emission_matrix_dict = alignment.Converter() #using Converter() method to convert numpy matrices from HMM() class to dict of dicts, because thats the data structure methods LongestPath() and backTracking() work with

	print(alignment.wrap()) #using wrap method from AlighnToHMMProfile() class to print out most likely hidden path

if __name__ == "__main__":
	main()