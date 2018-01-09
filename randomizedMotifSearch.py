#!/usr/bin/env python3
# Name: Ioannis Anastopoulos
# Date: 10/12/2017
# Group members: Alexis, Xiang, Balaji, 

"""Program accepts a fasta file that contains multiple fasta records as input (STDIN)
and outputs (STDOUT) the consensus sequence and associated profile score, by conducting randomized motif search. 
The score is the sum of entropies across each position in the final profile.
Pseudocounts are used for this assignment, which is the -p option.


Randomized mootif can be easily trapped in local minima, so we will need to iterate some number of times to find the best results. 
This is an option (-i) that establishes the iteration number.

-k option sepcifies length of motif to be searched
Usage looks like this:
python3 randomizedMotifSearch.py -i=100000 -p=1 -k=13 <input.fa >soutput.fa

Program has been used with fasta files that represent the upstream 50 bases from species of the Pyrobaculum group.

Program output:
Consensis sequence:
Entropy:
"""


import sys
import random
import math


class CommandLine() :
	"""CommandLine class to import args from commandline."""


	def __init__(self, inOpts=None) :
		import argparse

		self.parser = argparse.ArgumentParser(description="This program accepts mutliple Fasta files and outputs the consensus promoter motif sequence upstream of CRISPR arrays and associated entropy score based on the randomized motif search algorithm. Pseudocounts cannot be less than 1, because in case the probability of a base in a position is 0, entropy will be infinity")
		'''First argument sets number of iterations through the program'''
		self.parser.add_argument('-i', '--iterations', type=int, default=1000)
		'''Second argument sets the pseudocount number.'''
		self.parser.add_argument('-p', '--pseudocounts', type=int, default=1) #pseudocounts of less than 1 are not accepted. See program description
		'''Third arguments sets the length of the kmer.'''
		self.parser.add_argument('-k', '--kmer', type=int, default=0)

		if inOpts is None :
			self.args = self.parser.parse_args()
		else :
			self.args = self.parser.parse_args(inOpts)

class FastAreader :
	
	def __init__ (self, fname=''):
		'''contructor: saves attribute fname '''
		self.fname = fname
			
	def doOpen (self):
		if self.fname is '':
			return sys.stdin
		else:
			return open(self.fname)
 
	def readFasta (self):
		
		header = ''
		sequence = ''
		
		with self.doOpen() as fileH:
			
			header = ''
			sequence = ''
 
			# skip to first fasta header
			line = fileH.readline()
			while not line.startswith('>') :
				line = fileH.readline()
			header = line[1:].rstrip()

			for line in fileH:
				if line.startswith ('>'):
					yield header,sequence
					header = line[1:].rstrip()
					sequence = ''
				else :
					sequence += ''.join(line.rstrip().split()).upper()
						
		yield header,sequence


class consensusMotifSearch():
	'''Class calculates probability profile and associated entropy of motifs, and puts together a consesus sequence with teh lowerst entropy score'''


	def __init__(self, seq_list, k, p, i):
		self.k = k #option from command line for kmer length
		self.p = p #option from command line for pseudocounts
		self.i = i #option of iterations from command line
		self.seq_list = seq_list #storing fasta sequences in list

	def randomKmers(self):
		'''method stores random motifs from each sequence in the fasta file'''
		randomMotifs = list()

		for seq in self.seq_list:
			kmer_index = random.randint(0, len(seq)-self.k+1) #the length fo the is subtracted by the length of the kmer, so that we dont run out of index
			random_kmer =(seq[kmer_index:kmer_index+self.k])
			randomMotifs.append(random_kmer)
		return randomMotifs

	def profileCounts(self, kmer_list): #kmer_list paramter is the list of motifs(either random or new motifs from new profile)
		'''method counts the occurrence of each base in each position of the inputed kmer-length'''
		
		profile_counts_dict ={'A': self.k*[self.p], #pseudocounts stored beforehand for each base
							  'G': self.k*[self.p], #pseudocounts of less than 1 are not accepted. See program description
							  'T': self.k*[self.p], 
							  'C': self.k*[self.p]}
		for i in kmer_list:
			for j in range(len(i)):
				profile_counts_dict[i[j]][j] +=1 #updates counts of profile_counts_dict, adding to the pseudocounts
		return profile_counts_dict

	def probabilityProfile(self, kmer_lst):
		'''method to calculate probability porfile, of inputed kmer_lst(list of kmers)'''
		lst = kmer_lst

		profile_prob_dict = dict()
		#divide each count by the total to find the probability
		total = len(self.seq_list)+4*self.p

		for k,v in self.profileCounts(lst).items():
			profile_prob_dict[k] = [i/total for i in v] #updates profile_prob_dict dictionary with the probability of each base in position of the kmer
		return profile_prob_dict

	def Entropy(self, kmer_lst):
		'''method calculates entropy of inputed kmer_lst(kmer list)'''
		
		score = 0
		for k,v in self.probabilityProfile(kmer_lst).items():
			for i in range(self.k):
				score =score+(v[i]*math.log(v[i],2))
		return(-1*score) #summing all scores together

	def newKmerMotifs(self, probProfile):
		'''method creates a list of new kmer motifs based on a probability profile'''
		mostProbMotifs = list()

		for seq in self.seq_list: #for each sequence in the list of sequences from the fasta file
			
			high_prob=0
			high_motif = ''
			
			for i in range(len(seq)-self.k+1):
				slide_kmer = seq[i:i+self.k]

				new_kmer_prob = 1

				'''this is a little faster'''
				# c=0
				# for char in slide_kmer:
				# 	new_kmer_prob *= probProfile[char][c] #total probability of each slider slice (so each kmer in each sequence)
				# 	c+=1
				'''this is a little faster'''

				for v in range(0,self.k):
					new_kmer_prob *= probProfile[slide_kmer[v]][v] #total probability of each slider slice (so each kmer in each sequence)

				if new_kmer_prob > high_prob: #choosing highest probability kmer
					high_prob = new_kmer_prob
					high_motif = slide_kmer
			mostProbMotifs.append(high_motif) #appending higher probability kemrs in list

		return mostProbMotifs


	def randomizedMotifSearch(self):
		'''method returns the best profile with the lowest entropy'''

		motifs = self.randomKmers()

		bestMotifs=motifs
		###---------------profile of random motifs-----------------###
		profile_prob = self.probabilityProfile(bestMotifs)
		entropy = self.Entropy(bestMotifs)

		###---------------profile of random motifs-----------------###

		bestProfile=profile_prob
		bestEntropy=entropy

		# newMotifs = self.newKmerMotifs(bestProfile) #calling method to return new kmers based on previous profile

		while True:		
			newMotifs = self.newKmerMotifs(bestProfile) #calling method to return new kmerrs based on previous profile
			
			###---------------profile of new motifs-----------------###
			new_profile_prob = self.probabilityProfile(newMotifs)
			new_entropy = self.Entropy(newMotifs)
			###---------------profile of new motifs-----------------###

			if new_entropy<bestEntropy: #finding minimum entropy
				bestMotifs = newMotifs #updating bestMorifs with motifs associated with minimum entropy
				bestProfile = new_profile_prob #best profile based on lowest entropy
				bestEntropy = new_entropy
			else:
				return (bestProfile, bestEntropy)
		
	def iteration(self):
		'''method iterates through randomizedMotifSearch() method as many times as the user chooses
		returning profile with minimum entropy'''
		
		minimum_Entropy = None

		for i in range(self.i):
			Profile, Score = self.randomizedMotifSearch()	
			if minimum_Entropy is None or Score<minimum_Entropy:
				minimum_Entropy = Score
				best_Profile = Profile
		return (best_Profile, minimum_Entropy)

	def results(self):
		'''method puts together consensus sequence based on best profile and minimum entropy from iteration method'''
		best_Profile, minimum_entropy= self.iteration()

		consensus_sequence=''

		for j in range(self.k):
			consensus = list() #for each iteration a new list of the probability of each base in each position of the kmer is created
			for k,v in best_Profile.items():
				consensus.append((k,v[j])) #base,probability tuple in list
			t = (max(consensus, key=lambda x:x[1])) #picking highest probability base in the list of (base,probability)
			consensus_sequence += t[0] #building consensus sequence

		return (consensus_sequence, minimum_entropy)


def main():
	'''Calling CommandLine class'''
	command_line = CommandLine()

	'''Calling FastAreader class'''
	reading = FastAreader()
	'''Calling the readFasta method '''
	head_seq = reading.readFasta()

	'''storing all sequences in Fasta file in a list'''
	seq_list = list()
	for h,s in head_seq:
		clean_seq = s.replace('\\', "")
		seq_list.append(clean_seq)

	'''Setting conditions for when the program Errors out, exiting the program'''
	if command_line.args.pseudocounts == 0: 
		print("ERROR -p --pseudocounts cannot be 0. See program description")#pseudocounts cannot be less than 1. See program description for explanation
		sys.exit(-1)
	elif command_line.args.kmer < 1: 
		print("ERROR -k --kmer length cannot be less 1")
		sys.exit(-1)
	elif command_line.args.iterations < 1:
		print("ERROR -i --iterations cannot be less than 1")
		sys.exit(-1)
	else:
		res = consensusMotifSearch(seq_list, command_line.args.kmer, command_line.args.pseudocounts, command_line.args.iterations)
		consensus,entropy = res.results()
		print('Consensus Sequence: ', consensus)
		print('Entropy: ', entropy)


if __name__ == "__main__":
	main()


