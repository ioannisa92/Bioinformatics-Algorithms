#!/usr/bin/env python3
# Name: Ioannis Anastopoulos
# Date: 10/03/2017
# Members: Alexis Thornton, Lori Siao, Balaji Sundararaman, Noah Dove



"""This program  reads a fasta file from STDIN and ranks motifs based on how statistically underrepresented the specific motif is. 
Program considers sequence motifs that are up to 8nt in length (this is a program option). 
In order to use the equation developed in class, the minimum size should be at least 3 (this is a program option)
Program uses Z-scores, and since we are looking for underepresented motifs, we will be looking atm negative Z-scores (Z-score cutoff is also an option) 
Program sorts the output by z-score and print it out to STDOUT """

import itertools
import sys


class CommandLine() :
	"""CommandLine class to import args from commandline."""


	def __init__(self, inOpts=None) :
		import argparse

		self.parser = argparse.ArgumentParser("This program counts all the possible kmers in the genome and outputs the under-represented kmers. Default minimum, maximum kmer is 3 and 8, and default z-score of 0 is used")
		'''First argument is the minimum kmer length desired.'''
		self.parser.add_argument('--minMotif', '--minMotif', type=int, default=3)
		'''Second argument is the maximum kmer length desired.'''
		self.parser.add_argument('--maxMotif', '--maxMotif', type=int, default=8)
		'''Third argument is the z-score cutoff below which the calculations are going to done.'''
		self.parser.add_argument('--cutoff', '--cutoff', type=float, default=0)

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
			# initialize return containers
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


class Motif:
	"""Motif class counts kmers, compares counts to null model, calculates sd score, z score, and compiles results to be printed"""

	def __init__ (self, min_kmer, max_kmer, z_scoreCutoff):
		'''Instantiating dictionaries for data storage.'''
		self.kmer_observed = dict()
		self.expected_counts_dict = dict()
		self.sd_dict = dict()
		self.zScore_dict = dict()

		'''Instantiating tuple in which the output data are stored.'''
		self.resTuple = tuple()
		'''Instantiating list for storage of tupples.'''
		self.res_lst = list()

		
		'''Instantiating length of genome.'''
		self.length = 0

		'''Arguments from command line are instantiated here.'''
		self.min_kmer = min_kmer
		self.max_kmer = max_kmer
		self.z_Score_cutoff = z_scoreCutoff
	
	def lengthOfSequence(self, s):
		'''Method calculates length of genome'''

		'''Length of sequence will be calculated, by ingoring non AGCT chars.'''
		des_bases = ['A', 'G', 'C', 'T']
		'''Iteration through desired bases, checkingthe chars in sequence, and counting, storing result in self.length.'''
		for char in s:
			if char in des_bases:
				self.length += 1

	def kmerComposition(self, s):
		'''Method counts the observed kmers in the genome using a dict'''

		'''Iterating through all possible kmer values.'''
		for i in range((self.min_kmer-(self.min_kmer-1)),self.max_kmer+1):
			psbl_kmers = len(s)-i+1
			'''Iterating through all possible kmer sizes.'''
			for k in list(range(0,psbl_kmers)):
				'''Creating substring, meaning actual kmer, storing in dict, and counting.'''
				slide_k = s[k:k+i]
				if slide_k not in self.kmer_observed:
					self.kmer_observed[slide_k] = 1 
				elif slide_k in self.kmer_observed:
					self.kmer_observed[slide_k] += 1


	def addData(self):
		'''Method counts all the possible kmers,
		and checks if these kmers are in the observed kmer dict'''

		des_bases = ['A', 'G', 'C', 'T'] #list of desired bases


		'''Iterating through all possible kmer values.'''
		for i in range((self.min_kmer-(self.min_kmer-1)), self.max_kmer+1):
			'''All possible permutations of kmers using desired bases.'''
			possible_kmers = itertools.product(des_bases, repeat=i)
			for k in possible_kmers:
				seq = ''.join(k)
				'''If a kmer permutation is not in the observed count dict kmer_observed, the kmer
				is included, but count is set to 0, as a pseudocount.'''
				if seq not in self.kmer_observed:
					self.kmer_observed[seq] = 0
		
	def expectedCounts(self):
		'''Method calculates expected count of kmers'''

		'''Storing the keys from the observed dict, kmer_obserced, with default value of 0.'''
		'''kmers of size>3 are selected, because the expected counts cannot be calculated lower sizes.'''
		self.expected_counts_dict = {k:0 for k,v in self.kmer_observed.items() if len(k)>=3}

		'''Calculating exepcted counts of kmer using a null model.'''
		for k in self.expected_counts_dict:
			num_1 = k[:-1]
			num_2 = k[1:]
			denom_1 = k[1:-1]
			'''Populating the dictionary with the expected counts.'''
			try:
				self.expected_counts_dict[k] = (self.kmer_observed[num_1]*self.kmer_observed[num_2])/(self.kmer_observed[denom_1])
			#ignore division with 0
			except ZeroDivisionError:
				pass

	def sdScore(self):
		'''Method calculates sd score'''

		import math
		'''Calculating sd sdore using expected counts dictionary, expected_counts_dict'''
		for k,v in self.expected_counts_dict.items():
			sd = math.sqrt((v)*(1-(v/self.length)))
			self.sd_dict[k] = sd

	def zScore(self):
		'''Method calculates z-score'''

		'''Calculating z score using expected counts dictionary, expected_counts_dict'''
		for k,v in self.expected_counts_dict.items():
			'''Populating the zscore dict, zScore_dict, with z scores.'''
			try:
				self.zScore_dict[k] = ((self.kmer_observed[k]) -v)/(self.sd_dict[k])
			except ZeroDivisionError: #Ignoring division with 0
				pass

	def results(self):
		'''Method is used to create tuples of results containing (kmer, observed count, expected count, z-score, p value).'''

		'''iterating through zscore dictionary'''
		for k,v in self.zScore_dict.items():
			if len(k) >= self.min_kmer: #condition useful for minKmer option of >3

				if v<self.z_Score_cutoff: #choosing values below the z-score cutoff

					observed_count = self.kmer_observed[k]

					expected_count = self.expected_counts_dict[k]

					z_score = v

					#tuple of results
					self.resTuple = (k,observed_count,expected_count,z_score)
					self.res_lst.append(self.resTuple) #tuples are stored in a list to be sorted in main



def main():
	"""Main function to print out results in a sorted fashion"""

	'''Calling CommandLine class'''
	command_line = CommandLine()
	
	'''Calling FastAreader class'''
	reading = FastAreader()
	'''Calling the readFasta method '''
	head_seq = reading.readFasta()

	'''Passing command line arguments, minimum kmer size, maximum kmer size, and z-score cutoff, in the Motif class.'''
	res = Motif(command_line.args.minMotif, command_line.args.maxMotif, command_line.args.cutoff)

	'''Iteration through the head, sequence generator of the FastaReader class.'''
	for h, s in head_seq:
		#cleaning sequence of newline chars and N bases
		clean_seq = s.replace('\\', "")
		cleanDubBases = clean_seq.replace('N', "")# N bases are assumed to be sequencing errors, and therefore removed from sequence, and not counted in as part of the genome length
		#calling Motif class methods
		res.lengthOfSequence(cleanDubBases)
		res.kmerComposition(cleanDubBases)
	res.addData()
	res.expectedCounts()
	res.sdScore()
	res.zScore()
	res.results()

	#sorting results list based on z-score
	sort_1 = sorted(res.res_lst, key=lambda tup: tup[3])
	#sorting results list a second time based on motif size
	sort_2 = sorted(sort_1, key=lambda k: len(k[0]), reverse=True)

	#printing of headers
	print('Motif','\t','Actual_Count','\t','Expected_Count','\t','Zscore')
	#printing of sorted results
	for i in sort_2:
		print('{0:8}\t{1:0d}\t{2:0.2f}\t{3:0.2f}'.format(i[0], i[1], i[2], i[3]))



if __name__ == "__main__":
	main()