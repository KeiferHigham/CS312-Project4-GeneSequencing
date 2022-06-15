#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		self.back_pointers = {}

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):

		if banded == True:
			ban = BandedSequence(seq1,seq2,align_length)
			ban.align_sequence()
			score = ban.score
			first_alignment = ban.alignment_one
			second_alignment = ban.alignment_two
			return {'align_cost': score, 'seqi_first100': first_alignment, 'seqj_first100': second_alignment}


		self.banded = banded
		self.MaxCharactersToAlign = align_length

		numbers = []
		self.back_pointers = {}
		upper_bound_table_row = len(seq2) + 1
		upper_bound_table_column = len(seq1) + 1
		if upper_bound_table_row > align_length:
			upper_bound_table_row = align_length + 1

		if upper_bound_table_column > align_length:
			upper_bound_table_column = align_length + 1

		E = [[0 for x in range(upper_bound_table_row)] for y in range(upper_bound_table_column)]
		back_e = [[0 for x in range(upper_bound_table_row)] for y in range(upper_bound_table_column)]

		# Iterate through the first row and first column which takes M + N time where
		# M is the number of rows and N is the number of columns
		for i in range(upper_bound_table_column):
			E[i][0] = i * 5
		for j in range(upper_bound_table_row):
			E[0][j] = j * 5

		upper_bound_row = len(seq1) + 1
		upper_bound_column = len(seq2) + 1
		if upper_bound_row > align_length:
			upper_bound_row = align_length + 1

		if upper_bound_column > align_length:
			upper_bound_column = align_length + 1

			# This will take O(MN) time where M is the number of rows and N is the number of Columns
			# during each iteration we perform multiple constant time functions that compare the different numbers
			# so peforming MN constant time operations the total complexity comes out to be O(MN)
		for i in range(1,upper_bound_row):
			for j in range(1,upper_bound_column):
				# remember that letter indexs start one index later so must subtract 1 from both i and j when using them to access letters
				numbers = self.min_number(E[i-1][j] + 5,E[i][j-1] +5,self.get_num(E[i-1][j-1],seq1[i - 1],seq2[j - 1]))
				E[i][j] = numbers[0]
				back_e[i][j] = numbers[1]

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		list_modified_words = self.modify_words(upper_bound_row - 1,upper_bound_column - 1,seq1,seq2,back_e)
		seq1 = list_modified_words[0]
		seq2 = list_modified_words[1]
		score = E[upper_bound_row - 1][upper_bound_column - 1]
		alignment1 = seq1[:100]
		alignment2 = seq2[:100]
###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
      # complexity will just be O(1) here
	def get_num(self,first_val,second_val,third_val):
		if second_val == third_val:
			return first_val -3
		return first_val + 1
	# just comparing different values so time complexity will be O(1)
	def min_number(self,first_num,second_num,third_num):
		min_num = min(first_num,second_num,third_num)
		#second is left
		# first is above
		#third is from diagnal precedence is second, first, third
		# we pass in info such as where the min number came from and how it gets to the next box
		if second_num == first_num and second_num == third_num:
			return [second_num, 2]
		if second_num == third_num and second_num == min_num:
			return [second_num, 2]
		if second_num == first_num and second_num == min_num:
			return [second_num, 2]
		if first_num == third_num and first_num == min_num:
			return [first_num, 1]
		if first_num == min_num:
			return[first_num, 1]
		if second_num == min_num:
			return[second_num, 2]
		return[third_num, 3]


	# Start from optimal cost which in this case will be bottome right corner and
	# work our way back up, the total time complexity here will be m + n which is the worst case
	# where m is the number of rows and n is the number of columns
	# could be around m if all of the characters match and we just move along the diagnal
	def modify_words(self,i,j,seq1,seq2,back_e):
		#sec1 is row string
		#sec2 is column string
		# value is the an array where 0 is value it came from
		# 1 is the index it came from
		# 2 how it came from previous box 1 from above, 2 from left, 3 from diagnal
		while i != 0 and j != 0:
			if back_e[i][j] == 1:
				column_index = j
				seq2 = seq2[:column_index] + "-" + seq2[column_index:]
				i = i - 1
			if back_e[i][j] == 2:
				row_index = i
				seq1 = seq1[:row_index] + "-" + seq1[row_index:]
				j = j - 1
			if back_e[i][j] == 3:
				i = i - 1
				j = j - 1

		return [seq1,seq2]


class BandedSequence:
	seq1 = ""
	seq2 = ""
	align_length = 0
	score = 0
	alignment_one = ""
	alignment_two = ""
	def __init__(self,seq1,seq2, align_length):
		self.seq1 = seq1
		self.seq2 = seq2
		self.align_length = align_length
		self.score = 0
		self.alignment_one = ""
		self.alignment_two = ""



	def align_sequence(self):

		upper_bound_table_row = len(self.seq1) + 1
		upper_bound_column = len(self.seq2) + 1

		self.seq1 = self.seq1[:self.align_length]
		self.seq2 = self.seq2[:self.align_length]

		if upper_bound_table_row > self.align_length :
			upper_bound_table_row = self.align_length + 1


		if upper_bound_column > self.align_length:
			self.seq2 = self.seq2[0:self.align_length]

		# complexity here will be m + 7
		# where just going to fill in the first row which is m length
		# and feel up the first column with the length of the column being 7
		E = [[0 for x in range(7)] for y in range(upper_bound_table_row)]
		back_e = [[0 for x in range(7)] for y in range(upper_bound_table_row)]


		#Here I'm just looping through the first 4 rows and adding the values
		# just like we would in the unrestricted algorithm except we don't consider
		# the top when we get to the end because there is none
		# complexity here will be around (4 *7) because 4 rows around 7 columns for each row
		for i in range(4):
			E[i][0] = i * 5

		for j in range(4):
			E[0][j] = j * 5

		for i in range(1,4):
			j = 1
			while (j - i) <= 3:
				numbers = []
				if j - i == 3:
					# making the top number infinity so we don't even consider it like it doesn't exist
					numbers = self.min_number(math.inf, E[i][j - 1] + 5,
											  self.get_num(E[i - 1][j - 1], self.seq1[i - 1], self.seq2[j - 1]))
				else:
					numbers = self.min_number(E[i - 1][j] + 5, E[i][j - 1] + 5,
											  self.get_num(E[i - 1][j - 1], self.seq1[i - 1], self.seq2[j - 1]))
				E[i][j] = numbers[0]
				back_e[i][j] = numbers[1]
				j = j + 1

		k = 7
		last = 0
		last_row = 0

		# The total time complexity will be KN
		# where K represents the width of the band
		# which in this case will be 7 because
		# we have 3 to the left and 3 to the right of each value
		# numbers get shifted here so left is still left
		# up is now diagnal and diagnal to the right is now up
		# i + j -4 helps us find the correct character for the column so we can compare
		for i in range(4, upper_bound_table_row):

			for j in range(k):

				if (i +j) - 4 < len(self.seq2):
					numbers = []
					if j == 0:
						numbers = self.min_number(E[i - 1][j + 1] + 5,math.inf,self.get_num(E[i - 1][j],self.seq1[i - 1],self.seq2[(i +j) - 4]))
					if j == 6:
						numbers = self.min_number(math.inf,E[i][j - 1] + 5,self.get_num(E[i - 1][j],self.seq1[i - 1],self.seq2[(i +j) - 4]))
					if j != 0 and j != 6:
						numbers = self.min_number(E[i - 1][j + 1] + 5,E[i][j-1] + 5,self.get_num(E[i - 1][j],self.seq1[i - 1],self.seq2[(i + j) - 4]))
					E[i][j] = numbers[0]
					back_e[i][j] = numbers[1]
					last = j
					last_row = i

		# checking to see if there's a huge discrepancy here
		# if huge difference between numbers we have to account for that
		if math.fabs(len(self.seq2) - len(self.seq1)) > 1000:
			self.score = math.inf
			self.alignment_one = "No Alignment Possible"
			self.alignment_two = "No Alignment Possible"
		else:
			aligned_words = self.align_words(self.seq1,self.seq2,last_row,last,back_e)
			self.seq1 = aligned_words[0]
			self.seq2 = aligned_words[1]
			self.score = E[last_row][last]
			self.alignment_one = self.seq1[:100]
			self.alignment_two = self.seq2[:100]


		# start at the optimal alignment which will be the last entry
		# position i j in the cost matrix
		# work our way back to entry at i = 0 j = 0
		# complexity here will be O(m + K) where m is number of rows
		# and K is the number of columns which in this case will be 7
	def align_words(self,seq1,seq2,i,j,back_e):
		non = None

		while i != 0 and j != 0:
			if i > 3:
				if back_e[i][j] == 1:
					column_index = (i + j) - 4
					seq2 = seq2[:column_index] + "-" + seq2[column_index:]
					i = i - 1
					j = j + 1
				if back_e[i][j] == 2:
					row_index = i
					seq1 = seq1[:row_index] + "-" + seq1[row_index:]
					j = j - 1
				if back_e[i][j] == 3:
					i = i - 1
			else:
				if back_e[i][j] == 1:
					column_index = j
					seq2 = seq2[:column_index] + "-" + seq2[column_index:]
					i = i - 1
				if back_e[i][j] == 2:
					row_index = i
					seq1 = seq1[:row_index] + "-" + seq1[row_index:]
					j = j - 1
				if back_e[i][j] == 3:
					i = i - 1
					j = j - 1


		return [seq1,seq2]



	# just comparing two number O(1) complexity
	def get_num(self, first_val, second_val, third_val):
			if second_val == third_val:
				return first_val - 3
			return first_val + 1

	# comparing different values to determine tie breaker O(1) complexity
	def min_number(self, first_num, second_num, third_num):
		min_num = min(first_num, second_num, third_num)
		# second is left
		# first is above
		# third is from diagnal precedence is second, first, third
		# we pass in info such as where the min number came from and how it gets to the next box
		if second_num == first_num and second_num == third_num:
			return [second_num, 2]
		if second_num == third_num and second_num == min_num:
			return [second_num, 2]
		if second_num == first_num and second_num == min_num:
			return [second_num, 2]
		if first_num == third_num and first_num == min_num:
			return [first_num, 1]
		if first_num == min_num:
			return [first_num, 1]
		if second_num == min_num:
			return [second_num, 2]
		return [third_num, 3]









