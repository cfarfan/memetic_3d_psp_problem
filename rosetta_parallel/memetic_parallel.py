#!/usr/bin/env python
# -*- coding: utf-8 -*-

import textwrap
import os
import uuid
import time
import sys
import random
import time
import copy
import pdb
import traceback
import re
import string
import datetime
import logging
from math import *
from rosetta import *
from threading import Thread, Condition, Lock, Event
from Queue import PriorityQueue

#constants
NUM_POCKETS					= 5
NUM_AGENTS 					= 10
NUM_SUP 					= 2
MAX_GEN						= 30
USE_CHI_ANGLES				= False #When true, the alg uses the chi values
CUSTOMIZE_CHI_RANGE			= False #Flag that determines if chi interval comes from the input text file (true) or from functions (false)
USE_ANGLE_RANGE				= False #When true, use the range of the angles, when false, use probabilities
PROB_RADIUS 				= 0.5
FLOAT_PRECISION 			= 2 #Precision of angles - need to be > 1
GENERATE_RANDOM_SOLUTION	= False
TEST_DIV_DIFF 				= 2.0
TEST_LS_PROB				= 0.9
TEST_LS_FACT 				= 0.9
TEST_JUMP_PROB				= 0.3
TEST_JUMP_FACT				= 0.85
TEST_TEMP_INIT				= 2000 #??
TEST_JUMP_DIST				= 180
TEST_NOIMPROVE				= 50 #reset
CROSSOVER_PROB				= 0.4
TIME_LIMIT 					= datetime.timedelta(0,0,0,0,0,6,0)
TIME_WINDOW_STATUS			= datetime.timedelta(0,0,0,0,5,0,0)
IF_RESET					= True

#Constants for angle index
PHI  = 0
PSI  = 1
CHI1 = 2
CHI2 = 3
CHI3 = 4
CHI4 = 5

#Global parameters
primary_amino_sequence		= None
secondary_amino_sequence	= None
secondary_sequence_list		= []
maxmin_angles 				= []
input_file 					= None
hist_obj 					= None
scorefxn 					= None
workers = []
receivers = []

#path to histogram dir
PATH_GLOBAL = os.getcwd()+'/'

#Secondary structure constants
NUMBER_OF_SS = 7
SS_B = 0
SS_C = 1
SS_E = 2
SS_G = 3
SS_H = 4
SS_I = 5
SS_T = 6

#Amino constants
GLY = "G"
ALA = "A"
PRO = "P"
SER = "S"
CYS = "C"
THR = "T"
VAL = "V"
ILE = "I"
LEU = "L"
ASP = "D"
ASN = "N"
HIS = "H"
PHE = "F"
TYR = "Y"
TRP = "W"
MET = "M"
GLU = "E"
GLN = "Q"
LYS = "K"
ARG = "R"
ASX = "ASX"
UNK = "UNK"
GLX = "GLX"

sigla = {'A': "ALA",
		'R': "ARG",
		'N': "ASN",
		'D': "ASP",
		'R': "ARG",
		'B': "ASX",
		'C': "CYS",
		'E': "GLU",
		'Q': "GLN",
		'Z': "GLX",
		'G': "GLY",
		'H': "HIS",
		'I': "ILE",
		'L': "LEU",
		'K': "LYS",
		'M': "MET",
		'F': "PHE",
		'P': "PRO",
		'S': "SER",
		'T': "THR",
		'W': "TRP",
		'Y': "TYR",
		'V': "VAL"}
		
siglaSS = {"0" : "B", "1" : "C", "2" : "E", "3" : "G", "4" : "H", "5" : "I", "6" : "T"}

#Errors code
ERROR_CODE_fileopening = "Error while opening "
ERROR_CODE_fileclosing = "Error while closing"
ERROR_CODE_filewriting = "Error while writing "
ERROR_CODE_PDBwriting = "Error while writing PDB file "
ERROR_CODE_emptyfile = "Error: empty file"
ERROR_CODE_lessangles = "Error: file has less angles than expected"
ERROR_CODE_incompletefile = "Error: file is incomplete"
ERROR_CODE_mysteriousresidue = "Error: mysterious amino acid residue found"

class HistogramFiles:
	def __init__(self):
		self.histogram = {ALA:[], ARG:[], ASN:[], ASP:[], ASX:[],
					CYS:[], GLN:[], GLU:[], GLX:[], GLY:[],
					HIS:[], ILE:[], LEU:[], LYS:[], MET:[],
					PHE:[], PRO:[], SER:[], THR:[], TRP:[],
					TYR:[], UNK:[], VAL:[]}
					 
		self.prob_hist = {ALA:[], ARG:[], ASN:[], ASP:[], ASX:[],
				CYS:[], GLN:[], GLU:[], GLX:[], GLY:[],
				HIS:[], ILE:[], LEU:[], LYS:[], MET:[],
				PHE:[], PRO:[], SER:[], THR:[], TRP:[],
				TYR:[], UNK:[], VAL:[]}
		
		#define files names
		name = ""
		for key in self.histogram:
			for ss in ["B", "C", "E", "G", "H", "I", "T"]:
				if key == GLY:
					name = "GLY"
				elif key == ALA:
					name = "ALA"
				elif key == PRO:
					name = "PRO"
				elif key == SER:
					name = "SER"
				elif key == CYS:
					name = "CYS"
				elif key == THR:
					name = "THR"
				elif key == VAL:
					name = "VAL"
				elif key == ILE:
					name = "ILE"
				elif key == LEU:
					name = "LEU"
				elif key == ASP:
					name = "ASP"
				elif key == ASN:
					name = "ASN"
				elif key == HIS:
					name = "HIS"
				elif key == PHE:
					name = "PHE"
				elif key == TYR:
					name = "TYR"
				elif key == TRP:
					name = "TRP"
				elif key == MET:
					name = "MET"
				elif key == GLU:
					name = "GLU"
				elif key == GLN:
					name = "GLN"
				elif key == LYS:
					name = "LYS"
				elif key == ARG:
					name = "ARG"
				self.histogram[key].append(name+"_"+ss+"_histogram"+".dat")
		
		#probs list		
		for key in self.prob_hist:
			for ss in ["B", "C", "E", "G", "H", "I", "T"]:
				self.prob_hist[key].append([])
		
	def read_histograms(self):
		i=0
		while (i < len(primary_amino_sequence)):
			try:
				print "Opening "+str(self.histogram[primary_amino_sequence[i]][secondary_sequence_list[i]])
				hist_file = open(str(""+PATH_GLOBAL+"Histogramas/PhiPsi1-1/")+self.histogram[primary_amino_sequence[i]][secondary_sequence_list[i]])
			
			except:
				print("ERROR")
				sys.exit(ERROR_CODE_fileopening+self.histogram[primary_amino_sequence[i]][secondary_sequence_list[i]])
			
			hist_line = hist_file.readline()
			while (hist_line):
				if(hist_line.startswith("#") or hist_line.strip() == ''):
					hist_line = hist_file.readline()
					continue
				hist_data = re.findall('\[[^\]]*\]|\{[^}}]*\}|\"[^\"]*\"|\S+', hist_line)
				
				try:
					#read the angles values
					hist_phi = float(hist_data[0])
					hist_psi = float(hist_data[1])
					hist_prob = float(hist_data[2])
					hist_omega = ""
					hist_chi1 = ""
					hist_chi2 = ""
					hist_chi3 = ""
					hist_chi4 = ""
					num_chi = -1
					if len(hist_data) == 8:
						hist_omega = str(hist_data[3][1:-1])
						hist_chi1 = str(hist_data[4][1:-1])
						hist_chi2 = str(hist_data[5][1:-1])
						hist_chi3 = str(hist_data[6][1:-1])
						hist_chi4 = str(hist_data[7][1:-1])
						num_chi = 4
					elif len(hist_data) == 7:
						hist_omega = str(hist_data[3][1:-1])
						hist_chi1 = str(hist_data[4][1:-1])
						hist_chi2 = str(hist_data[5][1:-1])
						hist_chi3 = str(hist_data[6][1:-1])
						num_chi = 3
					elif len(hist_data) == 6:
						hist_omega = str(hist_data[3][1:-1])
						hist_chi1 = str(hist_data[4][1:-1])
						hist_chi2 = str(hist_data[5][1:-1])
						num_chi = 2
					elif len(hist_data) == 5:
						hist_omega = str(hist_data[3][1:-1])
						hist_chi1 = str(hist_data[4][1:-1])
						num_chi = 1
					elif len(hist_data) == 4:
						hist_omega = str(hist_data[3][1:-1])
						num_chi = 0	
										
				except:
					print("ERROR")
					sys.exit(ERROR_CODE_incompletefile+" "+self.histogram[amino_acid_sequence[i]][secondary_structure_sequence_i[i]]+ " "+str(hist_data))

				try:
					if (hist_prob != 0.0):
						if num_chi == -1:
							self.prob_hist[primary_amino_sequence[i]][secondary_sequence_list[i]].append((hist_phi, hist_psi, hist_prob))
						elif num_chi == 0:
							self.prob_hist[primary_amino_sequence[i]][secondary_sequence_list[i]].append((hist_phi, hist_psi, hist_prob, hist_omega))
						elif num_chi == 1:
							self.prob_hist[primary_amino_sequence[i]][secondary_sequence_list[i]].append((hist_phi, hist_psi, hist_prob, hist_omega, hist_chi1))
						elif num_chi == 2:
							self.prob_hist[primary_amino_sequence[i]][secondary_sequence_list[i]].append((hist_phi, hist_psi, hist_prob, hist_omega, hist_chi1, hist_chi2))
						elif num_chi == 3:
							self.prob_hist[primary_amino_sequence[i]][secondary_sequence_list[i]].append((hist_phi, hist_psi, hist_prob, hist_omega, hist_chi1, hist_chi2, hist_chi3))
						elif num_chi == 4:
							self.prob_hist[primary_amino_sequence[i]][secondary_sequence_list[i]].append((hist_phi, hist_psi, hist_prob, hist_omega, hist_chi1, hist_chi2, hist_chi3, hist_chi4))																		
				except:
					print("ERROR")
					sys.exit("Line "+": "+"Failed while loading data from histogram: "+ str(amino_acid_sequence[i])+ " " + str(secondary_structure_sequence_i[i]))

				hist_line = hist_file.readline()
				
			self.prob_hist[primary_amino_sequence[i]][secondary_sequence_list[i]].sort(key = lambda x: x[2], reverse=True) #Sort the prob list from the most probably to the less
			
			try:
				hist_file.close()
			except:
				print("ERROR")
				sys.exit(ERROR_CODE_fileclosing+Histogram[amino_acid_sequence[i]][secondary_structure_sequence_i[i]])
					
			i += 1
		
		i=0
		for aminoacid in primary_amino_sequence:
			if i == 0 or i == len(primary_amino_sequence)-1: 
				#jump the first and last residue
				i += 1 
				continue
			else:
				flag3Hist = True
				flag3Hist2 = True
				flag3Hist3 = True
				hist_file1 = ""
				hist_file2 = ""
				try:
					AA_Ant = sigla[primary_amino_sequence[i-1]]
					AA = sigla[primary_amino_sequence[i]]
					AA_Prox = sigla[primary_amino_sequence[i+1]]
					
					SS_Ant = siglaSS[str(secondary_sequence_list[i-1])]
					SS = siglaSS[str(secondary_sequence_list[i])]
					SS_Prox = siglaSS[str(secondary_sequence_list[i+1])]
										
					try:
						hist_file = open(str(""+PATH_GLOBAL+"Histogramas/HistogramasFinal/")+str(AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox).lower()+"_histogram.dat","r")
						print "Opening "+str(AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox).lower()+"_histogram.dat"
						
					except:
						try:
							hist_file1 = open(str(""+PATH_GLOBAL+"Histogramas/FinalDuploVai/")+str(AA_Ant+SS_Ant+AA+SS).lower()+"_histogram.dat","r")
							print "Opening "+str(AA_Ant+SS_Ant+AA+SS).lower()+"_histogram.dat"
							flag3Hist2 = False
							hist_file1.close()
							
							hist_file2 = open(str(""+PATH_GLOBAL+"Histogramas/FinalDuploVolta/")+str(AA+SS+AA_Prox+SS_Prox).lower()+"_histogram.dat","r")
							print "Opening "+str(AA+SS+AA_Prox+SS_Prox).lower()+"_histogram.dat"
							flag3Hist3 = False
							hist_file2.close()
						except:
							flag3Hist = False

				except:
					print("ERROR")
					traceback.print_exc()
					
				if flag3Hist:
					if hist_file1!="" and hist_file2!="":
						#Histograma Duplo Vai
						hist_file1 = open(str(""+PATH_GLOBAL+"Histogramas/FinalDuploVai/")+str(AA_Ant+SS_Ant+AA+SS).lower()+"_histogram.dat","r")
						hist_file = hist_file1
						hist_line = hist_file.readline()

						while (hist_line):
							if(hist_line.startswith('#') or hist_line.strip() == ''):
								hist_line = hist_file.readline()
								continue
							hist_data = string.split(hist_line)
							try:
								hist_phi = float(hist_data[0])
								hist_psi = float(hist_data[1])
								hist_prob = float(hist_data[2])
							except:
								print("ERROR")
								traceback.print_exc()
								sys.exit(ERROR_CODE_incompletefile+" "+Histogram[aminoacid][secondary_structure_sequence_i[i]])
							try:
								if (hist_prob != 0.0):
									if flag3Hist and flag3Hist2 and flag3Hist3:
										try:
											self.prob_hist[AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox].append((hist_phi, hist_psi, hist_prob))
										except:
											self.prob_hist[AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox] = []
											self.prob_hist[AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox].append((hist_phi, hist_psi, hist_prob))
											
									elif flag3Hist and (not flag3Hist2 or not flag3Hist3):
										try:
											self.prob_hist[AA_Ant+SS_Ant+AA+SS].append((hist_phi, hist_psi, hist_prob))
										except:
											self.prob_hist[AA_Ant+SS_Ant+AA+SS] = []
											self.prob_hist[AA_Ant+SS_Ant+AA+SS].append((hist_phi, hist_psi, hist_prob))
										try:
											self.prob_hist[AA+SS+AA_Prox+SS_Prox].append((hist_phi, hist_psi, hist_prob))
										except:
											self.prob_hist[AA+SS+AA_Prox+SS_Prox] = []
											self.prob_hist[AA+SS+AA_Prox+SS_Prox].append((hist_phi, hist_psi, hist_prob))
							except:
								print("ERROR")
								traceback.print_exc()
								sys.exit("Line "+": "+"Failed while loading data from histogram: "+aminoacid+" "+str(secondary_sequence_list[i]))
								
							hist_line = hist_file.readline()
							
						#Histograma Duplo Volta
						hist_file2 = open(str(""+PATH_GLOBAL+"Histogramas/FinalDuploVolta/")+str(AA+SS+AA_Prox+SS_Prox).lower()+"_histogram.dat","r")
						hist_file = hist_file2
						hist_line = hist_file.readline()
						
						while (hist_line):
							if(hist_line.startswith('#') or hist_line.strip() == ''):
								hist_line = hist_file.readline()
								continue
							hist_data = string.split(hist_line)
							try:
								hist_phi = float(hist_data[0])
								hist_psi = float(hist_data[1])
								hist_prob = float(hist_data[2])
							except:
								print("ERROR")
								traceback.print_exc()
								sys.exit(ERROR_CODE_incompletefile+" "+Histogram[aminoacid][secondary_structure_sequence_i[i]])
							try:
								if (hist_prob != 0.0):
									if flag3Hist and flag3Hist2 and flag3Hist3:
										try:
											self.prob_hist[AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox].append((hist_phi, hist_psi, hist_prob))
										except:
											self.prob_hist[AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox] = []
											self.prob_hist[AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox].append((hist_phi, hist_psi, hist_prob))
											
									elif flag3Hist and (not flag3Hist2 or not flag3Hist3):
										try:
											self.prob_hist[AA_Ant+SS_Ant+AA+SS].append((hist_phi, hist_psi, hist_prob))
										except:
											self.prob_hist[AA_Ant+SS_Ant+AA+SS] = []
											self.prob_hist[AA_Ant+SS_Ant+AA+SS].append((hist_phi, hist_psi, hist_prob))
										try:
											self.prob_hist[AA+SS+AA_Prox+SS_Prox].append((hist_phi, hist_psi, hist_prob))
										except:
											self.prob_hist[AA+SS+AA_Prox+SS_Prox] = []
											self.prob_hist[AA+SS+AA_Prox+SS_Prox].append((hist_phi, hist_psi, hist_prob))
							except:
								print("ERROR")
								traceback.print_exc()
								sys.exit("Line "+": "+"Failed while loading data from histogram: "+aminoacid+" "+str(secondary_sequence_list[i]))
								
							hist_line = hist_file.readline()
						
					else:
						hist_line = hist_file.readline()
						
						while (hist_line):
							if(hist_line.startswith('#') or hist_line.strip() == ''):
								hist_line = hist_file.readline()
								continue
							hist_data = string.split(hist_line)
							try:
								hist_phi = float(hist_data[0])
								hist_psi = float(hist_data[1])
								hist_prob = float(hist_data[2])
							except:
								print("ERROR")
								traceback.print_exc()
								sys.exit(ERROR_CODE_incompletefile+" "+Histogram[aminoacid][secondary_structure_sequence_i[i]])
							try:
								if (hist_prob != 0.0):
									if flag3Hist and flag3Hist2 and flag3Hist3:
										try:
											self.prob_hist[AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox].append((hist_phi, hist_psi, hist_prob))
										except:
											self.prob_hist[AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox] = []
											self.prob_hist[AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox].append((hist_phi, hist_psi, hist_prob))
									elif flag3Hist and (not flag3Hist2 or not flag3Hist3):
										try:
											self.prob_hist[AA_Ant+SS_Ant+AA+SS].append((hist_phi, hist_psi, hist_prob))
										except:
											self.prob_hist[AA_Ant+SS_Ant+AA+SS] = []
											self.prob_hist[AA_Ant+SS_Ant+AA+SS].append((hist_phi, hist_psi, hist_prob))
										try:
											self.prob_hist[AA+SS+AA_Prox+SS_Prox].append((hist_phi, hist_psi, hist_prob))
										except:
											self.prob_hist[AA+SS+AA_Prox+SS_Prox] = []
											self.prob_hist[AA+SS+AA_Prox+SS_Prox].append((hist_phi, hist_psi, hist_prob))
							except:
								print("ERROR")
								traceback.print_exc()
								sys.exit("Line "+": "+"Failed while loading data from histogram: "+aminoacid+" "+str(secondary_sequence_list[i]))
								
							hist_line = hist_file.readline()
						
					if flag3Hist and flag3Hist2 and flag3Hist3:
						self.prob_hist[AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox].sort(key = lambda x: x[2], reverse=True) #Sort the prob list from the most probably to the less
					elif flag3Hist and (not flag3Hist2 or not flag3Hist3):
						self.prob_hist[AA_Ant+SS_Ant+AA+SS].sort(key = lambda x: x[2], reverse=True)
						self.prob_hist[AA+SS+AA_Prox+SS_Prox].sort(key = lambda x: x[2], reverse=True)
					try:
						hist_file.close()
					except:
						try:
							hist_file1.close()
							hist_file2.close()
						except:
							print("ERROR")
							traceback.print_exc()
							sys.exit(ERROR_CODE_fileclosing+"Error 2")
							
				i+=1
		
		#Gets the angles
		for aminoacid in primary_amino_sequence:
			line = input_file.readline()
			if(line.startswith('#') or line.strip() == ''):
				line = input_file.readline()
			if not line:
				sys.exit(ERROR_CODE_lessangles)
				
			values = string.split(line)
			try:
				Phi_min  = float(values[0])
				Phi_max  = float(values[1])
				Psi_min  = float(values[2])
				Psi_max  = float(values[3])
				if CUSTOMIZE_CHI_RANGE:
					Chi1_min = float(values[4])
					Chi1_max = float(values[5])
					Chi2_min = float(values[6])
					Chi2_max = float(values[7])
					Chi3_min = float(values[8])
					Chi3_max = float(values[9])
					Chi4_min = float(values[10])
					Chi4_max = float(values[11])
				else:
					if USE_CHI_ANGLES:
						Chi1_min = self.min_rot_chi(sigla[aminoacid], 1)
						Chi1_max = self.max_rot_chi(sigla[aminoacid], 1)
						Chi2_min = self.min_rot_chi(sigla[aminoacid], 2)
						Chi2_max = self.max_rot_chi(sigla[aminoacid], 2)
						Chi3_min = self.min_rot_chi(sigla[aminoacid], 3)
						Chi3_max = self.max_rot_chi(sigla[aminoacid], 3)
						Chi4_min = self.min_rot_chi(sigla[aminoacid], 4)
						Chi4_max = self.max_rot_chi(sigla[aminoacid], 4)
					else:
						Chi1_min = 0.0
						Chi1_max = 0.0
						Chi2_min = 0.0
						Chi2_max = 0.0
						Chi3_min = 0.0
						Chi3_max = 0.0
						Chi4_min = 0.0
						Chi4_max = 0.0

			except:
				print("ERROR")
				traceback.print_exc()
				sys.exit(ERROR_CODE_incompletefile)
				
			if USE_ANGLE_RANGE:
				maxmin = [Phi_min, Phi_max, Psi_min, Psi_max,
						  Chi1_min, Chi1_max, Chi2_min, Chi2_max,
						  Chi3_min, Chi3_max, Chi4_min, Chi4_max]
			else:
				#Palliative measure for the local search when using probabilities instead of defined angle ranges
				if(hist_phi == -180.0):
					p_phi_min = hist_phi
				else:
					p_phi_min = hist_phi - PROB_RADIUS
				if(hist_phi == 180.0):
					p_phi_max = hist_phi
				else:
					p_phi_max = hist_phi + PROB_RADIUS
				#Adjust the psi range borders, can't be greater than 180 or smaller than -180
				if(hist_psi == -180.0):
					p_psi_min = hist_psi
				else:
					p_psi_min = hist_psi - PROB_RADIUS
				if(hist_psi == 180.0):
					p_psi_max = hist_psi
				else:
					p_psi_max = hist_psi + PROB_RADIUS
					
				maxmin = [-180.0, 180.0, -180.0, 180.0,
						  Chi1_min, Chi1_max, Chi2_min, Chi2_max,
						  Chi3_min, Chi3_max, Chi4_min, Chi4_max]
			maxmin_angles.append(copy.deepcopy(maxmin))

		try:
			input_file.close()
		except:
			print("ERROR")
			traceback.print_exc()
			sys.exit(ERROR_CODE_fileclosing+input_file_name)

		#print "Done"
	
	def use_histogram(self, maxmin_angles, prob_list, prob_radius, prob2_list, name):
		floatprecision=FLOAT_PRECISION

		if prob2_list == []:
			luck = random.uniform(0.0, 1.0)
			edge = 0.0
			for probs in prob_list:
				if(luck <= probs[2] + edge):
					#Adjust the phi and psi range borders, can't be greater than 180 or smaller than -180
					p_phi_min = probs[0] - prob_radius
					p_phi_min = -180.0 if p_phi_min < -180.0 else 180.0 if p_phi_min > 180.0 else p_phi_min
					p_phi_max = probs[0] + prob_radius
					p_phi_max = -180.0 if p_phi_max < -180.0 else 180.0 if p_phi_max > 180.0 else p_phi_max
					p_psi_min = probs[1] - prob_radius
					p_psi_min = -180.0 if p_psi_min < -180.0 else 180.0 if p_psi_min > 180.0 else p_psi_min
					p_psi_max = probs[1] + prob_radius
					p_psi_max = -180.0 if p_psi_max < -180.0 else 180.0 if p_psi_max > 180.0 else p_psi_max

					angles = self.get_angle_chis(probs)
					aa_angles = [round(random.uniform(p_phi_min, p_phi_max),floatprecision), round(random.uniform(p_psi_min, p_psi_max),floatprecision)] + angles
					break
				else:
					edge = edge + probs[2]
					p_backup = probs
			else: #for
				p_phi_min = p_backup[0] - prob_radius
				p_phi_min = -180.0 if p_phi_min < -180.0 else 180.0 if p_phi_min > 180.0 else p_phi_min
				p_phi_max = p_backup[0] + prob_radius
				p_phi_max = -180.0 if p_phi_max < -180.0 else 180.0 if p_phi_max > 180.0 else p_phi_max
				p_psi_min = p_backup[1] - prob_radius
				p_psi_min = -180.0 if p_psi_min < -180.0 else 180.0 if p_psi_min > 180.0 else p_psi_min
				p_psi_max = p_backup[1] + prob_radius
				p_psi_max = -180.0 if p_psi_max < -180.0 else 180.0 if p_psi_max > 180.0 else p_psi_max

				angles = self.get_angle_chis(p_backup)
				aa_angles = [round(random.uniform(p_phi_min, p_phi_max),floatprecision), round(random.uniform(p_psi_min, p_psi_max),floatprecision)] + angles
			return aa_angles
		else:
			luck = random.uniform(0.0, 1.0)
			edge = 0.0
			for probs in prob_list:
				if(luck <= probs[2] + edge):
					p_phi_min = probs[0] - prob_radius
					p_phi_min = -180.0 if p_phi_min < -180.0 else 180.0 if p_phi_min > 180.0 else p_phi_min
					p_phi_max = probs[0] + prob_radius
					p_phi_max = -180.0 if p_phi_max < -180.0 else 180.0 if p_phi_max > 180.0 else p_phi_max
					p_psi_min = probs[1] - prob_radius
					p_psi_min = -180.0 if p_psi_min < -180.0 else 180.0 if p_psi_min > 180.0 else p_psi_min
					p_psi_max = probs[1] + prob_radius
					p_psi_max = -180.0 if p_psi_max < -180.0 else 180.0 if p_psi_max > 180.0 else p_psi_max
					testboo = False
					angles = []
					for probs2 in prob2_list:
						if probs[0] == probs2[0]:
							if probs[1] == probs2[1]:
								angles = self.get_angle_chis(probs2)
								#print "Achou"
								testboo = True
								aa_angles = [round(random.uniform(p_phi_min, p_phi_max),floatprecision), round(random.uniform(p_psi_min, p_psi_max),floatprecision)] + angles
								break
					if not testboo:
						#print "Ferrou"
						#print name
						#print probs
						luck = random.uniform(0.0, 1.0)
						edge = 0.0
						for probs in prob2_list:
							if(luck <= probs[2] + edge):
								#Adjust the phi and psi range borders, can't be greater than 180 or smaller than -180
								p_phi_min = probs[0] - prob_radius
								p_phi_min = -180.0 if p_phi_min < -180.0 else 180.0 if p_phi_min > 180.0 else p_phi_min
								p_phi_max = probs[0] + prob_radius
								p_phi_max = -180.0 if p_phi_max < -180.0 else 180.0 if p_phi_max > 180.0 else p_phi_max
								p_psi_min = probs[1] - prob_radius
								p_psi_min = -180.0 if p_psi_min < -180.0 else 180.0 if p_psi_min > 180.0 else p_psi_min
								p_psi_max = probs[1] + prob_radius
								p_psi_max = -180.0 if p_psi_max < -180.0 else 180.0 if p_psi_max > 180.0 else p_psi_max

								angles = self.get_angle_chis(probs)
								aa_angles = [round(random.uniform(p_phi_min, p_phi_max),floatprecision), round(random.uniform(p_psi_min, p_psi_max),floatprecision)] + angles
								break
							else:
								edge = edge + probs[2]
								p_backup = probs
						else:
							p_phi_min = p_backup[0] - prob_radius
							p_phi_min = -180.0 if p_phi_min < -180.0 else 180.0 if p_phi_min > 180.0 else p_phi_min
							p_phi_max = p_backup[0] + prob_radius
							p_phi_max = -180.0 if p_phi_max < -180.0 else 180.0 if p_phi_max > 180.0 else p_phi_max
							p_psi_min = p_backup[1] - prob_radius
							p_psi_min = -180.0 if p_psi_min < -180.0 else 180.0 if p_psi_min > 180.0 else p_psi_min
							p_psi_max = p_backup[1] + prob_radius
							p_psi_max = -180.0 if p_psi_max < -180.0 else 180.0 if p_psi_max > 180.0 else p_psi_max

							angles = self.get_angle_chis(p_backup)
							aa_angles = [round(random.uniform(p_phi_min, p_phi_max),floatprecision), round(random.uniform   (p_psi_min, p_psi_max),floatprecision)] + angles

						return aa_angles
					break
				else:
					edge = edge + probs[2]
					p_backup = probs
			else:
				#print "Entrou Aqui-P_backup"
				p_phi_min = p_backup[0] - prob_radius
				p_phi_min = -180.0 if p_phi_min < -180.0 else 180.0 if p_phi_min > 180.0 else p_phi_min
				p_phi_max = p_backup[0] + prob_radius
				p_phi_max = -180.0 if p_phi_max < -180.0 else 180.0 if p_phi_max > 180.0 else p_phi_max
				p_psi_min = p_backup[1] - prob_radius
				p_psi_min = -180.0 if p_psi_min < -180.0 else 180.0 if p_psi_min > 180.0 else p_psi_min
				p_psi_max = p_backup[1] + prob_radius
				p_psi_max = -180.0 if p_psi_max < -180.0 else 180.0 if p_psi_max > 180.0 else p_psi_max
				testboo = False
				angles = []
				for probs2 in prob2_list:
					if p_backup[0] == probs2[0]:
						if p_backup[1] == probs2[1]:
							angles = self.get_angle_chis(probs2)
							#print "Achou2"
							testboo = True
							aa_angles = [round(random.uniform(p_phi_min, p_phi_max),floatprecision), round(random.uniform(p_psi_min, p_psi_max),floatprecision)] + angles
							break
							
				if not testboo:
					#print "Ferrou2"
					#print name
					#print probs
					luck = random.uniform(0.0, 1.0)
					edge = 0.0
					for probs in prob2_list:
						if(luck <= probs[2] + edge):
							#Adjust the phi and psi range borders, can't be greater than 180 or smaller than -180
							p_phi_min = probs[0] - prob_radius
							p_phi_min = -180.0 if p_phi_min < -180.0 else 180.0 if p_phi_min > 180.0 else p_phi_min
							p_phi_max = probs[0] + prob_radius
							p_phi_max = -180.0 if p_phi_max < -180.0 else 180.0 if p_phi_max > 180.0 else p_phi_max
							p_psi_min = probs[1] - prob_radius
							p_psi_min = -180.0 if p_psi_min < -180.0 else 180.0 if p_psi_min > 180.0 else p_psi_min
							p_psi_max = probs[1] + prob_radius
							p_psi_max = -180.0 if p_psi_max < -180.0 else 180.0 if p_psi_max > 180.0 else p_psi_max

							angles = self.get_angle_chis(probs)
							aa_angles = [round(random.uniform(p_phi_min, p_phi_max),floatprecision), round(random.uniform(p_psi_min, p_psi_max),floatprecision)] + angles
							break
						else:
							edge = edge + probs[2]
							p_backup = probs
					else:
						p_phi_min = p_backup[0] - prob_radius
						p_phi_min = -180.0 if p_phi_min < -180.0 else 180.0 if p_phi_min > 180.0 else p_phi_min
						p_phi_max = p_backup[0] + prob_radius
						p_phi_max = -180.0 if p_phi_max < -180.0 else 180.0 if p_phi_max > 180.0 else p_phi_max
						p_psi_min = p_backup[1] - prob_radius
						p_psi_min = -180.0 if p_psi_min < -180.0 else 180.0 if p_psi_min > 180.0 else p_psi_min
						p_psi_max = p_backup[1] + prob_radius
						p_psi_max = -180.0 if p_psi_max < -180.0 else 180.0 if p_psi_max > 180.0 else p_psi_max

						angles = self.get_angle_chis(p_backup)
						aa_angles = [round(random.uniform(p_phi_min, p_phi_max),floatprecision), round(random.uniform   (p_psi_min, p_psi_max),floatprecision)] + angles
					return aa_angles
			return aa_angles
	
	def get_angle_chis(self, probs):
		prob_radius = 0.5
		angles = []
		for i in range(3, len(probs)):
			aux = probs[i]
			anglesc = re.findall(r"[^[]*\[([^]]*)\]", aux)
			x1 = random.randint(0,len(anglesc)-1)
			x2 = random.randint(0,len(anglesc[x1].split(", "))-1)
			angle = float(anglesc[x1].split(", ")[x2])
			p_chi_min = angle - prob_radius
			p_chi_min = -180.0 if p_chi_min < -180.0 else 180.0 if p_chi_min > 180.0 else p_chi_min
			p_chi_max = angle + prob_radius
			p_chi_max = -180.0 if p_chi_max < -180.0 else 180.0 if p_chi_max > 180.0 else p_chi_max
			if i == 3:
				angles.append(angle)
			else:
				angles.append(round(random.uniform(p_chi_min, p_chi_max),FLOAT_PRECISION))
		return angles

	def min_rot_chi(self, amino, chi):
		if (amino == "SER"):
			if (chi == 1):
				return -39.0974526602
			else:
				return 0.0
	 
		if (amino == "CYS")or(amino == "CYX"):
			if (chi == 1): 
				return -157.751341619 
			else:
				return 0.0
				
		if (amino == "THR"): #only chi1
			if (chi == 1):
				return -153.24819924
			else:
				return 0.0
				
		if (amino == "VAL"):
			if (chi == 1): 
				return -37.1805850827 
			else:
				return 0.0
				
		if (amino == "ILE"): 
			if (chi == 1): 
				return -150.62880823
			if (chi == 2): 
				return -41.9536560708
			else:
				return 0.0
	
		if (amino == "LEU"):
			if (chi == 1): 
				return -158.45718524
			if (chi == 2): 
				return -41.1651784582 
			else:
				return 0.0

		if (amino == "ASP"):
			if (chi == 1):
				return -154.759146052
			if (chi == 2):
				return -46.0141277967
			else:
				return 0.0

		if (amino == "ASN"):
			if (chi == 1): 
				return -154.999658218
			if (chi == 2): 
				return -78.6804677964
			else:
				return 0.0

		if (amino == "HIS") or (amino == "HID"):
			if (chi == 1): 
				return -156.435615538 
			if (chi == 2):
				return -96.5817995178 
			else:
				return 0.0

		if (amino == "PHE"):
			if (chi == 1): 
				return -158.224844705
			if (chi == 2): 
				return -1.89539086499
			else:
				return 0.0

		if (amino == "TYR"):
			if (chi == 1): 
				return -158.518614718
			if (chi == 2): 
				return -3.56448935905
			else:
				return 0.0

		if (amino == "TRP"):
			if (chi == 1): 
				return -159.333604703
			if (chi == 2): 
				return -77.5000361533
			else:
				return 0.0

		if (amino == "MET"):
			if (chi == 1): 
				return -159.535936996
			if (chi == 2):
				return -95.9804879414
			if (chi == 3):
				return -113.617738242
			else:
				return 0.0

		if (amino == "GLU"):
			if (chi == 1): 
				return -154.3304485
			if (chi == 2): 
				return -114.381414072
			if (chi == 3): 
				return -46.2622279015
			else:
				return 0.0

		if (amino == "GLN"):
			if (chi == 1): 
				return -156.107286967
			if (chi == 2): 
				return -130.31218703
			if (chi == 3): 
				return -107.208042715
			else:
				return 0.0

		if (amino == "LYS"):
			if (chi == 1): 
				return -155.042685651
			if (chi == 2):
				return -131.521384644
			if (chi == 3):
				return -131.389541682
			if (chi == 4):
				return -131.389541682
	
		if (amino == "ARG"):
			if (chi == 1):
				return -151.830993739
			if (chi == 2): 
				return -105.206308323
			if (chi == 3): 
				return -136.745709869
			if (chi == 4): 
				return -136.745709869
		
		if (amino == "PRO") or (amino == "ALA") or (amino == "GLY"):
			if (chi == 1):
				return 0.0
			if (chi == 2): 
				return 0.0
			if (chi == 3): 
				return 0.0
			if (chi == 4): 
				return 0.0
				
	def max_rot_chi(self, amino, chi):
		if (amino == "SER"):
			if (chi == 1):
				return 159.364119327
			else:
				return 0.0
	 
		if (amino == "CYS")or(amino == "CYX"):
			if (chi == 1): 
				return 39.1513416194
			else:
				return 0.0

		if (amino == "THR"): #only chi1
			if (chi == 1):
				return 38.1815325735
			else:
				return 0.0
				
		if (amino == "VAL"):
			if (chi == 1): 
				return 156.980585083
			else:
				return 0.0

		if (amino == "ILE"): 
			if (chi == 1): 
				return 36.9621415634
			if (chi == 2): 
				return 158.553656071
			else:
				return 0.0
	
		if (amino == "LEU"):
			if (chi == 1): 
				return 34.5682963508
			if (chi == 2): 
				return 153.698511792
			else:
				return 0.0

		if (amino == "ASP"):
			if (chi == 1):
				return 35.1147016071
			if (chi == 2):
				return 46.9696833523
			else:
				return 0.0

		if (amino == "ASN"):
			if (chi == 1): 
				return 36.4441026627
			if (chi == 2): 
				return 102.880467796
			else:
				return 0.0

		if (amino == "HIS") or (amino == "HID"):
			if (chi == 1): 
				return 39.0133933155
			if (chi == 2):
				return 134.60402174
			else:
				return 0.0

		if (amino == "PHE"):
			if (chi == 1): 
				return 37.1915113716
			if (chi == 2): 
				return 90.4287241983
			else:
				return 0.0

		if (amino == "TYR"):
			if (chi == 1): 
				return 40.6852813843
			if (chi == 2): 
				return 90.3311560257
			else:
				return 0.0
				
		if (amino == "TRP"):
			if (chi == 1): 
				return 37.20027137
			if (chi == 2): 
				return 76.1667028199
			else:
				return 0.0

		if (amino == "MET"):
			if (chi == 1): 
				return 38.6100110704
			if (chi == 2):
				return 136.936043497
			if (chi == 3):
				return 127.810330834
			else:
				return 0.0

		if (amino == "GLU"):
			if (chi == 1): 
				return 36.0193373885
			if (chi == 2): 
				return 127.166599257
			if (chi == 3): 
				return 45.906672346
			else:
				return 0.0

		if (amino == "GLN"):
			if (chi == 1): 
				return 38.251731411
			if (chi == 2): 
				return 110.023298141
			if (chi == 3): 
				return 78.3524871589
			else:
				return 0.0

		if (amino == "LYS"):
			if (chi == 1): 
				return 45.9834263916
			if (chi == 2):
				return 107.881878471
			if (chi == 3):
				return 109.224109583
			if (chi == 4):
				return 109.224109583
	
		if (amino == "ARG"):
			if (chi == 1):
				return 55.4581542333
			if (chi == 2): 
				return 135.482851533
			if (chi == 3): 
				return 95.5111419674
			if (chi == 4): 
				return 95.5111419674
		
		if (amino == "PRO") or (amino == "ALA") or (amino == "GLY"):
			if (chi == 1):
				return 0.0
			if (chi == 2): 
				return 0.0
			if (chi == 3): 
				return 0.0
			if (chi == 4): 
				return 0.0

class Agent:
	def __init__(self, id):
		self.id 			= id
		self.id_supporters 	= [id*NUM_SUP+i for i in range(1,NUM_SUP+1) if id*NUM_SUP+i < NUM_AGENTS] if id*NUM_SUP+1 < NUM_AGENTS else None
		self.id_leader		= int((id-1)/NUM_SUP) if id > 0 else None
		self.pockets 		= [None for i in range(0,NUM_POCKETS)]
		self.pocket_leader	= None
		self.current 		= Solution()
		self.queue_leader 	= PriorityQueue(1)
		self.queue_pocket 	= PriorityQueue()
		self.generation 	= 0
		self.restarts		= 0
		self.time_ls		= datetime.timedelta(0)
		self.time_send		= datetime.timedelta(0)
		self.time_receive	= datetime.timedelta(0)
		self.trx_send		= 0
		self.trx_receive	= 0
		self.event_restart	= Event()
		self.status_log		= []
		
	def __str__(self):
		return  textwrap.dedent("""\
				Agent %d:
				--- leader: %s
				--- supporters: %s
				--- current_energy: %s
				--- pockets_energy: %s
				--- pocket_leader: %s""" % (
					self.id,
					str(self.id_leader),
					str(self.id_supporters),
					str(self.current),
					str(self.pockets),
					str(self.pocket_leader)
					)
				)

	def status_log_append(self, time):
		self.status_log.append([str(self.pockets[0]),self.trx_send,self.trx_receive,self.time_send.total_seconds(),self.time_receive.total_seconds(),self.generation,self.restarts, time.total_seconds()])

	def status_write(self, status_log_file_path):
		f = open(status_log_file_path, 'w')
		for line in self.status_log:
			for item in line:
				f.write(str(item)+'\t')
			f.write('\n')
		f.close()

	def diversity(self, agent_pocket_1, agent_pocket_2):
		#this function isnt considering the first and last residue
		aux=0
		primary_len=len(agent_pocket_1.pose.sequence())
		
		for i in range(1, primary_len):
			if (i+1)!=primary_len:
				#avalia apenas os angulos phi e psi
				for j in range(2):
					if j==0:
						#phi
						#print 'phi: %s,%s' % (agent_pocket_1.pose.phi(i+1), agent_pocket_2.pose.phi(i+1))
						diff = abs(agent_pocket_1.pose.phi(i+1) - agent_pocket_2.pose.phi(i+1))
						if diff > 180.0:
							diff = 360 - diff;
					else:
						#psi
						#print 'psi: %s,%s' % (agent_pocket_1.pose.psi(i+1), agent_pocket_2.pose.psi(i+1))
						diff = abs(agent_pocket_1.pose.psi(i+1) - agent_pocket_2.pose.psi(i+1))
						if diff > 180.0:
							diff = 360 - diff			
					aux += diff
		
		return aux/(primary_len-2)
		# result = CA_rmsd(agent_pocket_1.pose,agent_pocket_2.pose,3,primary_len-2)
		# return result

	def update(self, solution=None):
		if solution == None:
			solution = self.current
		pocket_worst = -1
		div_flag = True
		i = 0
		div_parameter = len(solution.pose.sequence())*6.0
		# div_parameter = 7
		for pocket in self.pockets:
			if pocket == None:
				if pocket_worst == -1 or div_flag:
					pocket_worst = i
				break
			else:
				# print '\033[33m ----- id: %d gen:%d --- diversity %f\033[39m' % (self.id, self.generation, div)
				if solution.energy_value <= pocket.energy_value:
					div = self.diversity(solution,pocket)
					# if div > TEST_DIV_DIFF and div_flag:
					if div > div_parameter and div_flag:
						pocket_worst = i
					# elif div < TEST_DIV_DIFF and div > 0:
					elif div < div_parameter and div > 0:
						div_flag = False
						pocket_worst = i
					else:
						pocket_worst = -1
						break
			i +=1

		if pocket_worst != -1:
			self.pockets[pocket_worst] = copy.deepcopy(solution)
			self.pockets.sort()
			#update executed
			return True
		#nothing to update
		return False

	def crossover(self, agent1_pocket, agent2_pocket, cross_prob):
		for i in range (len(agent1_pocket.pose.sequence())):
			random.seed()
			if random.random()<=cross_prob: #leader_agent
				phi = agent1_pocket.pose.phi(i+1)
				self.current.pose.set_phi(i+1,phi)
				psi = agent1_pocket.pose.psi(i+1)
				self.current.pose.set_psi(i+1,psi)
				omega = agent1_pocket.pose.omega(i+1)
				self.current.pose.set_omega(i+1,omega)
				#set chi angles
				for k in range(self.current.pose.residue(i+1).nchi()):
					try:
						chi = agent1_pocket.pose.chi(k+1, i+1) 
						self.current.pose.set_chi(k+1, i+1, chi)
					except:
						break;
							
			else:
				phi=agent2_pocket.pose.phi(i+1)
				self.current.pose.set_phi(i+1,phi)
				psi=agent2_pocket.pose.psi(i+1)
				self.current.pose.set_psi(i+1,psi)
				omega=agent2_pocket.pose.omega(i+1)
				self.current.pose.set_omega(i+1,omega)
				#set chi angles
				for k in range(self.current.pose.residue(i+1).nchi()):
					try:
						chi = agent2_pocket.pose.chi(k+1, i+1)
						self.current.pose.set_chi(k+1, i+1, chi)
					except:
						break;
				
		self.current.calculate_energy()
		#print 'New energy value from crossover: '+str(self.current.energy_value)

	def point_center(self, coord):
		int_aux = coord
		float_aux=int_aux
		if( int_aux > 0):
			float_aux = int_aux + 0.5

		if( int_aux < 0):
			float_aux = int_aux - 0.5
		
		return float_aux

	def local_search(self, prob_ls, prob_jump, radius_jump, hist):
		#LS only in agent's current pocket
		for i in range (1, len(self.current.pose.sequence())):
			if (i+1)!=len(self.current.pose.sequence()):
				name_res=self.current.pose.residue(i+1).name3()
				random.seed()
				if random.random()<=prob_ls:
					random.seed()
					control=0
					
					if random.random()<=prob_jump:
						#JUMP
						#Histogram based solutions
						aa_angles = []				
						try:
							AA_Ant = sigla[primary_amino_sequence[i-1]]
							AA = sigla[primary_amino_sequence[i]]
							AA_Prox = sigla[primary_amino_sequence[i+1]]
							SS_Ant = siglaSS[str(secondary_sequence_list[i-1])]
							SS = siglaSS[str(secondary_sequence_list[i])]
							SS_Prox = siglaSS[str(secondary_sequence_list[i+1])]
								
						except:
							print("ERROR")
							traceback.print_exc()
							sys.exit(ERROR_CODE_fileopening+" Error 3")
							
						if USE_ANGLE_RANGE:
							#Use ranges of max and min angles
							aa_angles = [round(random.uniform(maxmin_angles[i][0],maxmin_angles[i][1]),floatprecision),#PHI
										 round(random.uniform(maxmin_angles[i][2],maxmin_angles[i][3]),floatprecision),#PSI
										 round(random.uniform(maxmin_angles[i][4],maxmin_angles[i][5]),floatprecision),#CHI1
										 round(random.uniform(maxmin_angles[i][6],maxmin_angles[i][7]),floatprecision),#CHI2
										 round(random.uniform(maxmin_angles[i][8],maxmin_angles[i][9]),floatprecision),#CHI3
										 round(random.uniform(maxmin_angles[i][10],maxmin_angles[i][11]),floatprecision)]#CHI4
						else:
							#Use histogram
							proba = []
							proba2 = []
							name = ""
							try:
								proba = hist.prob_hist[AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox]
								proba2 = hist.prob_hist[primary_amino_sequence[i-1]][secondary_sequence_list[i]]
								name = AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox
							except:
								try:
									#50% for each combination
									if (random.randint(1,10) <=5):
										proba = hist.prob_hist[AA_Ant+SS_Ant+AA+SS]
										proba2 = hist.prob_hist[primary_amino_sequence[i-1]][secondary_sequence_list[i]]
										name = AA_Ant+SS_Ant+AA+SS
									else:
										proba = hist.prob_hist[AA+SS+AA_Prox+SS_Prox]
										proba2 = hist.prob_hist[primary_amino_sequence[i-1]][secondary_sequence_list[i]]
										name = AA+SS+AA_Prox+SS_Prox
								except:
									proba = hist.prob_hist[primary_amino_sequence[i-1]][secondary_sequence_list[i]]
									name = AA+SS
							
							aa_angles = hist.use_histogram(maxmin_angles[i], proba, PROB_RADIUS, proba2, name)

						#baseado nos dados do histograma
						phi = aa_angles[0] + random.random()	
						psi = aa_angles[1] + random.random()	
						
						#manhattan distance
						int_x = abs(self.current.pose.phi(i+1) - phi)
						int_y = abs(self.current.pose.psi(i+1) - psi)
						if(int_x > int_y):
							distance = int_x
						else:
							distance = int_y
						
						if (distance < radius_jump) and (radius_jump > 1.0):	
								
							self.current.pose.set_phi(i+1, phi)
							self.current.pose.set_psi(i+1, psi)
							self.current.calculate_energy()
					
					best_energy=self.current.energy_value
					energy_temp=best_energy
					
					backup_phi = self.current.pose.phi(i+1)
					backup_psi = self.current.pose.psi(i+1)
					backup_chi_list=[] #list of chi angles
					for k in range(self.current.pose.residue(i+1).nchi()):
						try:
							backup_chi_list.append(self.current.pose.chi(k+1, i+1))
						except:
							break;
								
					#phi
					bit_start = -1.0
					bit_end = 1.0

					if self.id_leader != None:
						temperature = TEST_TEMP_INIT #700?
					else:
						temperature = TEST_TEMP_INIT*2

					while(temperature > 0.1):
						bit = random.uniform(bit_start, bit_end)
						phi = self.current.pose.phi(i+1)
						phi_ant = copy.copy(phi) #valor anterior de phi
						
						phi += bit
						if (phi > self.point_center(backup_phi) - 0.5) and (phi < self.point_center(backup_phi) + 0.5):
							#alteracao temporaria
							self.current.pose.set_phi(i+1,phi)
							self.current.calculate_energy()
							energy_temp=self.current.energy_value
							if(energy_temp < best_energy) or (random.random() <= exp((best_energy-energy_temp)/temperature)):
								best_energy = energy_temp
								control = control - 1
							else:
								self.current.pose.set_phi(i+1,phi_ant)
								self.current.calculate_energy()
								phi = phi_ant
						else:
							phi=phi_ant
						temperature = temperature*TEST_LS_FACT #0.98
					
					#psi
					if self.id_leader != None:
						temperature = TEST_TEMP_INIT #700?
					else:
						temperature = TEST_TEMP_INIT*2

					while(temperature > 0.1):
						bit = random.uniform(bit_start, bit_end)
						psi = self.current.pose.psi(i+1)
						psi_ant=copy.copy(psi)
						
						psi += bit
						if (psi > self.point_center(backup_psi) - 0.5) and (psi < self.point_center(backup_psi) + 0.5):
							#alteracao temp
							self.current.pose.set_psi(i+1,psi)
							self.current.calculate_energy()
							energy_temp=self.current.energy_value
							if(energy_temp < best_energy) or (random.random() <= exp((best_energy-energy_temp)/temperature)):
								best_energy = energy_temp
								control = control - 1
							else:
								self.current.pose.set_psi(i+1,psi_ant)
								self.current.calculate_energy()
								psi = psi_ant
						else:
							psi=psi_ant
						temperature = temperature*TEST_LS_FACT
					
					#chi angles
					for n_chi in range(self.current.pose.residue(i+1).nchi()):
						try:
							control = 0
							best_energy = self.current.energy_value
							energy_temp = best_energy
							while(control != 9991):
								bit = random.random()
								chi_value=self.current.pose.chi(n_chi+1, i+1)
								chi_ant=copy.copy(chi_value)
								if control<=0:
									chi_value += bit
									if (chi_value > hist.min_rot_chi(name_res, n_chi+1)) and (chi_value < hist.max_rot_chi(name_res, n_chi+1)):
										self.current.pose.set_chi(n_chi+1, i+1, chi_value)
										self.current.calculate_energy()
										energy_temp=self.current.energy_value
										if energy_temp<best_energy:
											best_energy=energy_temp
											control -= 1
										else:
											if control==0:
												chi_value=backup_chi_list[n_chi]
												self.current.pose.set_chi(n_chi+1, i+1, chi_value)
												self.current.calculate_energy()
												control=999
											else:
												if control<0:
													chi_value=chi_ant
													self.current.pose.set_chi(n_chi+1, i+1, chi_value)
													self.current.calculate_energy()
													control=9991
									else:
										control=999
								if control==999:
									chi_value=chi_ant
									chi_value -= bit
									if (chi_value > hist.min_rot_chi(name_res, n_chi+1)) and (chi_value < hist.max_rot_chi(name_res, n_chi+1)):
										self.current.pose.set_chi(n_chi+1, i+1, chi_value)
										self.current.calculate_energy()
										energy_temp=self.current.energy_value
										if energy_temp<best_energy:
											best_energy=energy_temp
											
										else:
											self.current.pose.set_chi(n_chi+1, i+1, chi_ant)
											self.current.calculate_energy()
											chi_value=chi_ant
											control=9991
									else:
										control=9991
						except:
							break

	# def mutation(self):
	# 	pass

class Solution:
	def __init__(self, pose = None):
		if pose == None:
			self.pose = Pose()
			self.energy_value=None
		else:
			self.pose=copy.deepcopy(pose)
			self.calculate_energy()
		

	def __str__(self):
		return textwrap.dedent('%s' % (str(self.energy_value)))
		# return textwrap.dedent('%s\033[35m <%s>\033[39m' % (str(self.energy_value),str(hex(id(self)))))

	def __repr__(self):
		return textwrap.dedent('%s' % (str(self.energy_value)))
		# return textwrap.dedent('%s\033[35m <%s>\033[39m' % (str(self.energy_value),str(hex(id(self)))))
		# return textwrap.dedent('Solution<%s>: %s'  % (
		# 			str(id(self)) if not None else 'None',
		# 			str(self.energy_value)
		# 			))

	def __cmp__(self,other):
		if other == None or other.energy_value > self.energy_value:
			return -1
		elif other.energy_value == self.energy_value:
			return 0
		elif other.energy_value < self.energy_value:
			return 1

	def calculate_energy(self):
		self.energy_value = scorefxn(self.pose)

	def init_solution(self,hist):
		self.pose = Pose()
		self.pose = pose_from_sequence(primary_amino_sequence, 'fa_standard')
		self.generate_first_angles(hist)
		self.calculate_energy()

	def generate_first_angles(self, hist):
		if GENERATE_RANDOM_SOLUTION:
			random.seed()
			
			#first psi
			psi = random.uniform(-180,180)
			#set_psi(id_res, angle)
			self.pose.set_psi(1,psi)
			self.pose.set_phi(1,0)
			self.pose.set_omega(1, 180)		
			#set chi angles
			#name_res=self.pose.residue(1).name3()
			for k in range(self.pose.residue(1).nchi()):
				#min_chi=ag.min_rot_chi(name_res, k+1)
				#max_chi=ag.max_rot_chi(name_res, k+1)	
				#chi=random.uniform(min_chi, max_chi)
				chi=random.uniform(-180,180)
				self.pose.set_chi(k+1, 1, chi)
			
			for i in range(1, len(self.pose.sequence())):
				#last amino - jp psi
				#last amino - jp psi
				if (i+1) == len(self.pose.sequence()):
					break
				
				phi = random.uniform(-180,180)
				self.pose.set_phi(i+1,phi)
				psi = random.uniform(-180,180)
				self.pose.set_psi(i+1,psi)
				self.pose.set_omega(1, 180)
				
				#set chi angles
				#name_res=self.pose.residue(i+1).name3()			
				for k in range(self.pose.residue(i+1).nchi()):			
					#min_chi=ag.min_rot_chi(name_res, k+1)
					#max_chi=ag.max_rot_chi(name_res, k+1)	
					#chi=random.uniform(min_chi, max_chi)
					chi = random.uniform(-180,180)
					self.pose.set_chi(k+1, i+1, chi)

			#last phi
			phi = random.uniform(-180,180)
			self.pose.set_phi(len(self.pose.sequence()),phi)
			self.pose.set_psi(len(self.pose.sequence()),0)
			self.pose.set_omega(len(self.pose.sequence()), 180)
			#set chi angles
			#name_res=self.pose.residue(len(self.pose.sequence())).name3()			
			for k in range(self.pose.residue(len(self.pose.sequence())).nchi()):
				#min_chi=ag.min_rot_chi(name_res, k+1)
				#max_chi=ag.max_rot_chi(name_res, k+1)	
				#chi=random.uniform(min_chi, max_chi)
				chi = random.uniform(-180,180)
				self.pose.set_chi(k+1, len(self.pose.sequence()), chi)
			
		else:
			#Histogram based solutions
			i = 0
			aa_angles = []
			angles = []
			for aminoacid in primary_amino_sequence:
				if (i!=0) and (i!=len(primary_amino_sequence)-1):
					try:
						AA_Ant = sigla[primary_amino_sequence[i-1]]
						AA = sigla[primary_amino_sequence[i]]
						AA_Prox = sigla[primary_amino_sequence[i+1]]
						SS_Ant = siglaSS[str(secondary_sequence_list[i-1])]
						SS = siglaSS[str(secondary_sequence_list[i])]
						SS_Prox = siglaSS[str(secondary_sequence_list[i+1])]
							
					except:
						print("ERROR")
						traceback.print_exc()
						sys.exit(ERROR_CODE_fileopening+" Error 3")
						
					if USE_ANGLE_RANGE: #false
						#Use ranges of max and min angles
						aa_angles = [round(random.uniform(maxmin_angles[i][0],maxmin_angles[i][1]),floatprecision),#PHI
									 round(random.uniform(maxmin_angles[i][2],maxmin_angles[i][3]),floatprecision),#PSI
									 round(random.uniform(maxmin_angles[i][4],maxmin_angles[i][5]),floatprecision),#CHI1
									 round(random.uniform(maxmin_angles[i][6],maxmin_angles[i][7]),floatprecision),#CHI2
									 round(random.uniform(maxmin_angles[i][8],maxmin_angles[i][9]),floatprecision),#CHI3
									 round(random.uniform(maxmin_angles[i][10],maxmin_angles[i][11]),floatprecision)]#CHI4
					else:
						#Use histogram
						proba = []
						proba2 = []
						name = ""
						try:
							proba = hist.prob_hist[AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox]
							proba2 = hist.prob_hist[aminoacid][secondary_sequence_list[i]]
							name = AA_Ant+SS_Ant+AA+SS+AA_Prox+SS_Prox
						
						except:
							try:
								# #50% for each combination
								if (random.randint(1,10) <=5):
									proba = hist.prob_hist[AA_Ant+SS_Ant+AA+SS]
									proba2 = hist.prob_hist[aminoacid][secondary_sequence_list[i]]
									name = AA_Ant+SS_Ant+AA+SS
								else:
									proba = hist.prob_hist[AA+SS+AA_Prox+SS_Prox]
									proba2 = hist.prob_hist[aminoacid][secondary_sequence_list[i]]
									name = AA+SS+AA_Prox+SS_Prox
							except:
								proba = hist.prob_hist[aminoacid][secondary_sequence_list[i]]
								name = AA+SS
						
						aa_angles = hist.use_histogram(maxmin_angles[i], proba, PROB_RADIUS, proba2, name)
				else:
					proba = hist.prob_hist[aminoacid][secondary_sequence_list[i]]
					aa_angles = hist.use_histogram(maxmin_angles[i], proba, PROB_RADIUS, [], str(sigla[primary_amino_sequence[i]])+str(siglaSS[str(secondary_sequence_list[i])]))
				
				angles.append(copy.deepcopy(aa_angles))
				i+=1
			
			self.pose.set_phi(1, 0) #indiferente
			self.pose.set_psi(1,angles[0][1])
			self.pose.set_omega(1, angles[0][2])
			for kk in range(1,self.pose.residue(1).nchi()+1):
				try:
					self.pose.set_chi(kk, 1 ,angles[0][2+kk])
				except:
					break;
					
			i = 1
			for aminoacid in primary_amino_sequence:
				try:
					if (i+1) == len(self.pose.sequence()):
						#last amino		
						self.pose.set_omega(i+1, angles[i][2])
						self.pose.set_phi(i+1,angles[i][0])
						self.pose.set_psi(i+1, 0) #indiferente
						for kk in range(1,self.pose.residue(i+1).nchi()+1):
							try:
								self.pose.set_chi(kk, i+1,angles[i][2+kk])
							except:
								break;
						break;
						
					self.pose.set_phi(i+1,angles[i][0])
					self.pose.set_psi(i+1,angles[i][1])
					self.pose.set_omega(i+1, angles[i][2])
					
					for kk in range(1,self.pose.residue(i+1).nchi()+1):
						try:
							self.pose.set_chi(kk, i+1,angles[i][2+kk])
						except:
							break;
							
					i=i+1
				except:
					print("ERROR")
					traceback.print_exc()
					print("Line "+"- Error while creating geometry structure - Index "+str(i))


class WorkerThread(Thread):
	def __init__(self, id):
		Thread.__init__(self)
		self.id = id

	def select_rand_solution(self, solutions):
		index = 0
		counter_sol = 0
		for sol in solutions:
			if sol != None:
				counter_sol+=1
		counter_sol-=1 
		index = random.randint(0, counter_sol)
		return index
		

	def run(self):
		global agents
		global workers
		global hist_obj

		jump_radius_aux = TEST_JUMP_DIST
		agents[self.id].current.init_solution(hist_obj)
		agents[self.id].update()

		print 'WorkerThread %d: \n%s' % (self.id, agents[self.id])

		start_thread_time = datetime.datetime.now()
		agents[self.id].generation = 1

		best_energy = agents[self.id].pockets[0].energy_value
		gens_without_improve = 0
		gens_convergence = TEST_NOIMPROVE
		gens_start = 0
		start_log_time = start_thread_time

		while(datetime.datetime.now()-start_thread_time < TIME_LIMIT):

			#crossover it isn't allowed to execute on agent 0
			if agents[self.id].id_leader != None:
				if agents[self.id].pocket_leader != None:
					index_pocket_self_agent = self.select_rand_solution(agents[self.id].pockets)
					agents[self.id].crossover(
						agents[self.id].pocket_leader,
						agents[self.id].pockets[index_pocket_self_agent],
						CROSSOVER_PROB
						)
			else:
				index_pocket_self_agent = self.select_rand_solution(agents[self.id].pockets)
				agents[self.id].current = copy.deepcopy(agents[self.id].pockets[index_pocket_self_agent])

			#local search
			time_ls_start = datetime.datetime.now()
			agents[self.id].local_search(TEST_LS_PROB, TEST_JUMP_PROB, jump_radius_aux, hist_obj)
			agents[self.id].time_ls += datetime.datetime.now() - time_ls_start
			jump_radius_aux = jump_radius_aux * TEST_JUMP_FACT

			#update broadcast
			if agents[self.id].update():
				if agents[self.id].id_supporters:
					for id_supporter in agents[self.id].id_supporters:
						index_pocket_self_agent = self.select_rand_solution(agents[self.id].pockets)
						time_send_start = datetime.datetime.now()
						agents[id_supporter].queue_leader.put(copy.deepcopy(agents[self.id].pockets[index_pocket_self_agent]))
						agents[self.id].trx_send += 1
						agents[self.id].time_send += datetime.datetime.now() - time_send_start

				if agents[self.id].id_leader != None:
					time_send_start = datetime.datetime.now()
					agents[agents[self.id].id_leader].queue_pocket.put(copy.deepcopy(agents[self.id].pockets[0]))
					agents[self.id].trx_send += 1
					agents[self.id].time_send += datetime.datetime.now() - time_send_start

			#update pockets with supporter data 
			if agents[self.id].id_supporters:
				if not agents[self.id].queue_pocket.empty():
					time_receive_start = datetime.datetime.now()
					if agents[self.id].update(agents[self.id].queue_pocket.get()):
						pass
					agents[self.id].trx_receive += 1
					agents[self.id].time_receive += datetime.datetime.now() - time_receive_start
					agents[self.id].queue_pocket.task_done()

			#update pocket_leader with leader data 
			if agents[self.id].id_leader != None:
				if not agents[self.id].queue_leader.empty():
					time_receive_start = datetime.datetime.now()
					agents[self.id].pocket_leader = agents[self.id].queue_leader.get()
					agents[self.id].time_receive += datetime.datetime.now() - time_receive_start
					agents[self.id].trx_receive += 1
					agents[self.id].queue_leader.task_done()
			print '\033[33mGen %3d - WorkerThread %2d \033[36m- %s\033[39m' % (agents[self.id].generation, self.id, agents[self.id])
			agents[self.id].generation += 1

			if IF_RESET:
				if agents[self.id].id_leader == None:
					# reset control
					if best_energy == agents[self.id].pockets[0].energy_value:
						gens_without_improve += 1
					else:
						best_energy = agents[self.id].pockets[0].energy_value
						gens_without_improve = 0

					if gens_without_improve == gens_convergence:
						for agent in agents:
							agent.event_restart.set()

				#is event restart set?
				if agents[self.id].event_restart.is_set():
					if agents[self.id].id_leader == None:
						#only the root leader can keep the best solution
						agents[self.id].pockets = [agents[self.id].pockets[0]]+[None for i in range(1,NUM_POCKETS)]
					else:
						agents[self.id].pockets = [None for i in range(0,NUM_POCKETS)]
						agents[self.id].pocket_leader = None
					agents[self.id].queue_leader = PriorityQueue(1)
					agents[self.id].queue_pocket = PriorityQueue()
					agents[self.id].restarts += 1

					print '\033[35mRESTARTING %3d - WorkerThread %2d - %s\033[39m' % (agents[self.id].restarts, self.id, agents[self.id])
					agents[self.id].current.init_solution(hist_obj)
					agents[self.id].update()
					jump_radius_aux = TEST_JUMP_DIST
					gens_convergence = agents[self.id].generation - gens_start
					gens_start = agents[self.id].generation
					gens_without_improve = 0
					agents[self.id].event_restart.clear()
					print '\033[95mRESTARTED %3d - WorkerThread %2d - %s\033[39m' % (agents[self.id].restarts, self.id, agents[self.id])

			if (datetime.datetime.now()-start_log_time) > TIME_WINDOW_STATUS:
				agents[self.id].status_log_append( datetime.datetime.now()-start_thread_time )
				start_log_time = datetime.datetime.now()

		#update pockets with supporter data 
		if agents[self.id].id_supporters != None:
			supporters_dead = False
			while not supporters_dead:
				for id_supporter in agents[self.id].id_supporters:
					if workers[id_supporter].is_alive():
						supporters_dead = False
						break
					else:
						supporters_dead = True

				if not agents[self.id].queue_pocket.empty():
					if agents[self.id].update(agents[self.id].queue_pocket.get()):
						pass
					agents[self.id].queue_pocket.task_done()

		if agents[self.id].id_leader != None:
			agents[agents[self.id].id_leader].queue_pocket.put(copy.deepcopy(agents[self.id].pockets[0]))

		final_thread_time=datetime.datetime.now()-start_thread_time
		agents[self.id].status_log_append(final_thread_time)
		print '\n\033[34m************ WorkerThread %d done ************\033[39m\n' % (self.id)


def main():
	global agents
	global input_file
	global primary_amino_sequence
	global secondary_amino_sequence
	global hist_obj
	global scorefxn
	global workers
	global receivers
	global NUM_AGENTS
	global NUM_POCKETS
	global NUM_SUP
	global IF_RESET

	# sys.stdout = os.devnull
	# sys.stderr = os.devnull

	rosetta.init()
	scorefxn = create_score_function('talaris2013')
	agents = []
	protein=str(sys.argv[1]).strip().upper()
	NUM_AGENTS = int(sys.argv[2])
	NUM_POCKETS = int(sys.argv[3])
	NUM_SUP = int(sys.argv[4])
	IF_RESET = (str(sys.argv[5]) == 'True')

	input_file_name  = "TXT/%s.txt" % (protein)
	directory_results = "results_parallel/%s/%s" % (protein, 'log-'+str(uuid.uuid4()))
	
	try:
		input_file = open(input_file_name)
	except:
		print("ERROR")
		traceback.print_exc()
		sys.exit("Line "+str(lineno())+": "+ERROR_CODE_fileopening+input_file_name)
	
	#Ignores the comment lines in the beginning of the file
	line = input_file.readline()
	while(line.startswith('#')):
		line = input_file.readline()
		if not line:
			sys.exit(ERROR_CODE_emptyfile)
	amino_sequence = line.strip()
	line = input_file.readline()
	
	secondary_sequence = line.strip()
	
	primary_amino_sequence=amino_sequence
	secondary_amino_sequence=secondary_sequence

	print '\nPrimary Amino Sequence: %s' % (primary_amino_sequence)
	print 'Secondary Amino Sequence: %s\n' % (secondary_amino_sequence)

	for ss in secondary_amino_sequence:
		if   ss == "B":
			ssi = SS_B
		elif ss == "C":
			ssi = SS_C
		elif ss == "E":
			ssi = SS_E
		elif ss == "G":
			ssi = SS_G
		elif ss == "H":
			ssi = SS_H
		elif ss == "I":
			ssi = SS_I
		elif ss == "T":
			ssi = SS_T
		else:
			sys.exit("Line "+": "+"Secondary structure unknown: "+ss+" "+str(ssi))
		secondary_sequence_list.append(ssi)

	
	hist_obj=HistogramFiles()
	hist_obj.read_histograms()
	print

	for i in range(0,NUM_AGENTS):
		agents.append(Agent(i))
		workers.append(WorkerThread(i))
		workers[i].setName('Worker for Agent '+str(i))

	start_time = datetime.datetime.now()
	for worker in workers:
		worker.start()

	for worker in workers:
		worker.join()

	t_final=datetime.datetime.now()-start_time

	if not os.path.exists(directory_results):
		os.makedirs(directory_results)
	fout = open('%s/run-summary.txt'%(directory_results),'w')

	print 'Parameters'
	fout.write('Parameters\n')
	print '--- pockets: %d' % (NUM_POCKETS)
	fout.write('--- pockets: %d\n' % (NUM_POCKETS))
	print '--- agents: %d' % (NUM_AGENTS)
	fout.write('--- agents: %d\n' % (NUM_AGENTS))
	print '--- supporters per leader: %d' % (NUM_SUP)
	fout.write('--- supporters per leader: %d\n' % (NUM_SUP))
	print '--- do reset: %s' % (str(IF_RESET))
	fout.write('--- do reset: %s\n' % (str(IF_RESET)))
	print '--- max generations: %d' % (MAX_GEN)
	fout.write('--- max generations: %d\n' % (MAX_GEN))
	print '--- prob radius: %f' % (PROB_RADIUS)
	fout.write('--- prob radius: %f\n' % (PROB_RADIUS))
	print '--- diversity: %f' % (TEST_DIV_DIFF)
	fout.write('--- diversity: %f\n' % (TEST_DIV_DIFF))
	print '--- prob of ls: %f' % (TEST_LS_PROB)
	fout.write('--- prob of ls: %f\n' % (TEST_LS_PROB))
	print '--- simulated annealing decrease factor: %f' % (TEST_LS_FACT)
	fout.write('--- simulated annealing decrease factor: %f\n' % (TEST_LS_FACT))
	print '--- prob of jump before ls: %f' % (TEST_JUMP_PROB)
	fout.write('--- prob of jump before ls: %f\n' % (TEST_JUMP_PROB))
	print '--- jump decrease factor: %f' % (TEST_JUMP_FACT)
	fout.write('--- jump decrease factor: %f\n' % (TEST_JUMP_FACT))
	print '--- initial temperature for simulated annealing: %d' % (TEST_TEMP_INIT)
	fout.write('--- initial temperature for simulated annealing: %d\n' % (TEST_TEMP_INIT))
	print '--- initial max jump distance: %f' % (TEST_JUMP_DIST)
	fout.write('--- initial max jump distance: %f\n' % (TEST_JUMP_DIST))
	print '--- generations without improvements: %d' % (TEST_NOIMPROVE)
	fout.write('--- generations without improvements: %d\n' % (TEST_NOIMPROVE))
	print '--- prob of crossover: %f' % (CROSSOVER_PROB)
	fout.write('--- prob of crossover: %f\n' % (CROSSOVER_PROB))
	print '--- max time of execution: %s' % (str(TIME_LIMIT))
	fout.write('--- max time of execution: %s\n\n' % (str(TIME_LIMIT)))

	for i in range(0,NUM_AGENTS):
		if agents[i].pockets[0] != None:
			agents[i].pockets[0].pose.dump_pdb('%s/agent-%02d.pdb'%(directory_results,i))
			print '\n%s' % (agents[i])
			fout.write('%s\n' % (agents[i]))		

			print 'Total generation of agent_%02d: %d' % (i,agents[i].generation)
			fout.write('Total generation of agent_%02d: %d\n' % (i,agents[i].generation))

			print 'Total restarts of agent_%02d: %d' % (i,agents[i].restarts)
			fout.write('Total restarts of agent_%02d: %d\n' % (i,agents[i].restarts))

			print 'Total time LocalSearch of agent_%02d: %s' % (i,str(agents[i].time_ls))
			fout.write( 'Total time LocalSearch of agent_%02d: %s\n' % (i,str(agents[i].time_ls)))

			print 'Total time SEND of agent_%02d: %s' % (i,str(agents[i].time_send))
			fout.write( 'Total time SEND of agent_%02d: %s\n' % (i,str(agents[i].time_send)))

			print 'Total transactions SEND of agent_%02d: %s' % (i,str(agents[i].trx_send))
			fout.write( 'Total transactions SEND of agent_%02d: %s\n' % (i,str(agents[i].trx_send)))

			print 'Total time RECEIVE of agent_%02d: %s' % (i,str(agents[i].time_receive))
			fout.write( 'Total time RECEIVE of agent_%02d: %s\n' % (i,str(agents[i].time_receive)))

			print 'Total transactions RECEIVE of agent_%02d: %s\n' % (i,str(agents[i].trx_receive))
			fout.write( 'Total transactions RECEIVE of agent_%02d: %s\n\n' % (i,str(agents[i].trx_receive)))

			agents[i].status_write('%s/log-agent-%02d.txt'%(directory_results,i))

	fout.write('\nTotal time: %s\n' % (str(t_final)))
	fout.close()
	print 'Total time: '+str(t_final)
	print '**** end MainThread ****'	

if __name__ == '__main__':
	main()
