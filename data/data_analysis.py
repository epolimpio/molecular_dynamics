#!/usr/bin/python

# Declare libraries
import matplotlib as matl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab as P
from math import *
import csv

# Declare constants
PI = np.pi

LINES_HEADER = 5
NUM_FILES = 64
BOX_SIZE = 250

# Functions
def read_header(filename):

	f = open(filename)
	data = {}
	for i in range(0,LINES_HEADER):
		txt = f.readline().split("->")
		data[txt[0].strip()] = float(txt[1])
	f.close()
	return data

def mean_error(vector):

	num_boxes = len(vector)//BOX_SIZE
	data = np.split(vector, num_boxes)
	means = []
	for box in data:
		means.append(np.mean(box))

	boxes = np.array(means)

	return [np.mean(boxes), np.std(boxes)]

def significant_digit(number, error):

	power_10 = int(floor(log10(abs(error))))
	digit = round(error/pow(10,power_10),0)
	if digit == 10:
		digit = 1
		power_10 = power_10 + 1

	if power_10 > 0:
		digit = digit*pow(10,power_10)

	return [round(number, -power_10), int(digit)]

def calc_heat_cap(Ek, N):
# Calculate heat capacity from Kinetic Energy
	num_boxes = len(Ek)//BOX_SIZE
	data = np.split(Ek, num_boxes)
	ratio = []
	for box in data:
		ratio.append(np.var(box)/pow(np.mean(box),2))
	boxes = np.array(ratio)
	cv = 1/(1/(1.5*N) - np.mean(boxes))
	cv_err = abs(pow(cv,2)*np.std(boxes))

	return cv, cv_err

def calc_pressure(virial, virial_err, temp, 
					temp_err, cutoff, N, density):
# Calculate pressure from virial terms
	correction = 16*PI*density/(3*temp)*(2/(3*cutoff**9) - 1/(cutoff**3))
	correction_err = abs(correction*temp_err/temp)
	
	virial_term = virial/(3*N*temp)
	virial_term_err = abs(virial_term*sqrt((virial_err/virial)**2 + (temp_err/temp)**2))
	
	p = 1 + virial_term + correction
	p_err = virial_term_err + correction_err

	return [p, p_err]

# Main program
# calculates the physical quantities and write them in a file
# we also write a table in latex format to insert into the report
csvfile = open("consolidated2.dat", "w")
latex = open("latex_table.txt", "w")
fields = ["DENSITY", "T_SET", "T_OUT", "T_OUT_ERR", "CV", "CV_ERR", "P", "P_ERR"]
writer = csv.DictWriter(csvfile, fieldnames = fields, delimiter = "\t")
writer.writeheader()

for run in range(0, NUM_FILES):

	fileref = '%0*d' % (3, run+1)
	virial_vec = np.loadtxt("pressure"+fileref+".dat", skiprows=LINES_HEADER)
	energy = np.loadtxt("energy"+fileref+".dat", skiprows=LINES_HEADER)
	settings = read_header("pressure"+fileref+".dat")
	Ek = energy[:,0]
	Ep = energy[:,1]
	Etot = energy[:,2]
	temp_vec = energy[:,3]

	N = settings['N']
	density = settings['DENSITY']
	cutoff = settings['CUTOFF']
	temp = settings['TEMP']

	row = {}
	[t_out, t_err] = mean_error(temp_vec)
	[virial_out, virial_err] = mean_error(virial_vec)
	[cv_out, cv_err] = calc_heat_cap(Ek, N)
	[p_out, p_err] = calc_pressure(virial_out, virial_err, t_out, 
						t_err, cutoff, N, density)

	row["DENSITY"] = density
	row["T_SET"] = temp
	row["T_OUT"] = t_out
	row["T_OUT_ERR"] = t_err
	row["CV"] = cv_out
	row["CV_ERR"] = cv_err
	row["P"] = p_out
	row["P_ERR"] = p_err

	t = significant_digit(t_out, t_err)
	p = significant_digit(p_out, p_err)
	cv = significant_digit(cv_out, cv_err)

	latex_row = str(round(density,2)) + " & "
	latex_row += str(round(temp,2)) + " & "
	latex_row += str(t[0]) + "(" + str(t[1]) + ") & "
	latex_row += str(p[0]) + "(" + str(p[1]) + ") & "
	latex_row += str(cv[0]) + "(" + str(cv[1]) + ") \\\\ \n"

	latex.write(latex_row)

	writer.writerow(row)

csvfile.close()
latex.close()
		

