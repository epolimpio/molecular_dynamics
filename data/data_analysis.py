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
EXCLUSION = [1, 2, 3, 4, 9, 10, 11, 17, 25, 38]

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
	stds = []
	for box in data:
		means.append(np.mean(box))
		stds.append(np.std(box))

	return [np.mean(np.array(means)), np.mean(np.array(stds))/sqrt(num_boxes)]

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
	all_data = []
	var_data = []
	for box in data:
		mean = np.mean(box)
		ratio = np.power(np.divide(np.subtract(box, mean),mean),2)
		all_data.append(np.mean(ratio))
		var_data.append(np.std(ratio))

	cv = 1/(2.0/3.0 - N*np.mean(np.array(all_data))) - 1.5
	cv_err = abs(pow(cv,2)*N*np.mean(np.array(var_data))/sqrt(num_boxes))

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
# calculates the physical quantities
# we also write a table in latex format to insert into the report
densities = []
temps = []
temp_set = []
temps_err = []
cvs = []
cvs_err = []
press = []
press_err = []
en = []
en_err = []

file_run = open("resume.txt", "w")
for run in range(0, NUM_FILES):

	if (run+1 not in EXCLUSION):

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

		file_run.write(fileref+" T->"+str(round(temp,2))+ 
				" Dens->"+str(round(density,2))+"\n")

		# Calculate physical variables and its error
		[t_out, t_err] = mean_error(temp_vec)
		[virial_out, virial_err] = mean_error(virial_vec)
		[Epp, Epp_err] = mean_error(np.divide(Etot, N))
		[cv_out, cv_err] = calc_heat_cap(Ek, N)
		[p_out, p_err] = calc_pressure(virial_out, virial_err, t_out, 
							t_err, cutoff, N, density)

		# Append the data to the vectors for plotting
		densities.append(density)
		temps.append(t_out)
		temp_set.append(temp)
		temps_err.append(t_err)
		cvs.append(cv_out)
		cvs_err.append(cv_err)
		press.append(p_out)
		press_err.append(p_err)
		en.append(Epp)
		en_err.append(Epp_err)
file_run.close()

# Organize the data by temperature and density
np_dens = np.array(densities)
graph = {}
unique_dens = np.unique(np_dens).tolist()
for density in unique_dens:
	
	graph[density] = {}
	indexes = np.where(np_dens == density)[0].tolist()

	t = list(np.array(temps)[indexes])
	t_set = list(np.array(temp_set)[indexes])
	t_err = list(np.array(temps_err)[indexes])
	p = list(np.array(press)[indexes])
	p_err = list(np.array(press_err)[indexes])
	cv = list(np.array(cvs)[indexes])
	cv_err = list(np.array(cvs_err)[indexes])
	e = list(np.array(en)[indexes])
	e_err = list(np.array(en_err)[indexes])
	
	i = 0
	for temperature in t_set:
		graph[density][temperature] = {}
		
		# Pressure Data
		graph[density][temperature]["T"] = t[i]
		graph[density][temperature]["T_ERR"] = t_err[i]
		graph[density][temperature]["P"] = p[i]
		graph[density][temperature]["P_ERR"] = p_err[i]		
		graph[density][temperature]["CV"] = cv[i]
		graph[density][temperature]["CV_ERR"] = cv_err[i]
		graph[density][temperature]["E"] = e[i]
		graph[density][temperature]["E_ERR"] = e_err[i]
		i += 1	

# Write Latex table
latex = open("latex_table.txt", "w")
for density in graph.keys():
	for temperature in graph[density].keys():

		t_out = graph[density][temperature]["T"]
		t_err = graph[density][temperature]["T_ERR"]
		p_out = graph[density][temperature]["P"]
		p_err = graph[density][temperature]["P_ERR"]
		cv_out = graph[density][temperature]["CV"]
		cv_err = graph[density][temperature]["CV_ERR"]
		Epp = graph[density][temperature]["E"]
		Epp_err = graph[density][temperature]["E_ERR"]

		# Format the numbers to the correct digit for table
		t = significant_digit(t_out, t_err)
		p = significant_digit(p_out, p_err)
		cv = significant_digit(cv_out, cv_err)
		E = significant_digit(Epp, Epp_err)

		# Write in latex table
		latex_row = str(round(density,2)) + " & "
		latex_row += str(round(temperature,2)) + " & "
		latex_row += str(t[0]) + "(" + str(t[1]) + ") & "
		latex_row += str(E[0]) + "(" + str(E[1]) + ") & "
		latex_row += str(p[0]) + "(" + str(p[1]) + ") & "
		latex_row += str(cv[0]) + "(" + str(cv[1]) + ") \\\\ \n"
		latex.write(latex_row)

latex.close()

# Plot Pressure
plt.figure(1)
legend = []
for density in sorted(graph.keys()):
	if round(density,2) not in [0.30, 1.20]:
		legend.append(str(round(density,2)))
		x = []
		x_err = []
		y = []
		y_err = []
		for temperature in sorted(graph[density].keys()):
			x.append(graph[density][temperature]["T"])
			x_err.append(graph[density][temperature]["T_ERR"])
			y.append(graph[density][temperature]["P"])
			y_err.append(graph[density][temperature]["P_ERR"])

		plt.errorbar(x, y, xerr = x_err, yerr=y_err, fmt='-o')

plt.xlabel('T')
plt.ylabel(r'$\beta P / \rho$')
plt.legend(legend, loc='lower right', fontsize=12, title="Density")

# Plot Heat Capacity
plt.figure(2)
legend = []
for density in sorted(graph.keys()):
	
	legend.append(str(round(density,2)))
	x = []
	x_err = []
	y = []
	y_err = []
	for temperature in sorted(graph[density].keys()):
		x.append(graph[density][temperature]["T"])
		x_err.append(graph[density][temperature]["T_ERR"])
		y.append(graph[density][temperature]["CV"])
		y_err.append(graph[density][temperature]["CV_ERR"])

	plt.errorbar(x, y, xerr = x_err, yerr=y_err, fmt='-o')

plt.xlabel('T')
plt.ylabel(r'$C_{pot}/N$')
plt.legend(legend, loc='lower left', fontsize=12, title="Density")

plt.show()
