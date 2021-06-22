print("Composed by Dr Wen-Xuan Wang, CSU. version 3.1.2021.06.22. Xiangya School of Pharmaceutical Science.")
print("Please cite J. Org. Chem. 2020, 85, 17, 11350–11358")
print("作者：王文宣，中南大学，湘雅药学院")
print(" ")

#Parameters table
temp = 298 #termperature for Boltzmann population calculation (K)
threshold_geom = 0.1 #threshold to determine geometery duplicates (Angstrom)
threshold_energy = 0.2 #Gibbs free energy threshold to determine geometery duplicates (kcal/mol)

print("Termperature for Boltzmann population calculation (K) is ", temp)
print("threshold to determine geometery duplicates (Angstrom) is ", threshold_geom)
print("Gibbs free energy threshold to determine geometery duplicates (kcal/mol) is", threshold_energy)
print(" ")

import os
import glob
import numpy as np
import linecache
import pandas as pd
import csv
import math

filenames = []
gibbs_list = []
duplicates_list = []

class Conf_info():
  filename = ''
  gibbs = 0.0
  imaginary = 0
  error = ' '
  def __init__(self,filename):
    self.filename = filename

outfilenames = glob.glob('*.out') #obtain all out file names
logfilenames = glob.glob('*.log') #obtain all log file names
filenames.extend(outfilenames)
filenames.extend(logfilenames) #get the file list

#open files and read the information sequently
if len(filenames) == 0:
    print ("No out or log files found")
    quit()
else:
    conf_list = []
    for i in range(len(filenames)):
      atom_num = [] #for reading shielding tensors
      shielding = [] #for reading shielding tensors
      conf = Conf_info(filenames[i])
      with open(filenames[i], 'r') as f:
        gibbs_line_num = -1
        for line_number, key_word in enumerate(f.readlines()):
              if ' imaginary frequencies (negative Signs) ******' in key_word: 
                #read imaginary frequencies
                conf.imaginary = 1                             
              elif 'Sum of electronic and thermal Free Energies=' in key_word: 
                gibbs_line_num = line_number
                conf.gibbs = float(key_word.split( )[-1]) #build the Gibbs free energy list
              elif '-- Number of steps exceeded' in key_word:
                conf.error = 'Optimization failed because max steps exceeded'
              elif 'C    Isotropic =' in key_word:  #read the shielding tensor data of C                
                shielding_line_num = line_number
                shielding_line = linecache.getline(filenames[i], shielding_line_num+1)
                atom_num.append(shielding_line.split( )[0])
                shielding.append(shielding_line.split( )[4])
        if len(atom_num) > 0:     
          prefix = filenames[i][:-4]     
          shielding_file = open(prefix+'_C_shielding.txt', 'w')
          for j in range(len(atom_num)):
            print(atom_num[j], "  ", shielding[j], file = shielding_file)  #save shielding data into text file
          shielding_file.close()

      conf_list.append(conf)
      del conf
           
filenames_selected = []          
#Remove structures with imaginary frequencies or without gibbs free energy
for i in range(len(conf_list)):
  if conf_list[i].imaginary == 1:
    print(conf_list[i].filename, " has imaginary frequencies, and will not be considered")
    continue
  elif conf_list[i].gibbs == 0:
    print(conf_list[i].filename, " has no Gibbs free energy information", conf_list[i].error)
    continue
  else:
    filenames_selected.append(conf_list[i].filename)
    gibbs_list.append(conf_list[i].gibbs)


filenames_selected = np.array(filenames_selected)
gibbs_list = np.array(gibbs_list,dtype=float)


#sort the information matrix by Gibbs free energy
gibbs_order = np.argsort(gibbs_list)

filenames_sort = filenames_selected[gibbs_order]
gibbs_float_sort = gibbs_list[gibbs_order]

def dis_matrix_cp(file_nameA, file_nameB, threshold):
  #This function uses the coordinates for multiple steps calculation to identify duplicates
  #read the line number of coordinates
  with open(file_nameA, 'r') as A:
    line0 = -1
    line1 = -1
    for line_number, key_word in enumerate(A.readlines()):
      if 'Redundant internal coordinates found in file.  (old form).' in key_word:
            line0 = line_number
      if 'Recover connectivity data from disk.' in key_word:
            line1 = line_number
    if line0 == -1:
      print("There may be no coordinates for multiple steps calculation in ", file_nameA)
      return 0

  with open(file_nameB, 'r') as B:
    line2 = -1
    line3 = -1
    for line_number, key_word in enumerate(B.readlines()):
      if 'Redundant internal coordinates found in file.  (old form).' in key_word:
            line2 = line_number
      if 'Recover connectivity data from disk.' in key_word:
            line3 = line_number
    if line2 == -1:
      print("There may be no coordinates for multiple steps calculation in ", file_nameB) 
      return 0    

    #read the lines of coordinates and compare the distance between atoms
  k = 1 #record the linenumber of the coordinate line for processing
  distanceA = []
  distanceB = []
  for i in range(line0+2, line1-1):
      A_coordinate_line1 = linecache.getline(file_nameA, i)
      arrayA_atom1 = np.array(A_coordinate_line1.split(',')[2: ])
      arrayA_atom1_float = arrayA_atom1.astype(float)  

      k = k+1
      
      B_coordinate_line1 = linecache.getline(file_nameB, line2+k)
      arrayB_atom1 = np.array(B_coordinate_line1.split(',')[2: ])
      arrayB_atom1_float = arrayB_atom1.astype(float)

      
      
      l = 0 #record the linenumber of the coordinate line for processing
      for j in range(i, line1):
        A_coordinate_line2 = linecache.getline(file_nameA, j+1)
        arrayA_atom2 = np.array(A_coordinate_line2.split(',')[2: ])
        arrayA_atom2_float = arrayA_atom2.astype(float)   

        l = l + 1
        B_coordinate_line2 = linecache.getline(file_nameB, line2+k+l)
        arrayB_atom2 = np.array(B_coordinate_line2.split(',')[2: ])
        arrayB_atom2_float = arrayB_atom2.astype(float)

        distanceA.append(((arrayA_atom1_float[0]-arrayA_atom2_float[0]) ** 2 + (arrayA_atom1_float[1]-arrayA_atom2_float[1]) ** 2 + (arrayA_atom1_float[2]-arrayA_atom2_float[2]) ** 2) **0.5)
        distanceB.append(((arrayB_atom1_float[0]-arrayB_atom2_float[0]) ** 2 + (arrayB_atom1_float[1]-arrayB_atom2_float[1]) ** 2 + (arrayB_atom1_float[2]-arrayB_atom2_float[2]) ** 2) **0.5)
        
  distanceA = np.array(distanceA,dtype='float32')
  distanceB = np.array(distanceB,dtype='float32')
  distanceA_sort = np.sort(distanceA)
  distanceB_sort = np.sort(distanceB)

  for i in range(len(distanceA_sort)):
    if abs(distanceA_sort[i] - distanceB_sort[i]) > threshold: #Check if the distance matrices of conformers are the same
      return 0
  return 1

#check duplicates by energy and geometry
for i in range(len(filenames_sort)-1):
  gibbs_difference = gibbs_float_sort[i+1] - gibbs_float_sort[i]
  if gibbs_difference < (threshold_energy / 627.5094):
    if dis_matrix_cp(filenames_sort[i+1],filenames_sort[i],threshold_geom) == 1:
       print(filenames_sort[i+1], "is the same as", filenames_sort[i])
       duplicates_list.append([i+1])

#remove duplicates
if len(duplicates_list) > 0:
  duplicates_list = np.array(duplicates_list)
  duplicates_list = duplicates_list.astype(int)
  for i in range(len(duplicates_list)):    
    filenames_sort = np.delete(filenames_sort, (duplicates_list[i]-i)) 
    gibbs_float_sort = np.delete(gibbs_float_sort, (duplicates_list[i]-i))
  
    

#save energy information to excel file
gibbs_kcal = []
for i in range(len(gibbs_float_sort)):
  gibbs_kcal.append((gibbs_float_sort[i]-gibbs_float_sort[0])*627.5094) #Hatree converted into kcal/mol
gibbs_kcal = np.array(gibbs_kcal)
gibbs_kcal = gibbs_kcal.astype(float)
boltzmann_factor = 2.7182818284590452353602874**(-1*gibbs_kcal/temp/0.0019858775) #calculate Boltzmann factors
sum_boltzmann = 0
for i in range(len(boltzmann_factor)):
  sum_boltzmann = boltzmann_factor[i]+sum_boltzmann #sum Boltzmann factors
population = boltzmann_factor/sum_boltzmann #calculate population of each conformer

writer = pd.ExcelWriter(r'population.xlsx')
df1 = pd.DataFrame({'Name':filenames_sort, 'Gibbs free energy (hatree)':gibbs_float_sort, 'ΔG (kcal/mol)':gibbs_kcal, 'Population':population})
df1.to_excel(writer, sheet_name='sheet1') 
writer.save() #write energy data into excel file

#average 13C shielding tensors
shielding_average = [0.0]*len(atom_num) 
for i in range(len(filenames_sort)):
  with open(filenames_sort[i][:-4]+'_C_shielding.txt') as f:    
    shielding_atom = []
    for j in range(len(atom_num)):      
      shielding_line = linecache.getline(filenames_sort[i][:-4]+'_C_shielding.txt', j+1).split( )[1] #read the file of shielding tensors
      shielding_atom.append(float(shielding_line)*population[i]) #weight the shielding tensors with population value
    
    for k in range(len(atom_num)):
        shielding_average[k] = shielding_average[k] + shielding_atom[k] #accumulate the weighted shielding tensors

writer = pd.ExcelWriter(r'Averaged_C_shielding.xlsx')
df1 = pd.DataFrame({'Atom number':atom_num, 'Averaged shielding tensors':shielding_average})
df1.to_excel(writer, sheet_name='sheet1')
writer.save() #save the averaged shielding tensors into file






            




    
