#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 01:31:34 2019

@author: Ayca
"""
import numpy as np
#%% File Reading 

def fileReadProfile(filepath):
    text = open(filepath, "r") 
    text = text.read()
    alphabet = 'AGCT-\n'
    profile = ''.join([i for i in text if i in alphabet])
    profile = profile.split("\n")
    return profile

def fileReadSequence(filepath):
    seq = ""
    with open(filepath) as fp:
        line = fp.readline()
        while line:
            line = fp.readline()
            print(line)
            line = line.replace("\n","")
            seq += line 
    return seq

#profile = fileReadProfile("aligned_sequences.aln")
#sequence = fileReadSequence("seq.fasta")
#%%
def Matrix(row,col,gap):
    matrix = np.zeros((row+1,col+1),dtype= float)
    for i in range(1,row+1):
        matrix[i][0] = matrix[i-1][0]+gap
    for i in range(1,col+1):
        matrix[0][i] = matrix[0][i-1]+gap
    return matrix
#%%
def Score(p, char, j, match, mismatch):
    score = 0
    if char == 'A':
        score = match * p[0][j]
        for i in range(1,5):
            score += mismatch*p[i][j]
    elif char == 'C':
        ######################
#        print("matched with C")
#        print("mismatch * p[0][j]",mismatch * p[0][j])
#        print("match*p[1][j]",match*p[1][j])
        ######################
        score = mismatch * p[0][j] + match*p[1][j]
        for i in range(2,5):
            score += mismatch*p[i][j]
    elif char == 'G':
        ######################
#        print("matched with G")
        ######################
        for i in range(0,2):
            score += mismatch*p[i][j]
        score += match*p[2][j]
        for i in range(3,5):
            score += mismatch*p[i][j]
    elif char == 'T':
        ######################
#        print("matched with T")
        ######################
        for i in range(0,3):
            score += mismatch*p[i][j]
        score += match*p[3][j]+mismatch*p[4][j]
    elif char == '-':
        ######################
#        print("matched with -")
        ######################
        for i in range(0,4):
            score += mismatch*p[i][j]
        score += match*p[4][j]

    
#    print("Score", score)

    return score
#%%
def Align(matrix,p,s,match,mismatch,gap,i,j):
    align1 = ""
    while i >= 1 and j >= 1:
        current_score = matrix[i][j]
        if current_score == matrix[i-1][j-1] + Score(p,s[i-1],j-1,match,mismatch):
            align1 = align1+ s[i-1]
            i -=1
            j -=1
        elif current_score == matrix[i][j-1] + gap:
            #align1 = align1 + s1[i]
            align1 = align1 + "-"
            j -=1
        elif current_score == matrix[i-1][j] + gap:
            #align1 = align1 + "-"
            align1 = align1 + "-"
            i -=1
######################################################
#        if current_score == matrix[i][j-1] + gap:
#            #align1 = align1 + s1[i]
#            align1 = align1 + "-"
#            j -=1
#        elif current_score == matrix[i-1][j] + gap:
#            #align1 = align1 + "-"
#            align1 = align1 + "-"
#            i -=1
#        elif current_score == matrix[i-1][j-1] + Score(p,s[i-1],j-1,match,mismatch):
#            align1 = align1+ s[i-1]
#            i -=1
#            j -=1
#        

######################################################        
#        if current_score == matrix[i-1][j] + gap:
#            #align1 = align1 + "-"
#            align1 = align1 + "-"
#            i -=1
#        
#        elif current_score == matrix[i][j-1] + gap:
#            #align1 = align1 + s1[i]
#            align1 = align1 + "-"
#            j -=1
#        elif current_score == matrix[i-1][j-1] + Score(p,s[i-1],j-1,match,mismatch):
#            align1 = align1+ s[i-1]
#            i -=1
#            j -=1
    align1 = align1[::-1]
    #align2 = align2[::-1]
    #print("align1: {}, align2: {}".format(align1,align2))
    return align1

#%%
def Global(p, s,match,mismatch,gap):
    matrix = Matrix(len(s), len(p[0]),gap)
    row = len(s)+1
    col = len(p[0])+1
    for i in range(1,row):
        for j in range(1,col):
            diagonal = 0
            gap_left = 0
            gap_up = 0
            diagonal = matrix[i-1][j-1]+ Score(p,s[i-1],j-1,match,mismatch)
            gap_left = matrix[i][j-1] +gap
            gap_up = matrix[i-1][j]+gap
            matrix[i][j] = max(diagonal,gap_left,gap_up)
#            if i == 5 and j == 4:
#                print("diagonal",diagonal)
#                print("matrix[i-1][j-1]",matrix[i-1][j-1])
#                print("Score(p,s[i-1],j-1,match,mismatch)",Score(p,s[i-1],j-1,match,mismatch))
#                print("s[i-1]",s[i-1])
#                print("gap_left={}-2={}".format(matrix[i][j-1],gap_left))
#                print("gap_up ={}-2={}".format(matrix[i-1][j],gap_up))
#                print("matrix[i][j]",matrix[i][j])

    score = matrix[row-1,col-1]
    print(score)
    i = row-1
    j = col-1
    align = Align(matrix,p,s,match,mismatch,gap,i,j)
    print(align)
    return matrix, align,score
#%%
#profile = fileReadProfile("aligned_sequences.aln")
#sequence = fileReadSequence("seq.fasta")
#profile_scoring = np.zeros((5, len(profile[0])),dtype = float)
##ACGT-
#for i in profile:
#    for j in range(len(i)):
#        if i[j] == 'A':
#            profile_scoring[0][j] += 0.20
#        elif i[j] == 'C':
#            profile_scoring[1][j] += 0.20
#        elif i[j] == 'G':
#            profile_scoring[2][j] += 0.20
#        elif i[j] == 'T':
#            profile_scoring[3][j] += 0.20
#        elif i[j] == '-':
#            profile_scoring[4][j] += 0.20
#result = Global(profile_scoring,sequence,1,-1,-2)
#mmm = result[0]
#aaa = result[1]
#sss = result[2]

#%%
##careful with that axe, eugene

from shutil import copyfile
def Main():
    out_file = open("seq.aln","w+")
    print("Type the name of aln file which includes multiple aligned patterns:")
    file = str(input())
    copyfile(file,"seq.aln")
    profile = fileReadProfile(file)
    
    print("Type the name of sequence file:")
    file = str(input())
    sequence = fileReadSequence(file)
    print("Type match score")
    match = int(input())
    print("Type mismatch score")
    mismatch = int(input())
    print("Type gap score")
    gap = int(input())
    profile_scoring = np.zeros((5, len(profile[0])),dtype = float)
    #ACGT-
    for i in profile:
        for j in range(len(i)):
            if i[j] == 'A':
                profile_scoring[0][j] += 0.20
            elif i[j] == 'C':
                profile_scoring[1][j] += 0.20
            elif i[j] == 'G':
                profile_scoring[2][j] += 0.20
            elif i[j] == 'T':
                profile_scoring[3][j] += 0.20
            elif i[j] == '-':
                profile_scoring[4][j] += 0.20
    result = Global(profile_scoring,sequence,match,mismatch,gap)
    mmm = result[0]
    aaa = result[1]
    sss = result[2]

    out_str = "\nsequence  "+aaa
    with open('seq.aln', 'a') as writer:
        writer.write(out_str)

Main()