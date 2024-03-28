https://dl.dropboxusercontent.com/u/2838022/python_courses/eg_ap_2015/index.html

#Benchmarking (how long does something take)

#Unix time

#How long does our program take to run? On Linux/Mac we can do

#time somecommand

#In iPython, prefix shell commands with !

#!date

#Wed Jul 29 10:27:51 BST 2015

#Given output that looks like this:

#real    0m0.490s 
#user    0m0.457s 
#sys        0m0.032s

#real is the wallclock time (affected by busy systems and other programs)
#user is the time spent executing our code
#sys is the time spent waiting for system calls (file IO, memory, network)
#user+sys is probably the most useful.


#time.time() gives us current UNIX epoch (number of seconds since midnight January 1st 1970 (don't ask.)

import time 
start = time.time() 

# print the sum of the first million cube numbers
x = 0 
for i in range(1000000): 
    x = x + i ** 3 
print(x) 
 
end = time.time() 
print(end - start) 

#timing in python
#Python has a built in module for doing timing. From the command line:

!python -m timeit "4 ** 10"

# - automatically runs the code many times to get an accurate measurement
# - runs the whole thing three times and reports the best (accounts for other processes)
# - gives the answer in easy to read units


!python -m timeit "12345 in range(1000000)"

#In the above code, do we spend more time constructing the range list or checking if the number is in it? 

!python -m timeit "range(1000000)" 

#Yep, takes loads to time to construct the list. Seperate that bit out with a setup (-s) command:

!python -m timeit -s "r=range(1000000)" "12345 in r" 


#which is faster?
from __future__ import division 

def at_count(dna): 
    return (dna.count('a') + dna.count('t')) / len(dna) 
 
def at_iter(dna): 
    a_count = 0 
    t_count = 0 
    for base in dna: 
        if base == 'a': 
            a_count = a_count + 1 
        elif base == 't': 
            t_count = t_count + 1 
    return (a_count + t_count) / len(dna) 

test_dna = 'atcgatcgatcatgatcggatcgtagctagcatctagtc' 
assert(at_count(test_dna) == at_iter(test_dna)) 

import random
def random_dna(length):
    return "".join([random.choice(['A','T','G','C']) for _ in range(length)])
    
#Now we can compare the two functions:
%timeit at_count(random_dna(10000))
%timeit at_iter(random_dna(10000))

#but how much time is being used generating the random sequence?
d = random_dna(10000)
%timeit at_count(d)
%timeit at_iter(d)


#getting timing right is hard
#count() is faster than iteration (due to fast C code)


#Benchmarking memory - how much memory is your code using
pip install psutil

#in python:
import psutil, os

 #make a process object, then call os.getpid() on it
process = psutil.Process(os.getpid()) 

#then .get_memory_info().rss, divide by 1000 to get the kbs of memory being used by your script
mem = process.get_memory_info().rss / 1000 
print("Used this much memory: " + str(mem) + ' kb')

#This lets us investigate time/memory trade offs. We know that checking to see if a number is 
#in a set is faster than checking to see if it's in a list:

l = range(1000000)
s = set(l)
%timeit 12345 in l
%timeit 12345 in s

#1000 loops, best of 3: 154 µs per loop
#The slowest run took 19.15 times longer than the fastest. This could mean that an intermediate result is being cached 
#10000000 loops, best of 3: 112 ns per loop

%timeit list(range(1000000))
%timeit set(range(1000000))

#10 loops, best of 3: 56 ms per loop
#10 loops, best of 3: 113 ms per loop

!python list_mem.py
!python set_mem.py

#Used this much memory: 40349 kb
#Used this much memory: 65953 kb

#Conclusions:
# - if we need to create a list once then check membership many times, a set will be faster
# - if we need to create many lists, a set might be slower
# - a set will use more (x1.5) memory for these ranges


#Profiling is the process of taking an existing piece of code and identifying which bits are taking the time.

# create a random dna sequence
dna = random_dna(10000)

# create 100 random interesting motifs
motifs = [random_dna(4) for _ in range(100)]

%%timeit
# standard kmer counting code to identify frequent chunks
frequent_chunks = [] 
for start in range(len(dna) - 3): 
    chunk = dna[start:start + 4] 
    if dna.count(chunk) > 50: 
        frequent_chunks.append(chunk) 

# now check each chunk to see if it's in the list of motifs
for chunk in frequent_chunks: 
    if chunk in motifs: 
        print(chunk + " is frequent and interesting") 
    else: 
        print(chunk + " is frequent but not interesting")
        
#1 loops, best of 3: 451 ms per loop


#How can we speed this program up? 
#We know that checking to see if an element is in a list is slow, so let's change it to a set:

# create 100 random interesting motifs
motifs = set([random_dna(4) for _ in range(100)])

%%timeit
# standard kmer counting code to identify frequent chunks
frequent_chunks = [] 
for start in range(len(dna) - 3): 
    chunk = dna[start:start + 4] 
    if dna.count(chunk) > 50: 
        frequent_chunks.append(chunk) 

# now check each chunk to see if it's in the list of motifs
for chunk in frequent_chunks: 
    if chunk in motifs: 
        print(chunk + " is frequent and interesting") 
    else: 
        print(chunk + " is frequent but not interesting")
        
# loops, best of 3: 520 ms per loop

#Why did this not work? Probably because the line
#if chunk in motifs:
#doesn't actually get execute that often (maybe 600 times).

#cProfile is a built in module for profiling functions. 

import cProfile 

# we have to turn the code into a function so we can pass its name to run()
def classify_chunks():
    frequent_chunks = [] 
    for start in range(len(dna) - 3): 
        chunk = dna[start:start + 4] 
        if dna.count(chunk) > 50: 
            frequent_chunks.append(chunk) 

    for chunk in frequent_chunks: 
        if chunk in motifs: 
            print(chunk + " is frequent and interesting") 
        else: 
            print(chunk + " is frequent but not interesting")

cProfile.run("classify_chunks()")

#cProfile gives us tabular output:

#ncalls, which tells us how many times the function was called
#tottime, which tells us the total amount of time that was spent in that function (not including sub functions)
#percall, which tells us the amount of time that was spent in that function (not including sub functions) each time it was called
#cumtime, which is like tottime but does include sub functions
#another percall, which is like the first one except that it does include sub functions
#filename, which tells us the filename, line number, and name of the function or method

#cProfile only measures function calls, i.e. not stuff like
dna[start:start + 4] 
dna.count(chunk) > 50
chunk in motifs

#To use it well, we need structured code e.g. if we split our code into two functions:
def get_frequent_chunks(dna): 
    frequent_chunks = [] 
    for start in range(len(dna) - 3): 
        chunk = dna[start:start + 4] 
        if dna.count(chunk) > 50: 
            frequent_chunks.append(chunk) 
    return frequent_chunks 
 
def print_chunks(chunks): 
    for chunk in chunks: 
        if chunk in motifs: 
            print(chunk + " is frequent and interesting") 
        else: 
            print(chunk + " is frequent but not interesting") 
            
def classify_chunks(): 
    frequent_chunks = get_frequent_chunks(dna) 
    print_chunks(frequent_chunks) 

cProfile.run("classify_chunks()") 

#ncalls tottime percall   cumtime   percall filename:lineno(function) 
# 1    0.005    0.005      0.619     0.619  time_profile.py:14(get_frequent_chunks)
# 1    0.004    0.004      0.004     0.004  time_profile.py:22(print_chunks)


#Profiling with line_profiler

#line_profiler is a third part module that you have to install separately:
pip install line_profiler

#It measures execution time per line. To use it we add a @profile decorator to the function and run from the command line with:
kernprof -l -v slowfast.py

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    11                                           @profile
    12                                           def classify_chunks():
    13         1            3      3.0      0.0      frequent_chunks = [] 
    14      9998         8629      0.9      1.4      for start in range(len(dna) - 3): 
    15      9997        10739      1.1      1.7          chunk = dna[start:start + 4] 
    16      9997       601591     60.2     96.1          if dna.count(chunk) > 50: 
    17       318          394      1.2      0.1              frequent_chunks.append(chunk) 
    18                                           
    19       319          481      1.5      0.1      for chunk in frequent_chunks: 
    20       318         1043      3.3      0.2          if chunk in motifs: 
    21       109         1030      9.4      0.2              print(chunk + " is frequent and interesting") 
    22                                                   else: 
    23       209         2073      9.9      0.3              print(chunk + " is frequent but not interesting")

#"if dna.count(chunk) > 50" uses 96% of the time. Let's switch to a dict (rather than using count repeatidly) and keep a tally:

def classify_chunks(): 
    chunk_count = {}
    for start in range(len(dna) - 3):
        chunk = dna[start:start + 4]
        current_count = chunk_count.get(chunk, 0)
        new_count = current_count + 1 
        chunk_count[chunk] = new_count

    for chunk, count in chunk_count.items():
        if count > 50: 
            if chunk in motifs:
                print(chunk + " is frequent and interesting")
            else:
                print(chunk + " is frequent but not interesting")
                
Total time: 0.056549 s
File: chunks.py
Function: classify_chunks at line 11

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
    11                                           @profile
    12                                           def classify_chunks(): 
    13         1            2      2.0      0.0      chunk_count = {}
    14      9998         9893      1.0     17.5      for start in range(len(dna) - 3):
    15      9997        11790      1.2     20.8          chunk = dna[start:start + 4]
    16      9997        12889      1.3     22.8          current_count = chunk_count.get(chunk, 0)
    17      9997        10030      1.0     17.7          new_count = current_count + 1 
    18      9997        10919      1.1     19.3          chunk_count[chunk] = new_count
    19                                           
    20       257          341      1.3      0.6      for chunk, count in chunk_count.items():
    21       256          349      1.4      0.6          if count > 50: 
    22         9           35      3.9      0.1              if chunk in motifs:
    23         2          152     76.0      0.3                  print(chunk + " is frequent and interesting")
    24                                                       else:
    25         7          149     21.3      0.3                  print(chunk + " is frequent but not interesting")
    
    
#File and network IO is slow, so minimize it
#Existing modules are likely to be fast, so use them (scipy/numpy/pandas)
#Avoid calculating the same thing multiple times
#Functional structures (maps/comprehensions) tend to be faster than loops
#Pick data structures with the properties you want (list vs. set)

#and some advanced things to know about:
#write inline code in faster languages
#use parallel code
#let a database do the heavy lifting if you can
