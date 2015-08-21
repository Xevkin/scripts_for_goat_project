#https://dl.dropboxusercontent.com/u/2838022/python_courses/eg_ap_2015/index.html

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

Why did this not work? Probably because the line
if chunk in motifs:
doesn't actually get execute that often:
