#https://dl.dropboxusercontent.com/u/2838022/python_courses/eg_ap_2015/index.html

#Session 7 - Exception Handling

print("starting..")
print('abc' + 3)

dna = 'atctgcatattgcgtctgatg'
a_count = dna.count('A') 
print(a_count)

#^a bug - code works but is not donig what we want
#^these mistakes are all repeatable - we can't really recover from them

o = open('missingfile.txt')

#this error isn't determined by the code, but rather by the environment
#we could also write some code to try to fix this e.g. if the file doesn't exist we could
#ask the user to input a new file

#When something bad like this happens it's called an exception. 
#When writing code, we can decide what to do when an exception occurs.

#The default way to deal with exceptions is as above: do nothing and let Python print an error message. 
#If we want to actually do something based on the exception, wrap the code that might cause the exception 
#in a try block and put the exception-handling code in an except block.

try: 
    f = open('misssing.txt') 
    print('file contents: ' + f.read())
except: 
    print("Sorry, couldn't find the file you asked for") 
    
#Because we have handled or caught the exception, the program continues running

#Problem: except works for any error/exception

#we can catch multiple different types of exceptions with multiple blocks:
try: 
    f = open('my_file.txt') 
    my_number = int(f.read()) 
    print(my_number + 5) 
except IOError: 
    print("sorry, couldn't find the file") 
    # fix the problem somehow...
except ValueError: 
    print("sorry, couldn't parse the number") 
    # fix the problem somehow...
    
#For multiple exception types, use a tuple:
try: 
    f = open('my_file.txt') 
    my_number = int(f.read()) 
    print(my_number + 5) 
except (IOError, ValueError): 
    print("sorry, something went wrong") 

#we can find out more about the type of exceptions being raised:

try: 
    f = open('missing.txt') 
    my_number = int(f.read()) 
    print(my_number + 5) 
except IOError as ex: 
  #these exceptions are classes - not the capital letter
  #ex.strerror will contain the specific message produced by python for this error
  #not all exceptions have this .strerror - IOError does tho
    print("sorry, couldn't open the file: " + ex.strerror) 
except ValueError: 
    print("sorry, couldn't parse the number") 
    
#Where to put the line that prints the number (more generally: the code that relies on the lines that might raise an 
#exception)? With the version #above, print(my_number + 5) might also raise an IOError, so it's not a great idea to 
#have it inside the try block. Could put it outside:

try: 
    f = open('my_file.txt') 
    my_number = int(f.read()) 
except IOError as ex: 
    print("sorry, couldn't find the file: " + ex.strerror) 
except ValueError as ex: 
    print("sorry, couldn't parse the number: " +  ex.args[0]) 
print(my_number + 5) 

#but there's no point trying to print the number if it hasn't been sucessfully read from the file. 
#Solution: use an else block:

try: 
    f = open('my_file.txt') 
    my_number = int(f.read()) 
except IOError as ex: 
    print("sorry, couldn't find the file: " + ex.strerror) 
except ValueError as ex: 
    print("sorry, couldn't parse the number: " +  ex.args[0]) 
else:
    print(my_number + 5) 
    
#What if there's code that needs to run regardless of whether there was an exception or not? Consider this:
import os 

# write some temporary data to a file
t = open('temp.txt', 'w') 
t.write('some important temporary text') 
t.close() 

# do some other processing
f = open('my_file.txt') 
my_number = int(f.read()) 
print(my_number + 5) 

# delete the temporary file
os.remove('temp.txt') 

#When the exception is raised by int() the program exits and so the temp file does not get cleaned up. Where should we put 
#os.remove() if we want to make sure it always runs? Not using else, because else only runs in the absence of errors. 
#Also not at the end of the code, because it won't run if an exception is raised inside the try but not caught (anything 
#other than IOError or ValueError). Solution: finally blocks are always run:

import os 
t = open('temp.txt', 'w') 
t.write('some important temporary text') 
t.close() 
try: 
    f = open('my_file.txt') 
    my_number = int(f.read()) 
    print(my_number + 5) 
except IOError as ex: 
    print("sorry, couldn't find the file: " + ex.strerror) 
except ValueError as ex: 
    print("sorry, couldn't parse the number: " +  ex.args[0]) 
finally: 
    os.remove('temp.txt')
    
#finally blocks are useful for doing clean up code (files, network connections, database connections, logging, etc.).

try:
    # code in here will be run until an exception is raised
except ExceptionTypeOne:
    # code in here will be run if an ExceptionTypeOne
    # is raised in the try block
except ExceptionTypeTwo:
    # code in here will be run if an ExceptionTypeTwo 
    # is raised in the try block
else:
    # code in here will be run after the try block 
    # if it doesn't raise an exception
finally:
    # code in here will always be run
    
#you can find the exception types in the python docs

#Context managers - with(<my_file>) as file basically contains the IOError, try, finally etc (do not need to use try)

# Exceptions bubble up - handle exceptions in the place where your program can do something about it.
#lets say we have two functions, one which refers to the other
#we will have three opportunities to catch an exception in the first function - the first function itself,
#the second function and the main code

def function_one():
    try:
        # do some processing...
        return 5
    except SomeException:
        print("Handling exception")
        # handle the exception...
        
def function_two():
    my_number = function_one()
    return my_number + 2

print(function_two())


def function_one():
    # do some processing...
    return 5

def function_two():
    try:
        my_number = function_one()
        return my_number + 2
    except SomeException:
        print("Handling exception")
        # handle the exception...

        
print(function_two())


def function_one():
    # do some processing...
    return 5

def function_two():
    my_number = function_one()
    return my_number + 2
try:
    print(function_two())
except SomeException:
    print("Handling exception")
    # handle the exception
    
    
#the best place to have handle the exception is where it can be dealt with


#An exception is a signal that something has gone wrong. As well as responding to these signals, our code can create them. 

from __future__ import division
import re 
def get_at_content(dna): 
    if re.search(r'[^ATGC]', dna): 
        raise ValueError('Sequence cannot contain non-ATGC bases') 
    length = len(dna) 
    a_count = dna.count('A') 
    t_count = dna.count('T') 
    at_content = (a_count + t_count) / length 
    return at_content 

print(get_at_content('ATCGCTGTTATCGACTGACT'))
print(get_at_content('ATCGCTGANCGACTGATTCT'))

#Now we can make use of this. For example, given a large collection of sequences we don't want a single "bad" sequence to 
#cause the whole program to crash (as if we iterate the function over many lines, it will return an exception and
#exit once we hit an error (not ATGC base).
#We wrap the function in a try block:

for seq in sequences: 
    try: 
        print('AT content for ' + seq + ' is ' + str(get_at_content(seq)))
    except ValueError: 
        print('skipping invalid sequence '+ seq) 
        
#however, ValueError will be raised for lots of issues - we can look at the error message to distinguish:
for seq in sequences: 
    try: 
        print('AT content for ' + seq + ' is ' + str(get_at_content(seq)))
        #next line will create an error
        number = int('five')
    except ValueError as ex: 
        print('something went wrong with sequence '+ seq) 
        print("sorry, couldn't parse the number: " +  ex.args[0]) 

#But that doesn't really help us to recover from the error. We need a custom exception to signal a specific error:

class AmbiguousBaseError(Exception): 
    pass 
  
#in this case, the only point of the class is so we can recognize its name

def get_at_content(dna): 
    if re.search(r'[^ATGC]', dna): 
        raise AmbiguousBaseError('Sequence cannot contain non-ATGC bases') 
        
    #could have a pile more exceptions here
    
    length = len(dna) 
    a_count = dna.count('A') 
    t_count = dna.count('T') 
    at_content = (a_count + t_count) / length 
    return at_content 
 
sequences = ['ACGTACGTGAC', 'ACTGCTNAACT', 'ATGGCGCTAGC'] 
for seq in sequences: 
    try: 
        print('AT content for ' + seq + ' is ' + str(get_at_content(seq)))
    except AmbiguousBaseError: 
        print('skipping invalid sequence '+ seq) 
        
#Now we will only catch AmbiguousBaseError, any other exception can be dealt with separately.

#can put exception handling in the constructor of classes



#Session 8 - Packaging and Distribution
#packages are objects in python

import random

#we can find where the object is:
random.__file__

#.pyc means a compiled Python file, but if we drop the c from the end and open up the .py 
#we can prove that it's just normal Python code.

#Notice that the name of the module is just the name of the file without the extension.

#How does Python know what folders to look for the file in?
import sys
print sys.path

# pip, the Python package manager - remove old modules, updates them, ensures dependencies are met, installs
#extra command-line tools etc

pip install mymodule pip install --upgrade mymodule pip uninstall mymodule


#Creating a new module

from __future__ import division 
 
# calculate the AT content 
def calculate_at(dna): 
    length = len(dna) 
    a_count = dna.count('A') 
    t_count = dna.count('T') 
    at_content = (a_count + t_count) / length 
    return at_content 

#Move it to a new file in the current directory, call it at_calculator.py

import at_calculator 

at_calculator.__file__

#We could copy it to one of the path folders to use it universally (but we will see a better way).
 
# ask the user for a DNA sequence and filename 
dna = raw_input("Enter a DNA sequence:\n")
output_filename = raw_input("Enter a filename:\n")
 
# write the AT content to the output file 
with open(output_filename, "w") as out: 
    out.write(str(at_calculator.calculate_at(dna))) 



#Notice how we need to call modulename.functionname() i.e. at_calculator.calculate_at() . Annoying, but good because it 
#means each module has its own namespace - authors (including you!) don't have to worry about what names might have been 
#used by other people. #For long names this gets anoying:

import some_incredibly_long_awkward_to_type_module_name
some_incredibly_long_awkward_to_type_module_name.foo()

#We can get around it by using an alias:

import some_incredibly_long_awkward_to_type_module_name as bob
bob.foo()

#If we're really sure that names don't clash we can import the function directly:

from at_calculator import calculate_at
calculate_at("ACTGATCGTCGAT")

#If we know that names do clash, we can import functions with aliases:

from at_calculator import calculate_at as at1
from another_package import calculate_at as at2


#Python's docstring:  single string literal as the first line of a module or function, triple quotes let us 
#run over multiple lines.

"""Functions for calculating metrics of DNA sequences"""

from __future__ import division 

def calculate_at(dna): 
    """Return the AT content of the argument. 
    Only works for uppercase DNA sequences
    """
    
    length = len(dna)
    a_count = dna.count('A') 
    t_count = dna.count('T') 
    at_content = (a_count + t_count) / length 
    return at_content 
    
#Sometimes we want to have a Python file that can act either as a module, or as a program:
# - if there's currently only a single program that uses the module
# - if it's handy to have a small demonstration program for the module
# - to include a program which tests the module

#A Python file can access its own name using the __name__ variable.
# - If the file is being run as a program, then the __name__ variable is __main__
# - If the file is being imported, then the __name__ variable is the name of the file minus .py

if __name__ == "__main__":
    print("I am being run as a script!")   
    # demo program code goes here
    
    
#for example:    

from __future__ import division 

# calculate the AT content 
def calculate_at(dna): 
    length = len(dna) 
    a_count = dna.count('A') 
    t_count = dna.count('T') 
    at_content = (a_count + t_count) / length 
    return at_content 

if __name__ == "__main__": 
    dna = raw_input("Enter a DNA sequence:\n").rstrip("\n") 
    print("AT content is " + str(calculate_at(dna))) 
    

#Some modules require other code:

def translate_dna(dna): 
    """Return the translation of a DNA sequence"""

    # define a dict to hold the genetic code 
    gencode = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'} 

    last_codon_start = len(dna) - 2 
    protein = "" 

    # for each codon in the dna, append an amino acid to the protein
    for start in range(0,last_codon_start,3): 
        codon = dna[start:start+3] 
        amino_acid = gencode.get(codon.upper(), 'X') 
        protein = protein + amino_acid

    return protein 
    
#Inefficient - builds dic each time the funciton is used
#instead, we can define the variant outside the function within the module (initialize)
#this will be run when the module is imported

# define a dict to hold the genetic code 
gencode = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'} 

def translate_dna(dna): 
    last_codon_start = len(dna) - 2 
    protein = "" 

    # for each codon in the dna, append an amino acid to the protein
    for start in range(0,last_codon_start,3): 
        codon = dna[start:start+3] 
        amino_acid = gencode.get(codon.upper(), 'X') 
        protein = protein + amino_acid

    return protein 
    
#the gencode dict here could be manually accessed using module.dict
#We can imagine more sophisticated solutions here e.g. build the dict when first needed then reuse.


#A package gathers together related useful modules, just like a module gathers together related useful functions/classes. 
#Also, we need to build packages if we want to distribute code, even if it's just one module. 
#To make a package we create a folder with the package name, and add a special file called __init__.py inside it:

dnatools/
    __init__.py
    dna_translation.py
    
import dnatools.dna_translation
print(dnatools.dna_translation.translate_dna("ACTGTGAC"))

#Distributing a package:
# - to distribute a program that uses built in modules, just send the .py file
# - to distribute a program that uses third party modules, send the .py file and tell the end use what modules to install
# - to distribute a program plus separate module, send a zipped folder
# - for anything more complicated, use pip
