#https://dl.dropboxusercontent.com/u/2838022/python_courses/eg_ap_2015/index.html

#Session 5 - Functional Programming

#A variable is mutable if it can change state while a program is running:

x = 0
for i in range(10):
    x = x + i
print(x)

#x changes state while the programme is running
#but we may want to avoid changing state:

x = sum(range(10))
print(x)

#classes (yesterday's exercises) involve many changes of state


#Side Effects - when a function changes the variable that is passed in

def my_function(i): 
    i.append('a') 
    return(i) 

foo = [1,2,3] 
print(foo) 
bar = my_function(foo)

#foo is returned and has been altered - need to look through every line of the function
#to find the cause of the side effect.
#Version without a side effect:

def my_function(i): 
    return(i + ['a']) 

my_function([1,2,3])

#If we have functions referring to each other, it can be very complicated to
#find the source of the side effect.

#A function may also depend on a variable outside the function, which can be changed (this is a silly thing to do)
#This means that the same function may be given the same input but give a different output

#Functions which:
# - don't have any side effects and
# - always return the same output for the same input
#are called pure functions.

#Functions are OBJECTS in python.


def print_list_with_function(my_list, my_function): 
    # the my_function argument is the name of a function
    for element in my_list: 
        print(my_function(element)) 

input = ['abc', 'defhij', 'kl'] 
print_list_with_function(input, len) 


#we are passing a variable to the function which represents a function
#we only pass "len" rather than "len()" to the function

#For one-line functions we can use a shortcut: lambda expressions:
get_second = lambda(input) : input[1] 

print_list_with_function(input, get_second)

#Or even:
print_list_with_function(input, lambda(input) : input[1] )

#^In this case we are creating an annonymous function and using it as an argument

#In non-functional (procedural?) code, we describe how to obtain the answer.
total = 0 
for i in range(11): 
    total = total + i 
print(total) 

#In functional code, we decribe the result we want and let the computer figure out how to calculate it
print(sum(range(11))) 

#another example
input = ['hello', 'world']

# how to get the answer we want
result1 = []
for i in input:
    result1.append(i[1]) 

# describe the answer we want (don't need brackets for lambda)
result2 = map(lambda x : x[1], input)
result1,result2

#Higher order functions: take the name of a function as an input
#For the general pattern of transforming each element of a list, use map(). 
#First argument is the name of a transformation function, either built in or defined.

lengths = map(len, dna_list) 

#using lambda:

at_contents = map( 
    lambda dna : (dna.count('A') + dna.count('T')) / len(dna), 
    dna_list 
)
at_contents

#map works on any iterable type (eg File object, would iterate over lines):

map(lambda x: x.lower(), 'ABCDEF')



#filter() takes a function argument which returns True or False (has to return True or False)

def is_long(dna): 
    return len(dna) > 5 

def is_at_poor(dna): 
    at = (dna.count('A') + dna.count('T')) / len(dna) 
    return at < 0.6 

long_dna = filter(is_long, dna_list) 
at_poor_dna = filter(is_at_poor, dna_list) 
print(long_dna)
print(at_poor_dna)


#By default the sorted() function sorts alphabetically
#For custom sorting, we pass in a key keyword argument which is the name of a transformation function.

sorted(dna_list, key=len) 

sorted(dna_list, key=len, reverse=True) 

#sorted() works by returning a copy of the list - or other iterable

sorted('atcgatcg')

#There is also list.sort() which behaves the same way but works by mutating the original list.
#key can be arbitrarily complex:

import re 
def poly_a_length(dna): 
    poly_a_match = re.search(r'A+$', dna) 
    if poly_a_match: 
        return len(poly_a_match.group()) 
    else: 
        return 0 

poly_a_length('ACGTGC')

dna_list = ['ATCGA', 'ACGG', 'CGTAAA', 'ATCGAA']
sorted(dna_list, key=poly_a_length)

#Higher order functions 
#Rather than writing two functions...
def get_6mers_at(dna): 
    result = [] 
    for i in range(len(dna) - 5): 
        one_6mer = dna[i:i+6] 
        at = (one_6mer.count('a') + one_6mer.count('t')) / 6 
        result.append(at) 
    return result 

def get_6mers_cg(dna): 
    result = [] 
    for i in range(len(dna) - 5): 
        one_6mer = dna[i:i+6] 
        cg = one_6mer.count('cg') 
        result.append(cg) 
    return result

print(get_6mers_at(dna))
print(get_6mers_cg(dna))

#....we can figure out what is different between the two functions, take it out and
#use it as an argument (which will be a function we run on each argument)

def get_at(dna): 
    return (dna.count('a') + dna.count('t')) / len(dna) 

def get_6mers_f(dna, analyze_6mer): 
    result = [] 
    for i in range(len(dna) - 5): 
        one_6mer = dna[i:i+6] 
        result.append(analyze_6mer(one_6mer)) 
    return result 

get_6mers_f(dna, get_at)

#we can do the cg dinucleotide part using a lambda expression:
get_6mers_f(dna, lambda dna : dna.count('cg'))


#Session Six - List Comprehensions

#Python has a special syntax for defining lists called list comprehensions. 
#Here's the list of lengths of the DNA sequences in three ways:

# with a loop

l1 = []
for dna in dna_list:
    l1.append(len(dna))
    
# with a map
l2 = map(len, dna_list)

# as a list comprehension
l3 = [len(dna) for dna in dna_list]

assert l1 == l2
assert l1 == l3


#comprehensions can have conditions:

[len(dna) for dna in dna_list if get_at(dna) >= 0.5]

#or use an annonymous function

