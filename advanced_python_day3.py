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

#map works on any iterable type:

map(lambda x: x.lower(), 'ABCDEF')
