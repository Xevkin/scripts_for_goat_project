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
