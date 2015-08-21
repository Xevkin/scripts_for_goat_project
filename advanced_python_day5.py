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
