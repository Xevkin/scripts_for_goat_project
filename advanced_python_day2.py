#https://dl.dropboxusercontent.com/u/2838022/python_courses/eg_ap_2015/index.html

#Sessions 3 & 4 - Object-orientated programming

#Everything is an object
#Each object is an instance of a specific class
#The class definition tells Python how instances of a class behave

#Use type() to get the type of an object

#dir() to see all the properties of an object

#anything with double underscores are not supposed to be called directly by
#the end user (__future__ is an exception) - it's supposed to be used internally by python

#define a class using "class" which here is a type of "object"
#classes traditionally start with a capital letter

class DNARecord(object): 
    #the class contains individuals variable which can be defined
    #note that in this example the variables are constant in all instances of the class but can be changed
    sequence = 'ACGTAGCTGACGATC'
    gene_name = 'ABC1'
    species_name = 'Drosophila melanogaster'

    #methods are defined here - written like functions
    #self refers to itself - the object that is being 
    def complement(self): 
        replacement1 = self.sequence.replace('A', 't') 
        replacement2 = replacement1.replace('T', 'a') 
        replacement3 = replacement2.replace('C', 'g') 
        replacement4 = replacement3.replace('G', 'c') 
        return replacement4.upper() 

    def get_AT(self): 
        length = len(self.sequence) 
        a_count = self.sequence.count('A') 
        t_count = self.sequence.count('T') 
        at_content = (a_count + t_count) / length 
        return at_content 
#this has only DEFINED what a DNA record is
#we need to create an instance of the class
d = DNARecord()
print('Created a record for ' + d.gene_name + ' from ' + d.species_name)
print('AT is ' + str(d.get_AT()))
print('complement is ' + d.complement())

#One major problem: all the variables are set as part of the class definiion, so every object we create 
#will have the same sequence, etc. 
#We can change the object variables after creation just like any variable:

d1 = DNARecord() 
d1.sequence = 'ATATATTATTATATTATA' 
d1.gene_name = 'COX1' 
d1.species_name = 'Homo sapiens' 
 
d2 = DNARecord() 
d2.sequence = 'CGGCGGCGCGGCGCGGCG' 
d2.gene_name = 'ATP6' 
d2.species_name = 'Gorilla gorilla' 
 
for r in [d1, d2]: 
    print('Created ' + r.gene_name + ' from ' + r.species_name) 
    print('AT is ' + str(r.get_AT())) 
    print('complement is ' + r.complement())
    
#We can create a method that allows us to easily change the variables of the object 
    
class DNARecord(object): 
    sequence = 'ACGTAGCTGACGATC'
    gene_name = 'ABC1'
    species_name = 'Drosophila melanogaster'

    def complement(self): 
        replacement1 = self.sequence.replace('A', 't') 
        replacement2 = replacement1.replace('T', 'a') 
        replacement3 = replacement2.replace('C', 'g') 
        replacement4 = replacement3.replace('G', 'c') 
        return replacement4.upper() 

    def get_AT(self): 
        length = len(self.sequence) 
        a_count = self.sequence.count('A') 
        t_count = self.sequence.count('T') 
        at_content = (a_count + t_count) / length 
        return at_content 
    
    def set_variables(self, new_seq, new_gene_name, new_species_name): 
        self.sequence = new_seq 
        self.gene_name = new_gene_name 
        self.species_name = new_species_name     
Now we can do this:
d1 = DNARecord() 
d1.set_variables('ATATATTATTATATTATA','COX1','Homo sapiens') 

#we can now remove the preset variables that are part of the class definition
#this can be a problem if you forget to set the variables
#We need a constructor: a special method whose job is to create an object and set up the variables. 
#The constructor has a special name __init__() 
#short for initialize

class DNARecord(object): 
    
    def __init__(self, sequence, gene_name, species_name):
        self.sequence = sequence
        self.gene_name = gene_name
        self.species_name = species_name
        
    def complement(self): 
        replacement1 = self.sequence.replace('A', 't') 
        replacement2 = replacement1.replace('T', 'a') 
        replacement3 = replacement2.replace('C', 'g') 
        replacement4 = replacement3.replace('G', 'c') 
        return replacement4.upper() 

    def get_AT(self): 
        length = len(self.sequence) 
        a_count = self.sequence.count('A') 
        t_count = self.sequence.count('T') 
        at_content = (a_count + t_count) / length 
        return at_content        
        
d1 = DNARecord('ATATATTATTATATTATA', 'COX1', 'Homo sapiens')
print(d1.complement())

#__init__() will prevent us from setting up an instance of the class
#without the required variables
#we can set defaults, and allow them to be changed but the one to be changed/the optional
#argument must be the last variable when the instance is created

#we can also make classes which are derived from other classes (saves us from writing methods, etc repeatidy)
class SequenceRecord(object): 

    def __init__(self, sequence, gene_name, species_name): 
        self.sequence = sequence 
        self.gene_name = gene_name 
        self.species_name = species_name 

    def get_fasta(self): 
        safe_species_name = self.species_name.replace(' ','_') 
        header = '>' + self.gene_name + '_' + safe_species_name 
        return header + '\n' + self.sequence + '\n' 
        
        
class ProteinRecord(SequenceRecord): 
    
    def get_hydrophobic(self): 
        aa_list=['A','I','L','M','F','W','Y','V'] 
        protein_length = len(self.sequence) 
        total = 0 
        for aa in aa_list: 
            aa = aa.upper() 
            aa_count = self.sequence.count(aa) 
            total = total + aa_count 
        return total * 100 / protein_length  

class DNARecord(SequenceRecord): 

    def complement(self): 
        replacement1 = self.sequence.replace('A', 't') 
        replacement2 = replacement1.replace('T', 'a') 
        replacement3 = replacement2.replace('C', 'g') 
        replacement4 = replacement3.replace('G', 'c') 
        return replacement4.upper() 

    def get_AT(self): 
        length = len(self.sequence) 
        a_count = self.sequence.count('A') 
        t_count = self.sequence.count('T') 
        return (a_count + t_count) / length 
        
#We can also make constructors in methods which override the constructor in the parent method
class DNARecord(SequenceRecord): 
    
    def __init__(self, sequence, gene_name, species_name, genetic_code): 
        self.sequence = sequence 
        self.gene_name = gene_name 
        self.species_name = species_name 
        self.genetic_code = genetic_code 

    def complement(self): 
        replacement1 = self.sequence.replace('A', 't') 
        replacement2 = replacement1.replace('T', 'a') 
        replacement3 = replacement2.replace('C', 'g') 
        replacement4 = replacement3.replace('G', 'c') 
        return replacement4.upper() 

    def get_AT(self): 
        length = len(self.sequence) 
        a_count = self.sequence.count('A') 
        t_count = self.sequence.count('T') 
        return (a_count + t_count) / length 
        

#we can call methods (and constructors)  from the superclass eg
class SequenceRecord(object): 

    def __init__(self, sequence, gene_name, species_name): 
        if not re.match(r'[A-Z][a-z]+ [a-z]+', species_name): 
            raise ValueError(species_name + ' is not a valid species name!')
        self.sequence = sequence 
        self.gene_name = gene_name 
        self.species_name = species_name     
        
class DNARecord(SequenceRecord): 
    
    def __init__(self, sequence, gene_name, species_name, genetic_code): 
        # first call the SequenceRecord constructor to check the species name
        SequenceRecord.__init__(self, sequence, gene_name, species_name) 
        # now set the genetic code 
        self.genetic_code = genetic_code 
        
#this allows us to use the error checking from SequenceRecord that would be otherwise overwritten
#in DNARecord as it uses its own constructor (but instead we have called a method)

#we can also use the same name for different methods - this is called polymorphism

class ProteinRecord(SequenceRecord): 
    
    def get_protein_length(self): 
        return len(self.sequence) 
    
class DNARecord(SequenceRecord): 

    def get_protein_length(self): 
        return len(self.sequence) / 3    
        
        
for my_record in list_of_records:
# we don't care whether it's a DNA or protein record
    if my_record.get_protein_length() > 100:
        # do something with the record
        
#Project for the class:

from __future__ import division
import random

class Allele(object):
    
    
        def __init__(self, name, fitness):
            
            self.name = name
            self.fitness = fitness
            
            

class Locus(object):
    
    
    def __init__(self, name, alleles):
        
        self.name = name
        self.alleles = alleles
  
              
    def random_allele(self):
        
        return self.alleles[random.randint(0, len(self.alleles) -1)]
        
      
          
class Individual(object):
    
    
    def __init__(self, alleles):
        
        self.alleles = alleles
       
         
    def get_genotype(self):
        
        genotype = ""
        
        for i in self.alleles:
            genotype += i.name
            
        return genotype
    
    
    def get_fitness(self):
    
        fitness = 1
    
        for i in self.alleles:
    
            fitness = fitness * i.fitness
            
        return fitness
        
    
    
class Population(object):
    
    
    def __init__(self, name, individuals, loci_list):
        
        self.name = name
        
        self.individuals = individuals
        
        self.loci_list = loci_list
        
        
    def grow_pop(self, growth_size):
        
        for i in range(0, growth_size):
           
            self.individuals.append(random_ind(self.loci_list))
            
        return
        
        
    def alleles_in_pop(self):
        
        allele_list = []
        
        for i in self.loci_list:
            
            for j in i.alleles:
                
                if j not in allele_list:
                    
                    allele_list.append(j)
                    
        return allele_list
        
        
    def get_genotypes(self):
        
        population_genotypes = []
        
        for i in self.individuals:
            
            population_genotypes.append(i.get_genotype())
            
        return population_genotypes
       
         
    def get_genotypes_and_fitness(self):
        
        genotypes_and_fitness = []
        
        for i in self.individuals:
            
            genotypes_and_fitness.append([i.get_genotype(), i.get_fitness()])
            
        return genotypes_and_fitness
        
        
    def get_allele_freq(self, allele):
        
        population_counter = 0
        
        allele_counter = 0
        
        for i in self.individuals:
                        
            population_counter += 1
            
            if allele in i.alleles: 
           
                allele_counter += 1
        
        freq_in_pop = allele_counter / population_counter        
                                
        return freq_in_pop
        
    
    def allele_frequencies(self):
        
        allele_and_freq = []
        
        for i in self.alleles_in_pop():
            
            population_counter = 0
        
            allele_counter = 0
        
            for j in self.individuals:
                        
                population_counter += 1
            
                if i in j.alleles: 
           
                    allele_counter += 1
        
            freq_in_pop = allele_counter / population_counter
            
            allele_and_freq.append([i.name, freq_in_pop])
            
        return allele_and_freq
        
        
    def death(self):
        
        for i in self.individuals:
            
            
            #I am going to arbitrarily remove individuals here,
            #with a probability inversely related to fitness
            #I am going to add a random number between 0 and 9 to the inverse,
            #and if that number is greater than 10, I remove the individual
            #from the population
            if ((1 / i.get_fitness()) + random.randint(0,9) > 10):
                
                self.individuals.remove(i)
                
    
    def birth(self):
        
        
        alleles_in_individual = []
        
        for i in self.loci_list:
              
            x = i.random_allele()
              
            if x in self.alleles_in_pop():
                  
                  alleles_in_individual.append(x)
                  
            else:
                
                    i.alleles.remove(x)
                    
                    alleles_in_individual.append(i.alleles[0])
                  
        self.individuals.append(Individual(alleles_in_individual))
          
        return 
        
        
    def stabilize(self):
        
        while ( len(self.individuals) < 100):
            
            self.birth()
            
    
    def generation(self):
        
        self.death()
        
        self.stabilize()
            
            
                         
allele_A = Allele("A",1)
allele_a = Allele("a",0.75)

locus_A = Locus("locus_A", [allele_A, allele_a])
#we can create annonymous variables - they have no names
#in this case we are creating Allele objects
locus_B = Locus("locus_B", [Allele("B",1), Allele("b", 0.5)])
locus_C = Locus("locus_C", [Allele("C",1), Allele("c", 0.25)])

#we create allele objects once, but then refer to them repeatidly
#much more memory efficient
first_allele = locus_A.alleles[0]
second_allele = locus_A.alleles[1]
third_allele = locus_B.alleles[0]
fourth_allele = locus_B.alleles[1]
fifth_allele = locus_C.alleles[0]
sixth_allele = locus_C.alleles[1]

ind1 = Individual([first_allele,third_allele, fifth_allele])
ind2 = Individual([first_allele,third_allele, fifth_allele])

pop1 = Population("population_one", [ind1, ind2], [locus_A, locus_B, locus_C])


#create a function to return random alleles from a list of loci
def random_ind(loci_list):
    
    random_alleles = []
    
    for i in loci_list:
        
        random_alleles.append(i.random_allele())
    
    ind = Individual(random_alleles)
    
    return ind
    


    
ind3 = random_ind([locus_A, locus_B, locus_C])

pop2 = Population("population_two", [], [locus_A, locus_B, locus_C]) 

pop2.grow_pop(100)

def alleles_over_time(generations, output_file):
    
    file = open(output_file, "w")
        
    file.write("generation, A, a, B, b, C, c" + "\n")
        
    l = []
        
    for i in pop2.allele_frequencies():
           
        for j in i:
            
            if isinstance(j, float):    
                        
                l.append(str(j))
            
    
    file.write("0, " + ",".join(l) + "\n")
        
    for i in range(1, generations + 1):
        
        l = []    
                    
        pop2.generation()    
            
        for k in pop2.allele_frequencies():
            
           for j in k:
                               
                if isinstance(j, float):
                
                    l.append(str(j))
                    
        e = str(i) + "," + ",".join(l) + "\n"
        print e        
        file.write(e)
    
    file.close()
    
    
            
alleles_over_time(100, "test.csv")


#In OO code, we ask the object for the answer we want, and the object is responsible for figuring out how to calculate it
#the data is stored in the same place as the object definitions - useful for larger programmes/project
