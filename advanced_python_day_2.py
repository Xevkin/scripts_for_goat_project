#Sessions 3 & 4 - Object-orientated programming

#Everything is an object
#Each object is an instance of a specific class
#The class definition tells Python how instances of a class behave

#Use type() to get the type of an object

#dir() to see all the properties of an object

#anything with double underscores are not supposed to be called directly by
#the end user (__future__ is an exception) - it's supposed to be used internally by python

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
