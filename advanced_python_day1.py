#https://dl.dropboxusercontent.com/u/2838022/python_courses/eg_ap_2015/index.html

#Session 1 - Data Structures

#Data Structures
#Dicts
#Dicts of sets are an interesting data structure - can rapidly look up the key and value (eg condition plus genes expressed), \
#and also exploit the set nature e.g. check if one value is a  subset of another value

#Dicts of tuples - good for searching for specific records (list of tuples is fine for iteration but slow for the former), \
#as you can do something like this:

records = {
    'ABC123' : ('actgctagt', 1),
    'XYZ456' : ('ttaggttta', 1),
    'HIJ789' : ('cgcgatcgt', 5)
}

#iteration example

for accession, record in records.items():
    (sequence, code) = record
    print("looking at record " + accession + " with genetic code " + str(code))
    # do something with the record

#searching for something specific

my_record = records.get('XYZ456')
(this_sequence, this_code) = my_record
print("looking at record " + accession + " with genetic code " + str(code))

#^This works as the ascension number is being used as the key, rather than the record code - would have to \
#iterate if we were interested in that

#Collection data structure
import collections

dna = 'aattggaattggaattg'
base_counter = collections.Counter(dna)
print(base_counter)

collections.Counter is a special dict class 
 

#We may want to create dicts which have a default value i.e. an empty list
#We use collections.defaultdict to do this - it avoids us having to do something like creating an empty - {list_of_positions = kmer2list.get(kmer, [])}

dna = 'aattggaattggaattg'
k = 4 
kmer2list = {}
kmer2list = collections.defaultdict(list)
for start in range(len(dna) - k + 1): 
    kmer = dna[start:start + k] 
    kmer2list[kmer].append(start)
kmer2list



from __future__ import division

heavy_metals = {
    'arsenic' : {1,2,3,4,5,6,8,12},
    'cadmium' : {2,12,6,4},
    'copper' : {7,6,10,4,8},
    'mercury' : {3,2,4,5,1}
}


#input is  a dict of sets

def main(gene_sets):
    
    '''Gives the similarity score (pairwise similarity matrix) between gene expression conditions'''
    
    #create a  dictionary to store the output
    #each key is a condition
    #each value will be another dictionary
    #the keys of this dictionary will be every other condition
    #and the value for each is the similarity score for the two condition
    
    psm = {}
    
   #we have to iterate over all values in the dictionary,
   #and iterate again but ignore when the conditions are 
    
    for condition1, gene_set1 in gene_sets.items():
        
        #create a dictionary for condition1 in our final output
        psm[condition1] = {}
        
        #now iterate over every other conditions
        for condition2, gene_set2 in gene_sets.items():
            
            if condition1 != condition2:
    
                #calculate the similarity score
                #which is the intersect of the two sets divided by the union
                union_set = gene_set1.union(gene_set2)
               
                intersect_set = gene_set1.intersection(gene_set2)
    
                score = len(intersect_set) / len(union_set)
                
                #make a variable for our current condition, which is a
                #dictionary
                current_condition = psm[condition1]
                
                current_condition[condition2] = score
                
    
    return psm
    
    
main(heavy_metals)



##########################################################
#Session 2 - Recursions



#in writing
#Start with a list of all possible single-base sequences. Next, extend each sequence by adding each of the four 
#possible bases onto the end. Repeat this extension process as many times as necessary until the sequences are the length 
#you require

def generate_kmers(k): 
    # start with all the one-base strings
    result = ['A', 'T', 'G', 'C'] 

    # keep making the strings longer...
    for i in range(k - 1): 
        print("start of loop, result is", result)
        new_result = [] 
        
        
        # by taking each one and adding each possible base
        for kmer in result: 
            for base in ['A', 'T', 'G', 'C']: 
                new_result.append(kmer + base) 

        # the final list becomes the starting list for the next iteration
        ##note here - result is being altered as we build up the length of the mer, stopping
        ##when we reach k
        result = new_result 
        print("end of loop, result is", result)

    return result 

generate_kmers(3)

        

#Alternative way - see in writing
#To get a list of all kmers of a given length, start by checking the length. If the length is one then the result 
#is simply a list of the four bases. If the length is more than one, take the list of all possible sequences whose length 
#is one less that the length you're looking for, 
#and add each of the four possible bases to each of its elements to get the result.


def generate_kmers_rec(k): 
    # if k is one, then the result is all one-base strings
    if k == 1: 
        return ['A', 'T', 'G', 'C'] 
    
    # if k is bigger than one...
    else: 
        result = [] 
        
        # ...get a list of all kmers which are one base shorter...
        #this is the big line - the function references itself
        #it continually references itself until
        #note that the function will keep on reaching this part before calling itself
        #this means that the resuilts = [] will not clear what is returned
        #effectively, generate_kmers_rec(k - 1) will become a list that is what returned
        for seq in generate_kmers_rec(k - 1): 
        

            # ...and append each of the four possible bases
            for base in ['A', 'T', 'G', 'C']: 
                result.append(seq + base)
        return result 
    
generate_kmers_rec(3)



#So the function iteratively calls itself UNTIL IT GETS TO THE CONDITION THAT IT DOES NOT CALL ITSELF

tax_dict = { 
'Pan troglodytes' : 'Hominoidea',       'Pongo abelii' : 'Hominoidea', 
'Hominoidea' :  'Simiiformes',          'Simiiformes' : 'Haplorrhini', 
'Tarsius tarsier' : 'Tarsiiformes',     'Haplorrhini' : 'Primates',
'Tarsiiformes' : 'Haplorrhini',         'Loris tardigradus' : 'Lorisidae',
'Lorisidae' : 'Strepsirrhini',          'Strepsirrhini' : 'Primates',
'Allocebus trichotis' : 'Lemuriformes', 'Lemuriformes' : 'Strepsirrhini',
'Galago alleni' : 'Lorisiformes',       'Lorisiformes' : 'Strepsirrhini',
'Galago moholi' : ' Lorisiformes'
} 

def get_ancestors(taxon):
    if taxon == 'Primates':
        return []
    else:
        parent = tax_dict.get(taxon)
        parent_ancestors = get_ancestors(parent) 
        return [parent] + parent_ancestors
        #so each iteration will return all previous ancestors plus the parent taxon of the current taxon
        
        
new_tax_dict = { 
    'Primates': ['Haplorrhini', 'Strepsirrhini'], 
    'Tarsiiformes': ['Tarsius tarsier'], 
    'Haplorrhini': ['Tarsiiformes', 'Simiiformes'], 
    'Simiiformes': ['Hominoidea'], 
    'Lorisidae': ['Loris tardigradus'], 
    'Lemuriformes': ['Allocebus trichotis'], 
    'Lorisiformes': ['Galago alleni','Galago moholi'], 
    'Hominoidea': ['Pongo abelii', 'Pan troglodytes'], 
    'Strepsirrhini': ['Lorisidae', 'Lemuriformes', 'Lorisiformes'] 
} 

#two ways of writing a script to get all the descendents of a given taxon - iterative and recursive

def get_descendants(taxon): 
    result = [] 

    # the stack is the list of taxa whose children we need to include
    stack = [taxon] 

    while len(stack) != 0: 

        # remove the last element from the stack
        current_taxon = stack.pop() 

        # look up the children
        current_taxon_children = new_tax_dict.get(current_taxon, []) 
        
        # add the children onto the end of the stack
        stack.extend(current_taxon_children) 

        # add the children to the result list
        result.append(current_taxon) 

    # when the stack is empty we are done
    return result 

get_descendants("Simiiformes")




def get_descendants_rec(taxon): 

    # add the current taxon to the resul
    result = [taxon] 

    # look up the list of children for this taxon
    children = new_tax_dict.get(taxon, []) 

    #if the children list is empty, this step is skipped
    #result ill still contain the taxon so it will return that
    # for each child taxon...
    for child in children: 

        # ... add its decendnts to the result
        result.extend(get_descendants_rec(child)) 

    return result 

get_descendants_rec("Haplorrhini")

#LCA exercise using iterative answer:

tax_dict = { 
'Pan troglodytes' : 'Hominoidea',       'Pongo abelii' : 'Hominoidea', 
'Hominoidea' :  'Simiiformes',          'Simiiformes' : 'Haplorrhini', 
'Tarsius tarsier' : 'Tarsiiformes',     'Haplorrhini' : 'Primates',
'Tarsiiformes' : 'Haplorrhini',         'Loris tardigradus' : 'Lorisidae',
'Lorisidae' : 'Strepsirrhini',          'Strepsirrhini' : 'Primates',
'Allocebus trichotis' : 'Lemuriformes', 'Lemuriformes' : 'Strepsirrhini',
'Galago alleni' : 'Lorisiformes',       'Lorisiformes' : 'Strepsirrhini',
'Galago moholi' : ' Lorisiformes'
} 

#set up the function to get a list of the ancestors of a taxon
def get_ancestors(taxon):
    result = [taxon] 
    while taxon != 'Primates' and taxon != None:
        parent = tax_dict.get(taxon) 
        result.append(parent)
        taxon = parent
    return result


#input is a list of taxons
def get_LCA(taxon_list):
    
    while (len(taxon_list) >= 2) :
        
        ancestors1 = get_ancestors(taxon_list[0])
        
        ancestors2 = get_ancestors(taxon_list[1])
        
        taxon_list.pop(0)
        taxon_list.pop(0)
    
        counter = 0
    
        while (counter != 1):
    
            last_anc = ancestors1[0]
            
            ancestors1 = ancestors1[1:]
          
            if last_anc in ancestors2:

                    counter = 1
                
        taxon_list.append(last_anc)
        
    return last_anc

        
print get_LCA([ "Hominoidea", "Pan troglodytes","Tarsius tarsier"])


#recursive answer

tax_dict = { 
'Pan troglodytes' : 'Hominoidea',       'Pongo abelii' : 'Hominoidea', 
'Hominoidea' :  'Simiiformes',          'Simiiformes' : 'Haplorrhini', 
'Tarsius tarsier' : 'Tarsiiformes',     'Haplorrhini' : 'Primates',
'Tarsiiformes' : 'Haplorrhini',         'Loris tardigradus' : 'Lorisidae',
'Lorisidae' : 'Strepsirrhini',          'Strepsirrhini' : 'Primates',
'Allocebus trichotis' : 'Lemuriformes', 'Lemuriformes' : 'Strepsirrhini',
'Galago alleni' : 'Lorisiformes',       'Lorisiformes' : 'Strepsirrhini',
'Galago moholi' : ' Lorisiformes'
} 

#set up the function to get a list of the ancestors of a taxon
def get_ancestors(taxon):
    result = [taxon] 
    while taxon != 'Primates' and taxon != None:
        parent = tax_dict.get(taxon) 
        result.append(parent)
        taxon = parent
    return result

def get_LCA(taxa1, taxa2):
        
        ancestors1 = get_ancestors(taxa1)
        
        ancestors2 = get_ancestors(taxa2)
    
        counter = 0
    
        while (counter != 1):
    
            last_anc = ancestors1[0]
            
            ancestors1 = ancestors1[1:]
          
            if last_anc in ancestors2:

                    counter = 1
        
        return last_anc


#input is a list of taxons
def recursive(taxon_list):
    
    if (len(taxon_list) == 2):
    
        return get_LCA(taxon_list[0], taxon_list[1])
    
    else:
        
         return recursive(taxon_list[1:])
    
    
    
    
