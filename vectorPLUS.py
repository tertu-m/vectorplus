import textwrap

######
# START IMPORTANT METHODS
######

 # Returns a list (see below) of all matching sequences in the given scaffold
 # with a given initial search length.
 # [name, (start_index, end_index)]
def tryFindStartCodon(scaffold, sequence, initialSearchLength=10):
    searchLength = initialSearchLength - 1
    outputList = []
    workList = []
    while True:
        searchLength+=1
        outputList = workList
        workList = []
        for start in substrings(scaffold,sequence[:searchLength]):
            workList.append((start,start+searchLength))
        if len(workList)==0:
            outputList.insert(0,searchLength-1)
            return outputList

 # Returns a codon sequence for a given RNA sequence
def to_codon(seq):
    cDict = codon_dict()
    output = ""
    s = split_by_n(seq, 3)
    for codon in s:
        if cDict.get(codon) != None:
            if cDict[codon] == "Stop":
                break
            else:
                output += cDict[codon]
    return output

 # Transcribes DNA into RNA, returning an RNA seq
def transcribe_DNA(seq):
    rna = ''
    for c in seq:
        if c == 'T':
            rna += 'U'
        else:
            rna += c
    return rna

 # Transcribes RNA back into DNA
def RNA_to_DNA(seq):
    dna = ''
    for c in seq:
        if c == 'U':
            dna += 'T'
        else:
            dna += c
    return dna

 # Takes a list of (title, nucleotide) tuples and prints as a FASTA file to the given filename
def buildFastaFile(listOfTuples,fileName):
    output = ""
    if fileName.find(".txt",0) == -1:
        fileName += ".txt"
    for data in listOfTuples:
        output += (">" + data[0] + "\n")
        output += textwrap.fill(data[1], 80)
    with open(fileName, 'w') as outFile:
        outFile.write(output)

######
# START UTILITY METHODS
######

#Poor replacement for the BioPython generator written by Jack
def read_fasta(fileHandle):
    title = None
    data = None
    for line in fileHandle:
        if line[0]==">":
            if title != None:
                yield (title,data)
            title = line[1:]
            data = ''
        else:
            data += line
    if title == None:
        yield None
    yield (title,data)

 # Generates list of locations of all substrings in seq
def substrings(seq, sub):
    start = 0
    while True:
        start = seq.find(sub, start)
        if start == -1:
            break
        else:
            yield start
            start += 1

 # Generates sequence split by n units
def split_by_n(seq, n):
    while seq:
        yield seq[:n]
        seq = seq[n:]

 # Returns a codon dict, used in to_codon
def codon_dict(): return file_to_dict("codons.txt")

 # Creates dicts from two-cell tables in txt files
def file_to_dict(file):
    mDict = {}
    with open(file) as f:
        lines = f.readlines()
        for l in lines:
            mList = l.split()
            mDict[mList[0]] = mList[1]
    return mDict

######
# START I/O METHODS
######

 # Returns a single tuple (head, name) for a given fasta file
 # In other words, gets the first fasta entry from a fasta file
def fasta_to_strings(fileName):
    with open(fileName) as f:
        for name, seq in read_fasta(f):
            return (name[1:], seq)  

def main():
        # 
    is_DNA = input("Are you entering DNA or RNA? Type 'DNA' or 'RNA': \n")
    large_fasta = input("Input scaffold's fasta filename: \n")
    small_fasta = input("Input search sequence's filename: \n")
    scaffold = fasta_to_strings(large_fasta)
    search = fasta_to_strings(small_fasta)
    
        # 
    if is_DNA == "DNA":
        search_seq = transcribe_DNA(search[1].upper())
        scaffold_seq = transcribe_DNA(scaffold[1].upper())
    elif is_DNA == "RNA":
        search_seq = search[1].upper()
        scaffold_seq = scaffold[1].upper()
    start_codons = tryFindStartCodon(scaffold_seq, search_seq)

        # If only ONE result found, fasta outputted automatically
    if len(start_codons) == 2:
        start_index = start_codons[1][0]
        length = start_codons[0]
        print("Matching sequence of length " + str(length) +
              " found in scaffold at index " + str(start_index+1) + ".")
        fileName = input("What would you like to name the output fasta file? \n")
        bp = input("How many base pairs upstream would you like to find? \n")
        scaffold_seq = RNA_to_DNA(scaffold_seq)
        buildFastaFile([(bp + " base pairs upstream of " + search[0]
                        , scaffold_seq[start_index - int(bp):start_index])], fileName)
        print("File " + fileName + " saved.")
    #elif len(start_codons) > 2:           # If multiple results found, user chooses which fasta file to output
        #i = 0
        #for start, end in 
    

main()
