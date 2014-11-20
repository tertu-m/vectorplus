import textwrap

######
# START IMPORTANT METHODS
######

 # Returns a list (see below) of all matching sequences in the given scaffold
 # with a given initial search length.
 # Returns: [search_length, (start_index, end_index), ...]
def find_matching_seqs(scaffold, sequence, initialSearchLength=5):
    searchLength = initialSearchLength - 1
    outputList = []
    workList = []
    while True:
        searchLength += 1
        outputList = workList
        workList = []
        for start in substrings(scaffold, sequence[:searchLength]):
            workList.append((start, start+searchLength))
        if len(workList)==0:
            outputList.insert(0, searchLength-1)
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

######
# START UTILITY METHODS
######

 # generator for fasta files - returns (name, seq)
 # Call using 'with open(file) as f'
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
            data += line.strip()
    if title == None:
        yield None
    yield (title,data)

 # Takes a list of [(title, seq), ...] tuples and prints as a FASTA file
 # to the given filename
def write_fasta(listOfTuples, fileName):
    output = ""
    if fileName.find(".txt",0) == -1:
        fileName += ".txt"
    for data in listOfTuples:
        output += (">" + data[0] + "\n")
        output += textwrap.fill(data[1], 80)
    with open(fileName, 'w') as outFile:
        outFile.write(output)

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

 # Called when search returns only a single result. Asks for filename
 # and base pairs from user and outputs fasta file with seq X bp
 # upstream of found sequence
def single_result(length, start, searchName, scaffoldSeq, fileName, bp):
    print("Matching sequence of length " + str(length) +
          " found in scaffold at index " + str(start+1) + ".")
    seq = scaffoldSeq[start - int(bp):start]
    write_fasta([(bp + " base pairs upstream of " + searchName, scaffoldSeq)], fileName)
    print("File '" + fileName + "' saved.")


def main():
    
        # Gets search FASTA and scaffold to be searched in from user, adds .txt to name
    largeFasta = input("Input scaffold's fasta filename: \n")
    smallFasta = input("Input search sequence's filename: \n")
    if largeFasta.find(".txt",0) == -1: largeFasta += ".txt"
    if smallFasta.find(".txt",0) == -1: smallFasta += ".txt"

        # Builds a tuple of (name, seq) out of each fasta
    scaffold = fasta_to_strings(largeFasta)
    search = fasta_to_strings(smallFasta)
    
        # Changes seq to uppercase
    searchSeq = search[1].upper()
    scaffoldSeq = RNA_to_DNA(scaffold[1].upper())

        # Runs search
    startCodons = find_matching_seqs(scaffoldSeq, searchSeq)

        # Gets desired filename and base pairs to be returned
    fileName = input("What would you like to name the output fasta file? \n")
    bp = input("How many base pairs upstream would you like to find? \n")

        # Sets clear variables for search sequence's name and length
    searchName = search[0]
    length = startCodons[0]

        # Calls single_result
    if len(startCodons) == 2:
        start = startCodons[1][0]
        single_result(length, start, searchName, scaffoldSeq, fileName, bp)

        # If multiple results found, all are output in a single fasta file
    elif len(startCodons) > 2:
        print("Multiple matching sequences found. \n")
        listOfTuples = []
        for i in range(1, len(startCodons)):
            start = startCodons[i][0]
            title = bp + " base pairs upstream of " + searchName + " at index " + start
            seq = scaffoldSeq[start - int(bp):start]
            listOfTuples.append((title, seq))
        write_fasta(listOfTuples, fileName)

        # In other cases, return no results
    else:
        print("No results found.")
        
    

main()
