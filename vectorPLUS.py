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
        output += (">" + data[0])
        output += textwrap.fill(data[1], 70)
        output += "\n"
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
    with open(fileName, encoding="cp1252") as f:
        for name, seq in read_fasta(f):
            return (name, seq)


 # Called when search returns only a single result. Uses filename
 # and desired sequence length from user and outputs fasta file with
 # region upstream of found sequence
def single_result(length, start, searchName, scaffoldSeq, fileName, bp):
    print("Matching sequence of length " + str(length) +
          " found in scaffold at index " + str(start+1) + ".")
    x = start - int(bp)     # Makes sure the sequence doesn't wrap around the scaffold
    if x < 0:
        x = 0
    seq = scaffoldSeq[x:start]
    title = "Upstream of " + searchName + " from " + str(x) + " to " + str(start)
    title = title.replace("\n", "") + "\n"
    write_fasta([(title, seq)], fileName)
    print("File '" + fileName + "' saved.")


 # Called when sequenceSearch returns multiple results. Uses filename
 # and desired sequence length from user and outputs fasta file with 
 # *all* seq upstream of found sequence
def multiple_result(length, startCodons, searchName, scaffoldSeq, fileName, bp):
    print("Multiple matching sequences found of length " + str(length))
    listOfTuples = []
    for i in range(1, len(startCodons)):
        start = startCodons[i][0]
        x = start - int(bp)     # Makes sure the sequence doesn't wrap around the scaffold
        if x < 0:
            x = 0
        title = "Upstream of " + searchName + " from " + str(x) + " to " + str(start)
        title = title.replace("\n", "") + "\n"
        seq = scaffoldSeq[x:start]
        listOfTuples.append((title, seq))
    write_fasta(listOfTuples, fileName)
    print("File '" + fileName + "' saved.")


 # Searches in the scaffold for a given sequence. If it finds a sequence, allows
 # the user to output a sequence upstream of the start position in the scaffold
 # as a file. 
def sequence_search():

    print("Sequence search.")
    print()
    
        # Gets search and scaffold seq, adds .txt to names
    scaffoldFasta = input("Input scaffold's fasta filename: \n")
    searchFasta = input("Input search sequence's filename: \n")
    if scaffoldFasta.find(".txt",0) == -1: scaffoldFasta += ".txt"
    if searchFasta.find(".txt",0) == -1: searchFasta += ".txt"

        # Builds a tuple of (name, seq) out of each fasta w/ uppercase seq
    scaffold = fasta_to_strings(scaffoldFasta)
    search = fasta_to_strings(searchFasta)
    searchSeq = search[1].upper()
    scaffoldSeq = RNA_to_DNA(scaffold[1].upper())

        # Runs search
    startCodons = find_matching_seqs(scaffoldSeq, searchSeq)

        # Gets desired filename and length of upstream region to be returned
    fileName = input("What would you like to name the output fasta file? \n")
    bp = input("How many base pairs upstream would you like to find? \n")

        # Sets clear variables for search sequence's name and length
    searchName = search[0]
    length = startCodons[0]

        # If only one result found, output in fasta file
    if len(startCodons) == 2:
        start = startCodons[1][0]
        single_result(length, start, searchName, scaffoldSeq, fileName, bp)

        # If multiple results found, all are output in a single fasta file
    elif len(startCodons) > 2:
        multiple_result(length, startCodons, searchName, scaffoldSeq, fileName, bp)

        # In other cases, return no results
    else:
        print("No results found.")


    # Called in positionSearch to get correct list of positions from user.
    # If position is greater than given max integer or less than zero,
    # called recursively to get new list of integers. 
def get_positions(max): 
    inputPositions = input()
    positions = [int(x) for x in inputPositions.split(',')]
    for p in positions:
        if ((p < 0) or (p > max)):
            print("Position " + str(p) + " not in range of scaffold. Enter new position(s):")
            getPositions(max)
    return positions


 # Called in position_search. Uses filename
 # and desired sequence length from user and outputs fasta file with *all* seq X 
 # bp upstream of found sequence
def position_result(positions, scaffoldName, scaffoldSeq, fileName, bp):
    listOfTuples = []
    for i in range(0, len(positions)):
        start = positions[i]
        x = start - int(bp)     # Makes sure the sequence doesn't wrap around the scaffold
        if x < 0:
            x = 0
        title = scaffoldName + " from bp " + str(x) + " to " + str(start)
        title = title.replace("\n", "") + "\n"
        seq = scaffoldSeq[x:start]
        listOfTuples.append((title, seq))
    write_fasta(listOfTuples, fileName)
    print("File '" + fileName + "' saved.")


    # Searches in the scaffold for given positions. Allows the user
    # to output a sequence upstream of the positions as a file.
def position_search():

    print("Position search.")
    print()

        # Gets scaffold seq and adds .txt to name
    scaffoldFasta = input("Input scaffold's fasta filename: \n")
    if scaffoldFasta.find(".txt",0) == -1: scaffoldFasta += ".txt"

        # Builds (name, seq) tuple of scaffold and changes seq to uppercase
    scaffold = fasta_to_strings(scaffoldFasta)
    scaffoldName = scaffold[0]
    scaffoldSeq = RNA_to_DNA(scaffold[1].upper())

        # Gets scaffold position to find the upstream region of
    print("Enter position(s) to find, separated by commas: ")    
    positions = get_positions(len(scaffoldSeq))

        # Gets desired filename and length of upstream region to be returned
    fileName = input("What would you like to name the output fasta file? \n")
    bp = input("How many base pairs upstream would you like to find? \n")
        
        # Outputs all desired sequences in a single fasta file
    position_result(positions, scaffoldName, scaffoldSeq, fileName, bp)
    

def main():

        # Prints a nicely formatted introduction
    introText = '''Searches inside a scaffold DNA/RNA sequence for either positions
    in the scaffold or the start position of a given sequence in the
    scaffold. Outputs a region upstream of the found position(s) of a
    given length as a FASTA file.'''
    intro = textwrap.wrap(introText)
    for line in intro:
        print(line)

    print()
    print("Are you searching by position or with a sequence?")
    print("1: Position     2: Sequence")
    
    searchType = input("")
    if (searchType=="1") or (searchType=="Position") or (searchType=="position") or (searchType=="pos"):
        position_search()
    elif (searchType=="2") or (searchType=="Sequence") or (searchType=="sequence") or (searchType=="seq"):
        sequence_search()
        
    

main()
