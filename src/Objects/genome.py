import random
import numpy as np
import math
import csv

def CreateDiPWMfile(TFBS_list, idText, filePath):
    """
    Function that accepts Transcription Factor Binding sites, computes dinucleotide
    frequencies at each position (16 potential states), and writes results to 
    a file fitting the format accepted by the diPWM library
    
    This function does this by populating a list of dictionaries with the dinucleotide counts 
    at each position of the motif, converting counts to scores using the diChipMunk PWM formula,
    and writing to a .dpwm file path provided to the function in the parameters (The format for these files 
    entail a position x dinucleotide matrix separated by tabs)
    
    diChIPMunk formula: log(x + cq/((N + c)q)) where x = counts of a dinucleotide at a position, c = pseudocount (ln(N)),
    q = background frequency of nucleotide (uniform), N = number of sequences
    Parameters
    ----------
    TFBS_list: List
        List of Transcription Factor Binding sites.
    
    idText: str
        String to be written after carrot sign on first line of file
    
    filePath: str
        Path to a .dpwm file which will be written to
        
    Returns
    -------
    None.

    """
    
    #Initializes alphabet (16 dinucleotides)
    diNucAlphabet = ["AA", "AC", "AG", "AT", "CA", "CC", 
                     "CG", "CT", "GA", "GC", "GG", "GT",
                     "TA", "TC", "TG", "TT"]
    
    #Initailizes frequency list which contains dictionaries of diNuc freqs for each position
    diNucFreqs = []
    
    #Initializes variable holding length of binding sites
    lenBS = len(TFBS_list[0])
    
    #Populates with Counts
    for i in range(lenBS-1):
        newDict = {}
        for diNuc in diNucAlphabet:
            newDict[diNuc] = 0
        diNucFreqs.append(newDict)
    
    #Fills diNucFreqs with counts at each position
    for i in range(0, lenBS-1):
        
        #Iterates through each binding site in list
        for bsNum in range(0,len(TFBS_list)):
            site = TFBS_list[bsNum]
            
            #Adds count to diNucFreq
            diNucFreqs[i][site[i] + site[i+1]] = diNucFreqs[i][site[i] + site[i + 1]] + 1
    
    #Variables related to computing score in diChipMunk manner
    c_Munk = math.log(len(TFBS_list)) #Pseudocount parameter
    q_Munk = 1/16 #Background frequency (uniform)
    N_Munk = len(TFBS_list) #Number of sequences
    
    
    #Converts to scores
    for i in range(len(diNucFreqs)):
        
        #Iterates across alphabet to compute score at each position
        for diNuc in diNucAlphabet:
            diNucFreqs[i][diNuc] = math.log(
                (diNucFreqs[i][diNuc] + (c_Munk * q_Munk))/((N_Munk + c_Munk) * q_Munk)
                                            ,2)
            
            
    
    #Begins to write to file by creating lists for the scores of each dinuc at each positon
    scoreList = []
    
    #Iterates across position
    for i in range(len(diNucFreqs)):
        #Adds potential scores at a position to a list and puts it in scoreList
        positionScore = [score for score in diNucFreqs[i].values()]
        scoreList.append(positionScore)
    
    #Writes to file
    with open(filePath, "w", newline = '') as scoreFile:
        writer = csv.writer(scoreFile, delimiter = '\t', lineterminator = '\n')
        
        #Writes first line
        writer.writerow(">")
        
        #Writes subsequent lines (scores)
        for i in range(len(scoreList)):
            writer.writerow(scoreList[i])
        scoreFile.close()


class Genome:
    
    '''
    Defines Genome class initialized from random bases + organism's binding sites.
    Initialized with genomeString and returns circular genome via induceMethod()
    
    genomeString - str: Linear genome
    motifLength - int: Length of the motif. Used to create circularGenome
    '''
    def __init__(self, nonBSGenomeLength, bindingSites):
        #The following lines of code create the genome (Random letters + binding sites)
        genomeList = []

        self.positiveNegativeList = [] #List of positive and negatives for each genome index
        self.motifLength = len(bindingSites[0]) #Grabs motif length 
        self.nonBSGenomeLength = nonBSGenomeLength #Makes class attribute
        
        #Adds random letters ()
        for i in range(nonBSGenomeLength):
            genomeList.append(random.choice(["A", "C", "G", "T"]))
            self.positiveNegativeList.append(0)
            
        #Adds binding sites to end of genomeList
        for bs in bindingSites:
            genomeList = genomeList + list(bs)
            self.positiveNegativeList.append(1) #first position is positive (1)
            for i in range(self.motifLength-1):
                self.positiveNegativeList.append(0) #Following positions are negatives (0)
                
        #Initializes genome string
        self.genome = ''.join(genomeList)
        
    
    def getCircularGenome(self):
        #Grabs motif length for getCircularGenome method
        
        
        #Returns round genome (genome + first motifLength-1 letters)
        firstLetters = self.genome[0:(self.motifLength-1)]
        return self.genome + firstLetters
        
    def induceMutation(self, mut_rate):
        """
        Function that mutates the genome based on mut_rate, randomly sampling indicies 
        and randomly selecting letters for changing the nucleotide at select indicies
        
        More generally, this function allows us to degrade the information contained 
        within the sites in the genome, which we should be able to quantify through
        AUC. This function allows us to generate AUCs for each organism after each mutational event, 
        which should decrease as the binding sites become more similar to the rest of the genome.

        Parameters
        ----------
        genome : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """ 
        genomeLength = len(self.genome)
        #Computes expected mutations
        expected_mut = mut_rate * genomeLength
       
        #Generates the # of mutations based on expected value using the poisson distribution
        num_mut = np.random.poisson(expected_mut)
        
        #Indicates cutoff value for each base pair
        g_freq = 0.25
        c_freq = 0.5
        a_freq = 0.75 
        
        #Prepares to draw num_mut amount of numbers between 1 and the genome length (sites for mutation)
        collection_of_indicies = []
        sufficient_mut = 0
        
        #Iterative loop to ensure num_mut of unique mutations 
        while sufficient_mut < num_mut: 
            genomeIndex = random.randrange(0,genomeLength - 1)
            sufficient_mut += 1
            collection_of_indicies.append(genomeIndex)
        
        genomeList = list(self.genome)
        
        #Iterates through list and finds referenced position within genome
        for index in collection_of_indicies:
            
            #Generates random float between 0 and 1 used to mutate the positions stored in the above array
            rand_float = random.random()
            
            #Applies mutation
            if rand_float < g_freq:
                genomeList[index] = "G"
            
            elif rand_float < c_freq:
                genomeList[index] = "C"
            
            elif rand_float < a_freq:
                genomeList[index] = "A"
            
            else:
                genomeList[index] = "T"
        
        self.genome = "".join(genomeList)
    
    def induceMutationNonBindingSites(self, mut_rate):
        genomeLength = len(self.genome)
        #Computes expected mutations
        expected_mut = mut_rate * (genomeLength - self.nonBSGenomeLength)
       
        #Generates the # of mutations based on expected value using the poisson distribution
        num_mut = np.random.poisson(expected_mut)
        
        #Indicates cutoff value for each base pair
        g_freq = 0.25
        c_freq = 0.5
        a_freq = 0.75 
        
        #Prepares to draw num_mut amount of numbers between 1 and the genome length (sites for mutation)
        collection_of_indicies = []
        sufficient_mut = 0
        
        #Iterative loop to ensure num_mut of unique mutations 
        while sufficient_mut < num_mut: 
            genomeIndex = random.randrange(0,(genomeLength - self.nonBSGenomeLength) - 1)
            sufficient_mut += 1
            collection_of_indicies.append(genomeIndex)
        
        genomeList = list(self.genome)
        
        #Iterates through list and finds referenced position within genome
        for index in collection_of_indicies:
            
            #Generates random float between 0 and 1 used to mutate the positions stored in the above array
            rand_float = random.random()
            
            #Applies mutation
            if rand_float < g_freq:
                genomeList[index] = "G"
            
            elif rand_float < c_freq:
                genomeList[index] = "C"
            
            elif rand_float < a_freq:
                genomeList[index] = "A"
            
            else:
                genomeList[index] = "T"
            
            if index >= genomeLength:
                raise(ValueError("Big problem with induce Mutation Non Binding Sites operator"))
        
        self.genome = "".join(genomeList)
        
    def returnMotifSequence(self):
        """
        Returns the genome sequence corresponding to the motif's binding sites

        Returns
        -------
        None.

        """
        #The first "nonBSGenomeLength" letters are not binding sites. 
        #Returns starting from the first index where binding sites are placed
        return self.genome[self.nonBSGenomeLength:]
        
def scoreGenome(diPWM, motifLength, genome):
    """
    Scoring function that accepts diPWM and Genome object and returns a list of scores for each position
    """
    scoreList = []
    
    for i in range(len(genome) - motifLength + 1):
        scoreList.append(diPWM.set_score_from_str(genome[i:(i + motifLength)]))
    
    return scoreList