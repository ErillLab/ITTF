# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 15:03:32 2024

@author: thele
"""
import random
import csv
import math
import dipwmsearch as ds

def create_dipwm_list(TFBS_list):
    """
    Helper Function that accepts Transcription Factor Binding sites, computes dinucleotide
    frequencies at each position (16 potential states), and creates diNucleotide
    list appropriate for creating Kulaokovsky diPWM object
    
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
    
    #returns file
    return scoreList

def create_dipwm_object(TFBS_list):
    """
    Function used to convert list of transcription factor binding sites
    into a diPWM object as defined by the Kulakovsky package
    
    Parameters
    ----------
    scoreList : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    scoreList = create_dipwm_list(TFBS_list)
    return ds.diPWM(scoreList)
        
class coevGenome:
    
    '''
    Defines Genome class initialized with a motif-encoding gene and indices of uniformly spaced binding sites.
    The motif-encoding gene is at the beginning of the genome (index 0) and 
    is a 1D string of the motif's flattened columns (column 1 leters + column 2 letters + ...).
    As such, the gene encodes none of the information of the motif. This operation is reversed
    to transcribe the gene back into the motif for scanning and identifying binding sites.
    Used in coevolutionary simulation to generate motif embedded in genomes. 
    
    Contains methods for transcribing the gene, computing information of binding sites, and
    returning circular genome for scoring
    
    genomeSize (int) - Length of the genome. Must be sufficient for gene size and BS + spacing
    motif_numberBS (int) - Number of binding sites in the motif (# rows). Used to initialize the length of the gene
    lengthBS (int) - Length of binding sites in the motif (# columns). Used to initialize the length of the gene                                                            
    numberBS (int) - Number of binding sites in the genome. Used to initialize a list of indices
    
    Note: The length of the binding sites are correctly assumed to be the length of the sites in the motif
    '''
    def __init__(self, genomeSize, motif_numberBS, lengthBS, numBS):
        #Initializes attributes
        self.genomeSize = genomeSize
        self.motif_numberBS = motif_numberBS
        self.lengthBS = lengthBS
        self.indicesBS = []
                
        #The following lines of code create the genome (Gene + (binding sites + uniform space) + Random letters)
        genomeList = []
        
        #Computes length of gene and fills with random letters
        geneLength = self.motif_numberBS * self.lengthBS
        
        #Throws Value error if genomeSize isn't big enough to contain gene and uniformly spaced binding sites
        if self.genomeSize < (geneLength + (2 * lengthBS * numBS + lengthBS)):
            raise ValueError("Genome is too small for required elements! Increase genome size")
        
        #First geneLength letters of the genome are the gene 
        for i in range(geneLength):
            genomeList.append(random.choice(["A", "C", "G", "T"]))
        
        #Adds BS length space between end of gene and beginning of binding sites
        #Then adds letters for binding sites and binding site-sized spaces to genome
        for i in range(lengthBS + 2 * lengthBS * numBS):
            genomeList.append(random.choice(["A", "C", "G", "T"]))
        
        #Populates index list for Binding Sites
        firstPositionBS = geneLength + lengthBS #First position for binding sites
        #lastPositionBS = firstPositionBS + (2 * lengthBS * numBS)
        
        #Populates index List
        for i in range(0, (2 * lengthBS * numBS), (2 * lengthBS)):
            self.indicesBS.append(i + firstPositionBS)
        
        #Fills in the rest of the genome
        firstPositionGenome = firstPositionBS + 2 * lengthBS * numBS
        for i in range(firstPositionGenome, self.genomeSize):
            genomeList.append(random.choice(["A", "C", "G", "T"]))
            
        #Initializes genome
        self.genome = ''.join(genomeList)
        
    def transcribeGene(self):
        """
        The gene (at the beginning of the genome) is the flattened 1D version 
        of the motif so that it contains no information from the motif. To 
        obtain the motif, this "flattening" operation is reversed and columns/rows
        are reconstructed 

        Returns
        -------
        None.

        """
        #Grabs length of row/column for motif
        rowLength = self.motif_numberBS
        columnLength = self.lengthBS
        
        #Initializes list to be populated
        motiflist = []
        
        #Iterates through rows and columns to transcribe gene
        for i in range(rowLength):
            #Initializes list for binding site string
            BS = []     
            
            for j in range(i, rowLength * columnLength, rowLength):
                #Initializes empty binding site list and populates
                BS.append(self.genome[j])
            
            #Creates Binding site string and adds to list
            BS = ''.join(BS)
            motiflist.append(BS)
        
        return motiflist
                
                
    def getCircularGenome(self):
        """
        Used to generate a circular genome by wrapping the first letters 
        to the end. This becomes important when measuring the TF's ability 
        as a classifer. In order to compute AUC, the list of scored positions
        must equal the size of the genome. As such, we run the scoring function
        on circular genome rather than the attribute self.genome

        Returns
        circular genome (str)

        """
       
        #Returns round genome (genome + first motifLength-1 letters)
        firstLetters = self.genome[0:(self.lengthBS-1)]
        return self.genome + firstLetters
    
    def induceMutation(self):
        """
        Function that introduces 1 mutation in the genome at a random index with a random letter.

        In the coevolutionary simulation, this function will facilitate coevolution
        of binding sites and the motif to improve TF classification and distinction of binding sites from the rest of the genome.
        In doing this, binding sites will gain Information.
        
        In addition, in the coevolutionary simulation,this same method is used
        in experiments to move Genome organisms from MI strategies to IC strategies
        under different fitness criteria. 

        Returns
        -------
        None.

        """ 
        genomeLength = len(self.genome)
        
        #Stores number of mutations
        num_mut = 1
               
        #Indicates cutoff value for each base pair
        g_freq = 0.25
        c_freq = 0.5
        a_freq = 0.75 
        
        #Prepares to draw num_mut amount of numbers between 1 and the genome length (sites for mutation)
        collection_of_indicies = []
        sufficient_mut = 0
        
        #Iterative loop to ensure num_mut of unique mutations 
        while sufficient_mut < num_mut: 
            genomeIndex = random.randrange(0,genomeLength)
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
    
        
    def returnBindingSites(self):
        """
        Method used to return the binding sites stored in the Genome. These
        binding sites can then be fed into the "Organism" class in order to compute
        IC, MI, and TI
        
        The index of binding sites is stored as an attribute. To return the binding sites,
        this list is used in tandem with the length of the binding site
        
        Returns
        -------
        bindingSitelist (list) - List of binding sites in sequential order

        """
        #Initializes BS list
        bindingSitelist = []
        
        #Iterates through list and grabs binding sites
        for i in self.indicesBS:
            bindingSitelist.append(self.genome[i:(i+self.lengthBS)])
        
        return bindingSitelist
    
def scoreGenome(diPWM, motifLength, genome):
    """
    Scoring function that accepts diPWM and Genome object and returns a list of scores for each position
    """
    scoreList = []
    
    for i in range(len(genome) - motifLength + 1):
        scoreList.append(diPWM.set_score_from_str(genome[i:(i + motifLength)]))
    
    return scoreList
        