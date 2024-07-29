# -*- coding: utf-8 -*-

import random
from src.Objects.organism import MotifModel
class geneObject:
    """
    The motif gene is a 1D string of the motif's flattened columns (column 1 leters + column 2 letters + ...).
    By "transcribing" the gene, 
    information can be recovered from the motif.
    This algorithm setup permits mutation of the genome and simple updating
    of the information contained within the motif
    
    This object stores the motif gene object initialized in coevGenome
    Motif gene object. Initialized in coevGenome to keep track of the 
    recognizer. Contains methods for initializing the gene, transcribing the gene, 
    and computing information

    """
    def __init__(self, geneString, numBS, lengthBS, MI_mode = 1):
        """
        Initializes gene from geneString passed from coevGenome. Transcribes 
        the gene and stores as organism object. Computes information

        Parameters
        ----------
        geneString : str
            Motif gene string

        """
        #Saves gene as attribute
        self.gene = geneString
        self.numBS = numBS
        self.lenBS = lengthBS
        
        #Gets motif and saves as attribute
        self.recognizerObject = MotifModel(self.numBS, self.lenBS,0.5,0,"Deprecated", MI_mode=1, read_in = self.transcribeGene())
        
    def transcribeGene(self):
        """
        The gene is the flattened 1D version of the motif.
        To recover the information encoded in the motif, 
        this "flattening" operation is reversed and columns/rows are reconstructed 

        Returns
        -------
        motiflist: list
           List of TFBS 

        """
        #Grabs length of row/column for motif
        numRows = self.numBS
        numCols= self.lenBS

        #Initializes list to be populated
        motiflist = []
        
        #Iterates through rows and columns to transcribe gene
        for i in range(numRows):
            #Initializes list for binding site string
            BS = []     
            
            for j in range(i, numRows * numCols, numRows):
                #Initializes empty binding site list and populates
                BS.append(self.gene[j])
            
            #Creates Binding site string and adds to list
            BS = ''.join(BS)
            motiflist.append(BS)
        
        return motiflist
        
    def updateGene(self,newgene):
        """
        Updates gene and recognizer after mutation in coevGenome

        Parameters
        ----------
        newgene : str
            gene string

        """
        self.gene = newgene
        self.recognizerObject = Organism(self.numBS, self.lenBS,0.5,0,"Deprecated", MI_mode=1, read_in = self.transcribeGene())
        
    def returnMotifInformation(self):
        """
        Returns Motif Information (first element, second element, third element) = 
        (IC, MI, TI)

        Returns
        """
        return (self.recognizerObject.IC, self.recognizerObject.MI, self.recognizerObject.TI)
        
class coevGenome:
    
    '''
    Defines Genome class initialized a genome string with a motif-encoding gene 
    and indices of uniformly spaced binding sites.
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
                       
        #Computes length of gene and fills with random letters
        geneLength = self.motif_numberBS * self.lengthBS
        
        #Throws Value error if genomeSize isn't big enough to contain gene and uniformly spaced binding sites
        if self.genomeSize < (geneLength + (2 * lengthBS * numBS + lengthBS)):
            raise ValueError("Genome is too small for required elements! Increase genome size")
        
        #Initializes genome
        self.seq = ''.join(random.choices(["A","C","G","T"],k = self.genomeSize))

        #First position for binding sites
        firstPositionBS = geneLength + lengthBS 
        
        #Populates index List for binding sites (binding sites are followed by spaces)
        for i in range(0, (2 * lengthBS * numBS), (2 * lengthBS)):
            self.indicesBS.append(i + firstPositionBS)

        #Initializes gene object
        self.geneObj = geneObject(self.seq[:geneLength], self.motif_numberBS, self.lengthBS)
        
    def transcribeGene(self):
        """
        Calls geneObject method to return list of binding sites for recognizer
        """
        return self.geneObj.transcribeGene()
        
                
                
    def getCircularGenome(self):
        """
        Used to generate a circular genome by wrapping the first letters 
        to the end. This becomes important when measuring the TF's ability 
        as a classifer. In order to compute AUC, the list of scored positions
        must equal the size of the genome. As such, we run the scoring function
        on circular genome rather than the attribute self.seq

        Returns
        circular genome (str)

        """
       
        #Returns round genome (genome + first motifLength-1 letters)
        return self.seq + self.seq[:self.lengthBS-1]
    
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
        #Stores number of mutations
        num_mut = 1
               
        #Indicates cutoff value for each base pair
        g_freq = 0.25
        c_freq = 0.5
        a_freq = 0.75 
               
        #Iterative loop to ensure num_mut of unique mutations 
        collection_of_indicies = random.choices(range(genomeLength),k=num_mut)
        
        genomeList = list(self.seq)
        
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
        
        self.seq = "".join(genomeList)
        self.geneObj.updateGene(self.seq[:self.motif_numberBS * self.lengthBS])
    
    def induceMutation_freezeRecognizer(self):
        """
        Function that freezes the recognizer and introduces 1 mutation in the genome 
        at some other random index. 
        
        Used in the coevolutionary simulation to investigate the fitness decay
        of certain strategies 

        """ 
        genomeLength = len(self.seq)
        
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
            genomeIndex = random.randrange(self.motif_numberBS * self.lengthBS,genomeLength)
            sufficient_mut += 1
            collection_of_indicies.append(genomeIndex)
        
        genomeList = list(self.seq)
        
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
        
        self.seq = "".join(genomeList)
        
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
            bindingSitelist.append(self.seq[i:(i+self.lengthBS)])
        
        return bindingSitelist
    
def scoreGenome(diPWM, motifLength, genome):
    """
    Scoring function that accepts diPWM and Genome object and returns a list of scores for each position
    """
    scoreList = []
    
    for i in range(len(genome) - motifLength + 1):
        scoreList.append(diPWM.set_score_from_str(genome[i:(i + motifLength)]))
    
    return scoreList
        