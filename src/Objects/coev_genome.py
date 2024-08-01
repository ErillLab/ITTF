# -*- coding: utf-8 -*-

import random
       
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
        self.geneLength = self.motif_numberBS * self.lengthBS
            
        #Throws Value error if genomeSize isn't big enough to contain gene and uniformly spaced binding sites
        if self.genomeSize < (self.geneLength + (2 * lengthBS * numBS + lengthBS)):
            raise ValueError("Genome is too small for required elements! Increase genome size")
        
        #Initializes genome
        self.seq = ''.join(random.choices(["A","C","G","T"],k = self.genomeSize))
        
        #Stores gene sequence and recognizer object
        self.geneSeq = self.seq[:self.geneLength]
        
        #First position for binding sites
        firstPositionBS = self.geneLength + lengthBS 
        
        #Populates index List for binding sites (binding sites are followed by spaces)
        for i in range(0, (2 * lengthBS * numBS), (2 * lengthBS)):
            self.indicesBS.append(i + firstPositionBS)


        
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
        numRows = self.motif_numberBS
        numCols= self.lengthBS

        #Initializes list to be populated
        motiflist = []
        
        #Iterates through rows and columns to transcribe gene
        for i in range(numRows):
            #Initializes list for binding site string
            BS = []     
            
            for j in range(i, numRows * numCols, numRows):
                #Initializes empty binding site list and populates
                BS.append(self.geneSeq[j])
            
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
        on circular genome rather than the attribute self.seq

        Returns
        circular genome (str)

        """
       
        #Returns round genome (genome + first motifLength-1 letters)
        return self.seq + self.seq[:self.lengthBS-1]
    
    def induceMutation(self, mutateRecognizer = True):
        """
        Function that introduces 1 mutation in the genome at a random index with a random letter.

        In the coevolutionary simulation, this function will facilitate coevolution
        of binding sites and the motif to improve TF classification and distinction of binding sites from the rest of the genome.
        In doing this, binding sites will gain Information.
        
        In addition, in the coevolutionary simulation,this same method is used
        in experiments to move Genome organisms from MI strategies to IC strategies
        under different fitness criteria. 
        
        Parameters:
            mutateRecognizer: Boolean
                Set to true by default to mutate the whole genome (including the recognizer).
                Setting this parameter to false will mutate everything other than the recognizer (rest of the genome)
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
        if mutateRecognizer:
            collection_of_indicies = random.choices(range(self.genomeSize),k=num_mut)
        else:
            collection_of_indicies = random.choices(range(self.geneLength,self.genomeSize),k=num_mut)
            
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
        
        #Updates genome and gene
        self.seq = "".join(genomeList)
        self.geneSeq = self.seq[:self.geneLength]

        
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
        