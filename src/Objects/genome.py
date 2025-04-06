# -*- coding: utf-8 -*-

import random
import json

class Genome:
    
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
    readinBS(list) - List of binding sites to insert into the genome upon initialization. Used in genomeSearch with pregenerated motifs.
    Note: Length of binding sites and number of binding sites should match passed parameters lengthBS and numBS
        
    importJSON (str/Dict/None) - Used to initialize object from json file path, python dict object.
        Set to None to initialize normally. Other variables become unnecessary to specify
    Note: The length of the binding sites are correctly assumed to be the length of the sites in the motif
    '''
    def __init__(self, genomeSize = None, motif_numberBS = None, lengthBS = None, numBS = None, readinBS = None, importJSON = None):

        if importJSON is None:
            #Initializes attributes
            self.genomeSize = genomeSize
            self.motif_numberBS = motif_numberBS
            self.lengthBS = lengthBS
            self.numBS = numBS
            self.indicesBS = []
                           
            #Computes length of gene and fills with random letters
            self.geneLength = self.motif_numberBS * self.lengthBS
                
            #Throws Value error if genomeSize isn't big enough to contain gene and uniformly spaced binding sites
            if self.genomeSize < (self.geneLength + (2 * lengthBS * numBS + lengthBS)):
                raise ValueError("Genome is too small for required elements! Increase genome size")
            
            #Initializes random genome if BS aren't passed
            if readinBS is None:
                self.seq = ''.join(random.choices(["A","C","G","T"],k = self.genomeSize))
                
            #Otherwise fills genome sequence with random letters + binding sites
            else:
                self.seq = ''.join(random.choices(["A","C","G","T"],k = self.geneLength + self.lengthBS))
                
                #Adds binding sites and spacer in between each binding site
                self.seq += ''.join([bindingsite + ''.join(random.choices(["A","C","G","T"],k = self.lengthBS))
                                     for bindingsite in readinBS])
                
                #Fills in rest of genome
                self.seq += ''.join(random.choices(["A","C","G","T"],k = self.genomeSize - len(self.seq)))
            
            #Stores gene sequence and recognizer object
            self.geneSeq = self.seq[:self.geneLength]
            
            #First position for binding sites
            firstPositionBS = self.geneLength + lengthBS 
            
            #Populates index List for binding sites (binding sites are followed by spaces)
            for i in range(0, (2 * lengthBS * numBS), (2 * lengthBS)):
                self.indicesBS.append(i + firstPositionBS)
                
            #Populates list containing all index positions for BS
            self.allIndicesBS = []
            for i in self.indicesBS:
                self.allIndicesBS.append(i)
                
                #Adds other position indices
                for k in range(1,self.lengthBS):
                    self.allIndicesBS.append(i+k)
                
             
                    
        else:
            #Initializes from importJSON
            if type(importJSON) is str:
                self._import(filepath = importJSON)
            
            #Initializes from python object
            else:
                self._import(filepath = None, importObject = importJSON)


        
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
    """
    def induceMutation(self, mutateGenome = False,mutateRecognizer = True):
        '''
        Function that introduces 1 mutation in the genome at a random index with a random letter.

        In the coevolutionary simulation, this function will facilitate coevolution
        of binding sites and the motif to improve TF classification and distinction of binding sites from the rest of the genome.
        In doing this, binding sites will gain Information.
        
        In addition, in the coevolutionary simulation,this same method is used
        in experiments to move Genome organisms from MI strategies to IC strategies
        under different fitness criteria. 
        
        Parameters:
            mutateGenome: Boolean
                Set to false to not mutate the background genome (only BS and recognizer if other setting is True)
                Set to true to mutate everything if mutateRecognizer is True or everything
                everything except the recognizer
                
            mutateRecognizer: Boolean
                Set to true by default to mutate the whole genome (including the recognizer).
                Setting this parameter to false will mutate everything other than the recognizer (rest of the genome)
                
        Returns
        -------
        None.

        '''      
        #Stores number of mutations
        num_mut = 1
               
        #Indicates cutoff value for each base pair
        g_freq = 0.25
        c_freq = 0.5
        a_freq = 0.75 
        
        #Includes rest of genome for mutatable positions
        if mutateGenome:
            if mutateRecognizer:
                collection_of_indicies = random.choices(range(self.genomeSize),k=num_mut)
            else:
                collection_of_indicies = random.choices(range(self.geneLength,self.genomeSize),k=num_mut)
        else:
            if mutateRecognizer:
                #Picks viable location to be mutated
                num_positions_recognizer_BS = self.geneLength + self.lengthBS*self.numBS
                mutatePosition = random.choices(range(num_positions_recognizer_BS),k=num_mut)[0]
            else:
                mutatePosition = random.choices(range(self.lengthBS*self.numBS),k=num_mut)[0] + self.geneLength
                
            #Assigns to correct index
            if 0 <= mutatePosition <= (self.geneLength-1): #Recognizer
                collection_of_indicies = [mutatePosition]
            else: #Binding Site
                positionBS = mutatePosition-self.geneLength
                numBS = positionBS//self.lengthBS
                spotBS = positionBS %self.lengthBS
                
                collection_of_indicies = [self.indicesBS[numBS] + spotBS]
            
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
        """

    def induceMutation(self, mutateGenome = True, frozenRecognizer = False):
        """
        Function that introduces 1 mutation in the genome at a random index with a random letter.

        In the coevolutionary simulation, this function will facilitate coevolution
        of binding sites and the motif to improve TF classification and distinction of binding sites from the rest of the genome.
        In doing this, binding sites will gain Information.
        
        In addition, in the coevolutionary simulation,this same method is used
        in experiments to move Genome organisms from MI strategies to IC strategies
        under different fitness criteria. 
        
        Parameters:
            mutateGenome: Boolean
                Set to false to not mutate the background genome (only BS and recognizer if other setting is True)
                Set to true to mutate everything 
                
            frozenRecognizer: Booelean
                Set to false to mutate the recognizer. Set to true to prevent 
                mutations to the recognizer
                
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
        
        #Finds index for mutating
        if mutateGenome:
            collection_of_indicies = random.choices(range(self.genomeSize),k=num_mut)
            
        else:
            mutateIndex = random.choices(range(self.genomeSize),k=num_mut)[0]
            
            #Checks if either in gene range or in BS indices
            if (0 <= mutateIndex < self.geneLength) | (mutateIndex in self.allIndicesBS):
                collection_of_indicies = [mutateIndex] #Mutates recognizer/BS
            else:
                collection_of_indicies = [] #No indices are mutated
            
        #Removes index if in recognizer if frozenRecognizer is True
        if frozenRecognizer:
            if 0 <=collection_of_indicies[0] < self.geneLength:
                collection_of_indicies = []
                
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
    
    def _import(self,filepath = None, importObject = None):
        """
        Imports object from json file at filepath or from python object

        Parameters
        ----------
        filepath : str or None
            Path to json file for importing. If none, must specify python object in importObject
        importObject: Dict or None
            By default, none if filepath is specified. If filepath is not specified, should be set to python Dict
        Returns
        -------
        None.

        """
        #Sets import dict if reading from json or python object
        import_dic = {}
        
        if filepath is not None:
            with open(filepath, 'r') as f:
                import_dic = json.load(f)
        else:
            
            #Checks that importObject was passed
            if importObject is not None:
                import_dic = importObject
            else:
                raise ValueError("Since file path is not specified for coevGenome object, please specify importObject")
            
        #Imports attributes
        self.genomeSize = import_dic['genomeSize']
        self.motif_numberBS = import_dic['motif_numberBS']
        self.lengthBS = import_dic['lengthBS']
        self.indicesBS = import_dic['indicesBS']
        self.numBS = import_dic['numberBS']
        
        self.geneLength = import_dic['geneLength']
        self.seq = import_dic['genomeSeq']
        self.geneSeq = self.seq[:self.geneLength]
        
        #Populates list containing all index positions for BS
        self.allIndicesBS = []
        for i in self.indicesBS:
            self.allIndicesBS.append(i)
            
            #Adds other position indices
            for k in range(1,self.lengthBS):
                self.allIndicesBS.append(i+k)
        
    def _export(self,filepath = None):
        """
        Exports object to json file at filepath or to python Object

        Parameters
        ----------
        filepath : str or None
            Path to json file for exporting. If none, returns python object
            
        Returns
        -------
        None or export function
        """
        
        out_dict = {'genomeSize': self.genomeSize,
                    'motif_numberBS': self.motif_numberBS,
                    'lengthBS': self.lengthBS,
                    'numberBS': self.numBS,
                    'indicesBS': self.indicesBS,
                    'geneLength': self.geneLength,
                    'genomeSeq': self.seq,
#                    'geneSeq': self.geneSeq
                    }
        
        #Outputs to file
        if filepath is not None:
            with open(filepath, 'w') as f:
                json.dump(out_dict, f)
                
        #Otherwise returns object
        else:
            return out_dict
    
def scoreGenome(diPWM, motifLength, genome):
    """
    Scoring function that accepts diPWM and Genome object and returns a list of scores for each position
    """
    scoreList = []
    
    for i in range(len(genome) - motifLength + 1):
        scoreList.append(diPWM.set_score_from_str(genome[i:(i + motifLength)]))
    
    return scoreList

        