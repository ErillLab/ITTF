# -*- coding: utf-8 -*-

from Objects.organism import Organism
import random
import copy

"""
Population class for ITTF Genetic Algorithm - Finish Documentation of Script
    Information in Transcription Factor Binding Sites (TFBS) is typically conceived
    as being conveyed column-wise, measured through Information Content (IC).
    Information, however, can also be conveyed pairwise or between columns. Mutual 
    Information measures the information shared between two variables or columns. 
    
    The IT_TF project seeks to explore questions related to evolutionary preferences
    between these two ways of storing information and what this reveals about MI.
    This algorithm attempts to provide insight into this question by 
    modeling evolution in TFBS, fitting populations of organisms 
    based on Total Information (The sum of IC and MI).
TODO
Work on the following :
    Documentation that


"""

class Population():
    """
    Population Class for IT_TF Genetic Algorithm
       
    The Population class acts as the backbone for the Genetic Algorithm, using the 
    Organism class to simulate evolution of TFBS. This class accomplishes this through the following:
        
        1. Initializes a population list of organisms of size num_organisms with passed
        organism parameters OR passed fasta files
        2. Uses a reproduce method to fill up a list of parent organisms
        and mutated offspring organisms 
        3. Uses a select method to compare parents with their respective offspring,
        adding the more fit organism to the population list for the next generation 
        
    Note: Any of the parameters marked _org are not used if read_in_org_list is populated
    with fasta file names. 
        
    Parameters:
        num_organisms - Int - Number of organisms
        
        N_org - Number of binding sites
        
        L_org - length of binding sites
        
        GC_cont_org - GC content of the organism's genome 
        
        mut_rate_org - Probability of a base pair mutating
            Note: Currently used as the probability that a base could change into any base
            in the organism's alphabet, meaning 'mutations' that do not change 
            the base are possible. 
        
        Goal_IC - Goal IC (bits)
        
        Goal_MI - Goal MI (bits)
        
        constrictMode - 1 of 3 options (["MI","IC", ""]). Parameter for the fitness function
        of the genetic algorithm. Either selects exclusively for MI, IC, or both ("IC"). 
        
        MI_mode - Choice between computing MI based on adjacent pairs (MI_mode = 2) or overlapping pairs (MI_mode = 1) 
        
        readin - Creates organism object with a list of TFBS read from a fasta file,
        overwriting N and L parameters
        
    """
    
    def __init__(self, num_organisms, N_org, L_org, GC_cont_org, mut_rate_org, Goal_IC, Goal_MI, constrictMode, Pseudocount_org = None, MI_mode_org = 2, read_in_org_list = []):
        
        #Initializes population list and scrap list
        self.population_list = []
        self.scrap_list = []
        
        #Initializes MSE list
        self.MSE_list = []

        #Initializes variables holding goal IC and MI of population
        self.Goal_IC = Goal_IC
        self.Goal_MI = Goal_MI
        
        #Initializes variables referring to distance between beginning of pairs
        self.MI_mode_org = MI_mode_org
        
        self.constrictMode = constrictMode
        
        #Sets up population list using passed Organism Class variables if no fasta list is provided
        if read_in_org_list == []:
            
            #Appends initialized organism using passed values
            for i in range(num_organisms):
                self.population_list.append(Organism(N_org,L_org,GC_cont_org,mut_rate_org,i,Pseudocount_org,MI_mode_org))
                
        else:
            #Initializes fasta list
            self.fasta_list = read_in_org_list
            count_list = []
                       
            #Creates counts list (number of organisms dedicated to each motif)
            counts = num_organisms//len(self.fasta_list)
            remainder = num_organisms%len(self.fasta_list)
            print(counts,remainder)
            for f in range(len(self.fasta_list)):
                count_list.append(counts)
            for i in range(remainder):
                count_list[i] += 1
            
            #Creates num organisms from fasta file using count_list 
            for i in range(len(count_list)):
                
                #Number of organisms to create for each fasta file
                for c in range(count_list[i]):
                    self.population_list.append(Organism(N_org,L_org,GC_cont_org,mut_rate_org,i,Pseudocount_org,MI_mode_org,self.fasta_list[i]))

                    
                
            
            
            
                
        
    
    def Reproduce(self):
        '''
        Fills up scrap list serially with parents and mutated offspring
        
        
        '''
        
        #Iterates through population list (parent generation) to create scrap list
        for i in range(len(self.population_list)):
            
            #Grabs parent organism object
            parent_organism = self.population_list[i]
            
            #Adds parent organism to scrap list
            self.scrap_list.append(parent_organism)
            
            #Creates offspring object and mutates
            offspring_organism = copy.deepcopy(parent_organism)
            offspring_organism.Mutate()
            
            #Adds offspring organism to scrap list as a copy + mutated ver of the parent
            self.scrap_list.append(offspring_organism)
            
            
    def DetermineFitness(self,IC,MI):
        '''
        Fitness function. Used on each organism of the population to compute the 
        fitness of the organism depending on self.constrictMode
        
        After reproduction, DetermineFitness is run on each organism to determine 
        the more fit organism for continuing on to the next generation. Each 
        organism's IC and MI is passed to this function and
        the mean square error from the population's Goal IC/MI is returned

        '''
        
        #Depending on self.constrictMode, returns fitness value
        if (self.constrictMode == ""):
            return ((((self.Goal_IC-IC)**(2)) + ((self.Goal_MI - MI)**2))**(1/2))
        
        #When only restricting for IC
        elif (self.constrictMode == "IC"):
            return (((self.Goal_IC-IC)**(2))**(1/2))
        
        #When only restricting for MI
        else:
            return (((self.Goal_MI-MI)**(2))**(1/2))
    
    def Select(self):
        '''
        Run after reproduction. Bottlenecks population back to num_organisms through selection
    
        Reproduction leads to one parent having one mutated offspring. This
        method selects the more fit organism (based on the method DetermineFitness)
        for proceeding to the next generation.
        
        Since a parent and their offspring were serially added to the 
        scrap list, the fitness of adjacent elements in scrap list are compared
        and the more fit organism is added back to the population list.
        '''
                
        #Zeros population list
        self.population_list = []
        #Zeros MSE list
        self.MSE_list = []
        
        #Iterates through scrap list, comparing fitness of consecutive pairs
        for i in range(0,len(self.scrap_list)-1,2):
            
            #Grabs fitness value (TI) of consecutive organisms
            fitness_organism1 = self.DetermineFitness(self.scrap_list[i].IC, self.scrap_list[i].MI)
            fitness_organism2 = self.DetermineFitness(self.scrap_list[i+1].IC, self.scrap_list[i+1].MI)
            
            #Adds organism with higher TI back to population list
            if fitness_organism1 < fitness_organism2:
                self.population_list.append(self.scrap_list[i])
                self.MSE_list.append(fitness_organism1)
                
            elif fitness_organism2 < fitness_organism1:
                self.population_list.append(self.scrap_list[i+1])
                self.MSE_list.append(fitness_organism2)
            
            #Accounts for case where TI values are equal by randomly selecting one organism
            else:
                rndm_num = random.choice([0,1])
                rndm_index = [i,i+1][rndm_num]
                self.population_list.append(self.scrap_list[rndm_index])
                self.MSE_list.append([fitness_organism1, fitness_organism2][rndm_num])

        #Resets scrap_list
        self.scrap_list = []
    
    def Update_PopList(self):
        '''
        Function to verify diversity is not lost

        Returns
        -------
        None.

        '''
        cur_organism_list = []
        self.Unique_Org = 0
        
        #Iterates through population list looking for unique organisms
        for i in range(len(self.population_list)):

            #Unique organism ID
            if not(self.population_list[i].ID in cur_organism_list):
                cur_organism_list.append(self.population_list[i].ID)
                self.Unique_Org += 1
        print(self.Unique_Org)
        self.Unique_Org = 0
        
    
    def Lose_Information(self):
        '''
        Function used to degenerate the information present in the motifs.
        Used in the Genetic Algorithm to explore how fast certain strategies
        lose information (such as in Figure 1)
        
        Simply iterates through the population list and induces each organism's 
        mutate method. Mutated population is kept in self.population_list


        '''
        #Iterates through population, mutating organism line
        for i in range(len(self.population_list)):
            
            #Creates offspring object and mutates
            self.population_list[i].Mutate()
            
    def Lose_Information2(self):
        '''
        Function used to degenerate the information present in the motifs.
        Used in the Genetic Algorithm to explore how fast certain strategies
        lose information (such as in Figure 1)
        
        Ignores mutation rate and simply induces 1 letter mutations in each organism.

        '''
        #Grabs object parameters
        lenBSorganism = self.population_list[0].length_bs
        numBSorganism = self.population_list[0].number_bs
        
        numPositions = numBSorganism * lenBSorganism
                                            
        #Iterates through population, mutating organism line
        for i in range(len(self.population_list)):
            
            #Mutates organism (Redefines organism mutation by inducing 1 mutation in organism)     
            
            #Randomly chooses position
            randPosition = random.randrange(0,(numPositions) - 1)
    
            #Gets row and column
            row = randPosition // lenBSorganism
            col= randPosition % lenBSorganism
            
            #Picks random base pair
            rand_float = random.random()
            if rand_float < 0.25: 
                #Changes to G
                self.population_list[i].TFBS[row] = self.population_list[i].TFBS[row][:col] + 'G' + self.population_list[i].TFBS[row][(col+1):]
                
            elif rand_float < 0.5:
                #Changes to C
                self.population_list[i].TFBS[row] = self.population_list[i].TFBS[row][:col] + 'C' + self.population_list[i].TFBS[row][(col+1):]
                
            elif rand_float < 0.75:
                #Changes to A
                self.population_list[i].TFBS[row] = self.population_list[i].TFBS[row][:col] + 'A' + self.population_list[i].TFBS[row][(col+1):]
    
            else:
                #Changes to T
                self.population_list[i].TFBS[row] = self.population_list[i].TFBS[row][:col] + 'T' + self.population_list[i].TFBS[row][(col+1):]
            
            self.population_list[i].Compute_Information()
            
if __name__ == "__main__":
    '''
    Run this script directly if you want to generate your own test population of motif organisms
    '''
    #Sample parameters
    Population = Population(num_organisms = 10, N_org = 32,L_org = 4, GC_cont_org = .50, mut_rate_org = 0.01,Goal_IC = 6, Goal_MI = 0, constrictMode = "IC", MI_mode_org = 1)
    
    #Boolean variable if the user wants just one organism to reach goal IC/MI
    oneOrganism = True
    
    #How close the motifs have to be to goal IC/MI
    euclideanDistance = 0.05
    
    #Keeps running until either one organism or whole population is close to goal information
    notFinished = True
    while notFinished:
        
        #Population goes through one generation of algorithm
        Population.Reproduce()
        Population.Select()
        
        #Checks if finished
        for i in range(len(Population.population_list)):
            if Population.MSE_list[i] < euclideanDistance:
                
                #Done if one organism is close to goal
                notFinished = False

            else:
                #If one organism isn't close and oneOrganism == False, continues algorithm
                if oneOrganism == False:
                    notFinished = True
            
                
                    
                    
                
        



