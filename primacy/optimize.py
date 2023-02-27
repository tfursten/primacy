import logging
import copy
import pandas as pd
import numpy as np
import itertools as it
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import Seq
from collection import calculate_hybrid_score

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('primacy')

class Primer(object):
    def __init__(self, seq, name=None, amplicon=None, flank=None):
        self.seq = seq
        self.name = name
        self.amplicon = amplicon
        self.flank = flank

class PrimerPair(object):
    def __init__(self, amplicon_table, max_sz, min_sz):
        self.max_sz = max_sz
        self.min_sz = min_sz
        self.amplicon_table = amplicon_table
        self.forward_primer, self.reverse_primer = self.get_primer_pair(
            f=None, r=None)
    
    def mutate(self, flank=None):
        if flank == None:
            flank = np.random.choice(['F', 'R'])
        if flank == "R":
            self.forward_primer, self.reverse_primer = self.get_primer_pair(
                f = self.forward_primer
            )
        else:
            self.forward_primer, self.reverse_primer = self.get_primer_pair(
                r = self.reverse_primer
            )

    def get_primer_pair(self, f=None, r=None):
        """
        f and r are primer objects
        """
        loop_iter = 0
        if f==None and r==None:
            while True:
                loop_iter += 1
                forward = self.amplicon_table[self.amplicon_table['Flank'] == 'F'].sample(1).iloc[0]
                reverse = self.amplicon_table[(self.amplicon_table['Flank'] == 'R')]
                if self.max_sz != None:
                    reverse = reverse[reverse['Position'] <= (forward['Position'] + self.max_sz)]
                if self.min_sz != None:
                    reverse = reverse[reverse['Position'] >= (forward['Position'] + self.min_sz)]
                if len(reverse):
                    reverse = reverse.sample(1).iloc[0]
                    break
                if loop_iter > 1000:
                    raise ValueError("Having trouble finding amplicons of the correct size, please adjust parameters.")
        elif f != None:
            forward = self.amplicon_table.loc[f.name]
            reverse = self.amplicon_table[(self.amplicon_table['Flank'] == 'R')]
            if self.max_sz != None:
                reverse = reverse[reverse['Position'] <= (forward['Position'] + self.max_sz)]
            if self.min_sz != None:
                reverse = reverse[reverse['Position'] >= (forward['Position'] + self.min_sz)]
            if len(reverse):
                reverse = reverse.sample(1).iloc[0]
            else:
                # Run again replacing both primers
                return self.get_primer_pair(f=None, r=None)
        elif r != None:
            reverse = self.amplicon_table.loc[r.name]
            forward = self.amplicon_table[(self.amplicon_table['Flank'] == 'F')]
            if self.max_sz != None:
                forward = forward[forward['Position'] <= (reverse['Position'] + self.max_sz)]
            if self.min_sz != None:
                forward = forward[forward['Position'] >= (reverse['Position'] + self.min_sz)]
            if len(forward):
                forward = forward.sample(1).iloc[0]
            else:
                return self.get_primer_pair(f=None, r=None)
        return [Primer(seq=forward['Seq'], name=forward.name), Primer(seq=reverse['Seq'], name=reverse.name)]



class Panel(object):
    def __init__(self, chromosome):
        self.chromosome = chromosome

    def __repr__(self):
        return "\n".join([",".join([p.forward_primer.name, p.reverse_primer.name]) for p in self.chromosome])
    
    def hybrid_score(self, required_primers=[]):
        primers = self.get_all_primer_seqs(required_primers)
        return np.max(np.array(
            [calculate_hybrid_score(Seq(primers[i]), Seq(primers[j]))
             for i in range(len(primers) - 1) for j in range(i+1, len(primers))]))
        
    def get_all_primer_seqs(self, required_primers=[]):
        primers = []
        for locus in self.chromosome:
            primers.append(locus.forward_primer.seq)
            primers.append(locus.reverse_primer.seq)
        for primer in required_primers:
            primers.append(primer.seq)
        return primers

    def mutate(self):
        np.random.choice(self.chromosome).mutate()

    def crossover(self, start, stop):
        return self.chromosome[start: stop]
    

class Population(object):
    def __init__(self, primer_options, required_primers, pop_size, iterations, mut_rate, cross_rate, max_sz, min_sz):
        self.pop_size = pop_size
        self.iter = iterations
        self.mut_rate = mut_rate
        self.cross_rate = cross_rate
        self.max_sz = max_sz
        self.min_sz = min_sz
        self.required_primers = required_primers
        self.primer_options = primer_options
        self.targets = self.primer_options['Amplicon'].unique()
        self.population = self.initialize_population()
    def initialize_population(self):
        pop = []
        for i in range(self.pop_size):
            chromosome = []
            for amp in self.targets:
                amp_table = self.primer_options[self.primer_options['Amplicon'] == amp]
                chromosome.append(PrimerPair(amp_table, self.max_sz, self.min_sz))
            pop.append(Panel(chromosome))
        return pop
    
    def mating(self):
        offspring = []
        for i in range(self.pop_size):
            mate1, mate2 = np.random.choice(self.population, size=2, replace=False)
            mate1 = Panel(copy.deepcopy(mate1.chromosome))
            mate2 = Panel(copy.deepcopy(mate2.chromosome))
            self.mutate(mate1)
            self.mutate(mate2)
            self.crossover(mate1, mate2)
            offspring.append(mate1)
            offspring.append(mate2)
        self.population += offspring


    def mutate(self, panel):
        if np.random.random() <= self.mut_rate:
            panel.mutate()
            
    def crossover(self, mate1, mate2):
        # don't do cross over if there is only one amplicon
        if len(self.targets) == 1:
            return
        if np.random.random() <= self.cross_rate:
            pos = np.random.randint(1, len(self.targets))
            mate1.chromosome = mate1.crossover(0, pos) + mate2.crossover(pos, len(self.targets))
            mate2.chromosome = mate2.crossover(0, pos) + mate1.crossover(pos, len(self.targets))
    

    def selection(self):
        fitness = [
            panel.hybrid_score(required_primers=self.required_primers)
            for panel in self.population]
        logger.info(
            "Population Min Hybrid: {0} Mean Hybrid: {1:.2} Median Hybrid: {2} Max Hybrid {3}".format(
            min(fitness), np.mean(fitness), np.median(fitness), max(fitness)))
        # get only the panels with the lowest max hybrid scores
        fitness_indx = np.argsort(fitness)[:self.pop_size]
        self.population = list(np.take(self.population, fitness_indx))

    def get_minimum(self):
        hybrid_scores = [
            panel.hybrid_score(required_primers=self.required_primers)
            for panel in self.population]
        best_idx = np.argmin(hybrid_scores)
        return self.population[best_idx]
        
    def generation(self):
        self.mating()
        self.selection()

    def run_optimization(self):
        for gen in range(self.iter):
            logger.info("Generation {}".format(gen))
            self.generation()


def get_result_table(result, primer_options):
    primers = []
    for amp in result.chromosome:
        primers.append(amp.forward_primer.name)
        primers.append(amp.reverse_primer.name)
    return primer_options.loc[primers]    

# TODO: report amplicon size and max hybrid score for each primer


def run_optimization(
    primer_options, required_primers,
    max_amp_len, min_amp_len,
    pop_size, iterations, mutation_rate, crossover_rate):

    primers = []
    for primer_set in primer_options:
        primers.append(pd.read_csv(primer_set, sep="\t", index_col='PrimerName'))
    primers = pd.concat(primers)
    if required_primers:
        req_primers = list([
            Primer(name=n, seq=r.iloc[0]) for n, r in 
            pd.read_csv(required_primers, sep="\t", header=None, index_col=0, comment="#").iterrows()])
    else:
        req_primers = []
    pop = Population(
        primers, req_primers, pop_size,
        iterations, mutation_rate, crossover_rate, max_amp_len, min_amp_len)
    pop.run_optimization()
    result = get_result_table(pop.get_minimum(), primers)
    return result

















