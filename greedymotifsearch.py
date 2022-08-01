def ProbableKmer(string, matrix): ## matrix is a list of 4 lists
    probable = 1 ## probability of the k-mer
    for i in range(len(string)): ## for each nucleotide in the k-mer
        if string[i] == 'A': ## if the nucleotide is A
            probable *= matrix[0][i] ## multiply the probability by the probability of the nucleotide
        if string[i] == 'C': ## if the nucleotide is C
            probable *= matrix[1][i] ## multiply the probability by the probability of the nucleotide
        if string[i] == 'G': ## if the nucleotide is G
            probable *= matrix[2][i] ## multiply the probability by the probability of the nucleotide
        if string[i] == 'T': ## if the nucleotide is T
            probable *= matrix[3][i] ## multiply the probability by the probability of the nucleotide
    return probable ## return the probability of the k-mer

# Profile-most probable k-mer in the i-th string in Dna
def FindProfileMostProbableKmer(string, k, matrix):  ## matrix is a list of 4 lists
    seq = {} ## dictionary of k-mers and their probabilities
    for i in range(len(string) - k + 1): ## for each k-mer in the string
        seq[string[i:i + k]] = ProbableKmer(string[i:i + k], matrix) ## add the k-mer and its probability to the dictionary
    max_key = sorted(seq.items(), key=lambda x:x[1], reverse=True)[0][0] ## get the k-mer with the highest probability
    return max_key ## return the k-mer with the highest probability

# Score(Motifs)
def Score(Motifs): ## Motifs is a list of strings
    score = 0 ## score of the motifs
    for i in range(len(Motifs[0])): ## for each nucleotide in the motif
        j = [motif[i] for motif in Motifs] ## get the nucleotide in each motif
        score += (len(j) - max(j.count("A"), j.count("C"), j.count("T"), j.count("G")))   ## add the number of nucleotides in the motif that are not A, C, T, or G
        return score ## return the score of the motifs
 
def GreedyMotifSearch(Dna, k, t): ## Dna is a list of strings, k is the length of the k-mers, t is the number of strings in Dna
    # BestMotifs ← motif matrix formed by first k-mers in each string from Dna
    BestMotifs = [dna[:k] for dna in Dna]   ## create a list of the first k-mers in each string in Dna
    # for each k-mer Motif in the first string from Dna 
    for k_mer in [Dna[0][i:i+k] for i in range(len(Dna[0])-k+1)]: ## for each k-mer in the first string from Dna
        # Motif1 ← Motif
        Motifs = [k_mer] ## create a list of the first k-mer in the first string from Dna
        # for i = 2 to t
        for i in range(1, t):   ## for each string in Dna
            # form Profile from motifs Motif1, …, Motifi - 1
            motifs = Motifs[:i] ## create a list of the motifs Motif1, …, Motifi - 1
            # Motifi ← Profile-most probable k-mer in the i-th string in Dna
            matrix = [] ## create a list of 4 lists
            for nar in ["A", "C", "G", "T"]: ## for each nucleotide in the alphabet
                mat = [] ## create a list of the probabilities of the nucleotides
                for j in range(k): ## for each nucleotide in the k-mer
                    mm = [m[j] for m in motifs] ## create a list of the nucleotides in the motifs
                    mat.append(mm.count(nar)/len(motifs)) ## add the probability of the nucleotide to the list
                matrix.append(mat) ## add the list of the probabilities of the nucleotides to the list of 4 lists
            # Motifs ← (Motif1, …, Motift)    
            Motifs.append(FindProfileMostProbableKmer(Dna[i], k, matrix)) ## add the profile-most probable k-mer in the i-th string in Dna to the list of motifs
        # print(Motifs) 
        # if Score(Motifs) < Score(BestMotifs), BestMotifs ← Motifs
        if Score(Motifs) < Score(BestMotifs):  ## if the score of the motifs is less than the score of the best motifs
            BestMotifs = Motifs ## set the best motifs to the motifs
    return BestMotifs ## return the best motifs

if __name__ == "__main__":  ## if the program is run directly
    with open(r"D:\Python\rosalind_ba2d.txt", "r") as f: ## open the file
        k, t = map(int, f.readline().strip().split()) ## get the length of the k-mers and the number of strings in Dna
        Dna = [line.strip() for line in f] ## create a list of the strings in Dna
    # print(k,t,Dna) 
    BestMotifs = GreedyMotifSearch(Dna, k ,t) ## get the best motifs   
    print("\n".join(BestMotifs)) ## print the best motifs
