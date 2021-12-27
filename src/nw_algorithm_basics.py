from enum import Enum

class op(Enum):
    GAP_A = 1
    GAP_B = 2
    MA_MM = 3

def parseNucleotide(nucleotide):
    switcher = {
        'A': 0,
        'T': 1,
        'G': 2,
        'C': 3,
        'N': 4
    }
    return switcher.get(nucleotide,"Invalid nucleotide")

def get_score(score_mtx,characterA,characterB):
        indexA = parseNucleotide(characterA)
        indexB = parseNucleotide(characterB)
        return score_mtx[indexA,indexB]

def traceback(ops_matrix, seqA, seqB):
    seqA_alineada = []
    seqB_alineada = []
    i = len(seqA)
    j = len(seqB)
    while i > 0 or j > 0:
        if ops_matrix[i, j] == op.GAP_A.value:
            seqA_alineada.append("-")
            seqB_alineada.append(seqB[j -1])
            j -= 1
        elif ops_matrix[i, j] == op.GAP_B.value:
            seqB_alineada.append("-")
            seqA_alineada.append(seqA[i -1])
            i -= 1
        else: # es decir, match / mismatch
            i -= 1
            j -= 1
            seqA_alineada.append(seqA[i])
            seqB_alineada.append(seqB[j])
            
    seqA_alineada = seqA_alineada[::-1]
    seqB_alineada = seqB_alineada[::-1]
    
    return seqA_alineada, seqB_alineada

def init_ops_matrix(sizeA,sizeB):
    matrix = np.repeat(0, (sizeA+1) * (sizeB+1)).reshape(sizeA+1, sizeB+1)
    
    for i in range(1,sizeB+1):
        matrix[0][i] = op.GAP_A.value
    for j in range(1,sizeA+1):
        matrix[j][0] = op.GAP_B.value
    return matrix

def nw_basic(seqA_str, seqB_str, score_mtx, gap_score = 0):
    sizeA = len(seqA_str)
    sizeB = len(seqB_str)
    sol_matrix = np.repeat(0, (sizeA+1) * (sizeB+1)).reshape(sizeA+1, sizeB+1)
    ops_matrix = init_ops_matrix(sizeA,sizeB)
    
    for i in range(1,sizeA+1):
        for j in range(1,sizeB+1):
            
            diagonal_score = get_score(score_mtx,seqA_str[i-1],seqB_str[j-1])
            diagonal_result = sol_matrix[i-1][j-1] + diagonal_score
            left_result = sol_matrix[i][j-1] + gap_score
            above_result = sol_matrix[i-1][j] + gap_score
            results = [diagonal_result,left_result,above_result]
            result_index = np.argmax(results)
            if result_index == 0:
                sol_matrix[i][j] = diagonal_result
                ops_matrix[i][j] = op.MA_MM.value
            elif result_index == 1:
                sol_matrix[i][j] = left_result
                ops_matrix[i][j] = op.GAP_A.value
            else:
                sol_matrix[i][j] = above_result
                ops_matrix[i][j] = op.GAP_B.value

    aln_a, aln_b = traceback(ops_matrix,seqA_str,seqB_str)

    profile = Profile()
    profile.add_seq(aln_a,sol_matrix[sizeA, sizeB])
    profile.add_seq(aln_b,sol_matrix[sizeA, sizeB])
    
    return profile


class Profile:
    def __init__(self):
        self.alignment_dashboard = {"A":{},"T":{},"G":{},"C":{},"N":{},"-":{}}
        self.score=0
        self.elems_per_col=0
        self.size = 0
        self.seqs = []

    def calculate_score_at(self,index,score_matrix,gap_score):
        nucleotides = []
        score = 0
        for seq in self.seqs:
            try:
                nucleotides.append(seq[index])
            except IndexError:
                print("seq", len(seq))
                print("index", index)
        for x in range(len(nucleotides)):
            for y in range(x+1,len(nucleotides)):
                first= nucleotides[x]
                second= nucleotides[y]
                if not (first == "-" and second == "-"):
                    if first == "-" or second == "-":
                        score+=gap_score
                    else:
                        score+= get_score(score_matrix,first,second)
        return score
                        
                
    def columnAt(self,index):
        a = self.alignment_at("A",index)
        t = self.alignment_at("T",index)
        g = self.alignment_at("G",index)
        c = self.alignment_at("C",index)
        n = self.alignment_at("N",index)
        gap = self.alignment_at("-",index)
        return {"A":a,"T":t,"G":g,"C":c,"N":n,"-":gap}
    def alignment_at(self,nucleotide,index):
        res = 0
        try:
            res =self.alignment_dashboard[nucleotide][index]
        except KeyError:
            return res
        return res
    def traceback(self,ops_matrix,seq):
        gaps_profile=[]
        seq_alineada = []
        i = len(seq)
        j = self.size
        while i > 0 or j > 0:
            if ops_matrix[i, j] == op.GAP_A.value:
                seq_alineada.append("-")
                j -= 1
            elif ops_matrix[i, j] == op.GAP_B.value:
                gaps_profile.append(i -1)
                seq_alineada.append(seq[i -1])
                i -= 1
            else: # es decir, match / mismatch
                i -= 1
                j -= 1
                seq_alineada.append(seq[i])
            
        seq_alineada = seq_alineada[::-1]

        return seq_alineada,gaps_profile
        
    def apply_gaps(self,gaps_positions):
        self.alignment_dashboard = {"A":{},"T":{},"G":{},"C":{},"N":{},"-":{}}
        self.elems_per_col=0
        new_aligned_seq=[]
        for seq in self.seqs:
            aligned_seq=[]
            for index in range(len(seq)):
                if index in gaps_positions:
                    aligned_seq.append("-")
                aligned_seq.append(seq[index])
            
            new_aligned_seq.append(aligned_seq) 
        self.seqs=[]
        for seq in new_aligned_seq:
            self.add_seq(seq,self.score)
    def add_seq(self,aligned_sequense,score):
        for index in range(len(aligned_sequense)):
            actual_character = aligned_sequense[index]
            
            try:
                actual_charater_count = self.alignment_dashboard[actual_character][index]
            except KeyError:
                self.alignment_dashboard[actual_character][index] = 0
                actual_charater_count = self.alignment_dashboard[actual_character][index] 
            self.alignment_dashboard[actual_character][index] = actual_charater_count + 1
            self.size = max(self.size, len(aligned_sequense))
        self.seqs.append( aligned_sequense)   
        self.elems_per_col+=1
        self.score = score
    