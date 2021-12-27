from nw_algorithm_basics import *

def get_score_prof(score_mtx,gap_score,char,column,total_per_column):
    score =0
    for k,v in column.items():
        if k!= "-":
            value= (v/total_per_column)*get_score(score_mtx,char,k)
        else:
            value=(v/total_per_column)*gap_score
        score += value
    return score

def get_score_prof_gap(gap_score,char,column,total_per_column):
    score = 0
    if column == "-":
        return total_per_column * gap_score
    for k,v in column.items():
        if k != "-" and char == "-":
            score += (v/total_per_column)*gap_score
    return score

def nw_msa_util(seq, profile,score_mtx,gap_score):
    sizeA = len(seq)
    sizeB = profile.size
    sol_matrix = np.repeat(0.0, (sizeA+1) * (sizeB+1)).reshape(sizeA+1, sizeB+1)
    ops_matrix = np.repeat(0, (sizeA+1) * (sizeB+1)).reshape(sizeA+1, sizeB+1)
    
    for i in range(1,sizeB+1):
        ops_matrix[0][i] = op.GAP_A.value
    for j in range(1,sizeA+1):
        ops_matrix[j][0] = op.GAP_B.value
        
    for i in range(1,sizeA+1):
        for j in range(1,sizeB+1):
            col=profile.columnAt(j-1)
            nucleotide = seq[i-1]
            
            diagonal_score = get_score_prof(score_mtx,gap_score,nucleotide,col,profile.elems_per_col)
            diagonal_result = sol_matrix[i-1][j-1] + diagonal_score
            
            left_result = sol_matrix[i][j-1] + get_score_prof_gap(gap_score,"-",col,profile.elems_per_col)
            
            above_result = sol_matrix[i-1][j] + get_score_prof_gap(gap_score,nucleotide,"-",profile.elems_per_col)
            
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
 
    aligned_seq,gaps_profile = profile.traceback(ops_matrix,seq)

    return profile, aligned_seq, sol_matrix[sizeA, sizeB],gaps_profile

def nw_msa(seqs,score_mtx, gap_score = 0):
    if len(seqs) == 1 or len(seqs) == 0:
        return
    seq_a = seqs.pop(0)
    seq_b = seqs.pop(0)
    profile = nw_basic(seq_a,seq_b,score_mtx, gap_score)
    for seq in seqs:
        profile,aligned_seq, score,gaps_profile =nw_msa_util(seq,profile,score_mtx,gap_score)
        profile.apply_gaps(gaps_profile)
        profile.add_seq(aligned_seq,score)
    return profile