import random
from nw_algorithm import *

def generate_random_solution_with_profile(seqs,profile,score_matrix,gap_score=0):
    while len(seqs) != 0:
        size = len(seqs)
        solutions=[]
        for seq in seqs:
            prof, aligned_seq, score,gaps_profile=nw_msa_util(seq, profile,score_matrix,gap_score)
            solutions.append((aligned_seq, score,gaps_profile,seq))
        solutions.sort(key=lambda x:x[1], reverse=True)
        first=0
        last = size-1
        mid=(size//2) if size>1 else 0
        aligned_seq, score,gaps_profile,original_seq = solutions[random.choice([first,mid,last])]
        profile.apply_gaps(gaps_profile)
        profile.add_seq(aligned_seq,score)
        seqs.remove(original_seq)
    return profile

def local_max_from_neighbourhood(seq_index,profile,score_matrix,gap_score):
   
    actual_seq = profile.seqs[seq_index]
    changes = []
    for index in range(len(actual_seq)-1):
                
            first = actual_seq[index]
            second = actual_seq[index+1]
            if ((first == "-" or second == "-") and first!=second):
                new_seq = actual_seq[:index]
                new_seq.append(second)
                new_seq.append(first)
                new_seq += actual_seq[index+2:]
                first_score = profile.calculate_score_at(index,score_matrix,gap_score)
                second_score = profile.calculate_score_at(index+1,score_matrix,gap_score)
                    
                profile.seqs[seq_index] = new_seq
                new_first_score = profile.calculate_score_at(index,score_matrix,gap_score)
                new_second_score = profile.calculate_score_at(index+1,score_matrix,gap_score)
                profile.seqs[seq_index] = actual_seq
                diff= new_first_score + new_second_score - first_score -second_score
                if (diff>0):
                    changes.append((index,index+1,new_seq,diff))
                    
                    
    if len(changes)>0:
        change = max(changes, key=lambda p: p[3])
        profile.seqs[seq_index] = change[2]
        profile.score = profile.score +  change[3]
        return (change[3],profile.score,profile)
    else:
        return (0,profile.score,profile)

def generate_random_solution(seqs,score_matrix,gap_score=0):
    size = len(seqs)
    solutions = []
    for row in range(size):
        for col in range(row+1,size):
            seq_a = seqs[row]
            seq_b = seqs[col]
            solutions.append((nw_basic(seq_a,seq_b,score_matrix,0),seq_a,seq_b))
    solutions.sort(key=lambda x:x[0].score, reverse=True)
    first=0
    last = size-1
    mid= (size//2)  if size>1 else 0
    random_index = random.choice([first,mid,last])
    profile,seq_a,seq_b = solutions[random_index]
    seqs.remove(seq_a)
    seqs.remove(seq_b)
    return generate_random_solution_with_profile(seqs,profile,score_matrix,gap_score)    

def look_for_better_solution(profile,score_matrix,gap_score,max_iterations):

    max_results_per_it = {}
    max_results_per_it[0]= (0,profile.score,profile)
    for actual_it in range(1,max_iterations+1):
        local_max_profiles= []
        for index in range(len(profile.seqs)):
            iter_changes = local_max_from_neighbourhood(index,profile,score_matrix,gap_score)
            local_max_profiles.append(iter_changes)
        max_results_per_it[actual_it] = local_max_profiles[len(local_max_profiles)-1]   
    return max_results_per_it


def grasp_nw(seqs,score_matrix,gap_score,max_iterations):
    p=generate_random_solution(seqs,score_matrix,gap_score)
    return look_for_better_solution(p,score_matrix,gap_score,max_iterations)


    