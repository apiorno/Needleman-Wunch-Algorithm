import numpy as np
from Bio import SeqIO
from grasp_nw import *
from nw_pdf_creator import PDF

def score_matrix():
    score_matrix = np.repeat(-1, (5) * (5)).reshape(5, 5)
    for i in range(0,5):
        score_matrix[i][i]=1
    return score_matrix

def get_seqs(file_name,file_type):
    seqs = []
    file_align = SeqIO.parse(file_name, file_type)

    while(True):
        try:
            elem = next(file_align)
            seq = str(elem.seq)
            seqs.append(seq)
        except StopIteration:
            break
def main():
    score_mtx = score_matrix()
    gap_score = 0
    seqs = get_seqs("x.fasta",'fasta')
    itarations = 20
    result=grasp_nw(seqs,score_mtx,gap_score,itarations)

    iterations = []
    incs = []
    scores = []
    for k,t in result.items():
        (inc,score,prof) = t
        iterations.append(k)
        scores.append(score)
        incs.append(inc)
        profile = prof

    with open('seqs.txt', 'w') as f:
        index = 1
        for line in profile.seqs:
            str_line= ''.join(line)
            f.write(f'Aligned sequense {index}:')
            f.write('\n')
            f.write(str_line)
            f.write('\n')
            index += 1
                
    plt.plot(iterations, incs, label = 'Increments')
    plt.xlabel('iteration')
    plt.ylabel('score increment')
    plt.legend()
    plt.savefig('chart1.png')

    plt.plot(iterations, scores, label = 'Scores')
    plt.xlabel('iteration')
    plt.ylabel('score')
    plt.legend()
    plt.savefig('chart2.png')

    pdf = PDF(orientation='P', unit='mm', format='A4')

    pdf.add_page()
    pdf.lines()
    pdf.titles()
    pdf.charts('chart1.png',40.0,25.0)
    pdf.charts('chart2.png',40.0,100.0)
    pdf.add_page()
    pdf.texts('seqs.txt')
    pdf.output('out.pdf','F')



if __name__ == "__main__":
    main()
