# Michał Borowski 406105
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbitblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO, Seq, SeqRecord
from Bio.SearchIO import HmmerIO
import networkx
import obonet
import numpy as np
import os
from goatools.anno.gaf_reader import GafReader
from scipy.stats import fisher_exact
from statsmodels.sandbox.stats.multicomp import multipletests

def create_db():
    makeblastdb_cline = NcbimakeblastdbCommandline(input_file="genes_e_coli_new.fa_.txt", dbtype="nucl", out="DB")
    makeblastdb_cline()

def search_in_db():
    tblastn_cline = NcbitblastnCommandline(db='DB', query="protein_fragments.fa_.txt", out="out.xml", outfmt=5)
    tblastn_cline()

def blast_work():    
    create_db()
    search_in_db()

    handle = open("out.xml")
    blast_records = NCBIXML.parse(handle)
    blast_records = list(blast_records)

    all_records = []
    with open("genes_e_coli_new.fa_.txt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            all_records.append(record)

    result_file = open("answer1.txt", "w")

    best_records = []
    """m = 1e-200"""
    for record in blast_records:
        best = record.alignments[0]
        gene_id = best.title.split(' ')[1]
        protein_id = record.query.split(' ')[0]
        result_file.write(protein_id + " " + gene_id + "\n")

        for r in all_records: #find full record
            if r.id == gene_id:
                r.seq = r.seq.translate()
                best_records.append(r)

        """m = max(m, record.alignments[0].hsps[0].expect)
        print("===========================")
        print(protein_id)
        print("===========================")
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                print(alignment.title.split(' ')[1], hsp.expect)"""
    """print(m)"""

    result_file.close()

    with open("bests.fasta", "w") as output_handle:
        SeqIO.write(best_records, output_handle, "fasta")

    return best_records

def hmm_work():
    os.system('hmmscan Pfam-A.hmm bests.fasta > hmm_output.txt')

    result_file = open("answer2.txt", "w")

    with open('hmm_output.txt') as handle:
        for qresult in HmmerIO.hmmer3_text.Hmmer3TextParser(handle):
            query_id = qresult.id
            best_hit = qresult[0]
            result_file.write(query_id + " " + best_hit.id + "\n")

            """print("===========================")
            print(query_id)
            print("===========================")
            for hit in qresult:
                print(hit.evalue, hit.id)"""

    result_file.close()

def gene_ontology_work(best_records):
    gaf = GafReader("ecocyc.gaf")
    nts = sorted(gaf.associations, key=lambda nt:nt.DB_ID)
    gafs = dict()
    for record in best_records:
        gafs[record.id] = []

    for nt in nts:
        for record in best_records:
            if nt[2] == record.id:
                gafs[record.id].append(nt[4])

    graph = obonet.read_obo("go.obo")

    # Mapping GO to term.
    id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}

    result_file = open("answer3.txt", "w")

    """counts = dict()"""
    for record in best_records:
        terms = set()
        for g in gafs[record.id]:
            terms = terms | networkx.descendants(graph, g) # descendants from current GO

        result_file.write(record.id + "\n")
        terms = sorted(list(terms))
        l = len(terms) - 1
        it = 0
        for term in terms:
            if it != l: 
                result_file.write(id_to_name[term] + ", ")
            else: 
                result_file.write(id_to_name[term])
            it += 1

            """if term in counts:
                counts[term] += 1
            else:
                counts[term] = 1"""
        
        result_file.write("\n")

    """for c in counts:
        print(c, counts[c])"""

    result_file.close()

def pvalue(dictA, dictB, feature):
    Af = 0
    An = 0
    Bf = 0
    Bn = 0

    for a in dictA:
        if feature in dictA[a]:
            Af += 1
        else:
            An += 1

    for b in dictB:
        if feature in dictB[b]:
            Bf += 1
        else:
            Bn += 1

    _, pvalue = fisher_exact([[Af, An], [Bf, Bn]], alternative='greater')
    return pvalue

def fisher_bonferroni(dictA, dictB, features):
    pvalues = np.array([pvalue(dictA, dictB, feature) for feature in features])
    return multipletests(pvalues, method='bonferroni', alpha=0.05) 

def fisher_bonferroni_work():
    handle1 = open("answer1.txt", "r")
    handle2 = open("answer2.txt", "r")

    features = set()
    lines1 = handle1.readlines()
    lines2 = handle2.readlines()
    dictA = dict()
    dictB = dict()
    for i in range(len(lines1)):
        l1 = lines1[i].split(' ')
        l2 = lines2[i].split(' ')
        l2[1] = l2[1][:-1]
        if l1[0][5] == 'A':
            dictA[l1[0]] = {l2[1]}
        else:
            dictB[l1[0]] = {l2[1]}
        features.add(l2[1])

    handle1.close()
    handle2.close()

    pfam_testA = fisher_bonferroni(dictA, dictB, features)[0]
    pfam_testB = fisher_bonferroni(dictB, dictA, features)[0]

    print("Pfam test A:")
    it = 0
    for f in features:
        if pfam_testA[it] == True:
            print(f)
        it += 1

    print("Pfam test B:")
    it = 0
    for f in features:
        if pfam_testB[it] == True:
            print(f)
        it += 1

    handle3 = open("answer3.txt", "r")
    features = set()
    dictA = dict()
    dictB = dict()
    for i in range(len(lines1)):
        l1 = lines1[i].split(' ')
        handle3.readline()
        l3 = handle3.readline().split(',')
        if l1[0][5] == 'A':
            dictA[l1[0]] = set(l3)
        else:
            dictB[l1[0]] = set(l3)
        features = features | set(l3)

    handle3.close()

    go_testA = fisher_bonferroni(dictA, dictB, features)[0]
    go_testB = fisher_bonferroni(dictB, dictA, features)[0]

    print("GO test A:")
    it = 0
    for f in features:
        if go_testA[it] == True:
            print(f)
        it += 1

    print("GO test B:")
    it = 0
    for f in features:
        if go_testB[it] == True:
            print(f)
        it += 1

if __name__ == '__main__':
    best_records = blast_work()
    hmm_work()
    gene_ontology_work(best_records)
    fisher_bonferroni_work()
