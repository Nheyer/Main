import sys
import os
debug = False


# translates starting at first M, assuming the 3n+0 ORF, also some matching logic
def translate(dna_seq, name_base, make_dna, pos, min_length, match_seq):
    poly_peptide_chain = ""
    edit_dna_seq = ""
    num_start = 0
    in_mrna = False
    transcription_initiator = "ATG"
    name_chain_llist = []
    seqs_match = match_seq.split(",")

    trna = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': 'Ochre', 'TAG': 'Amber',
        'TGC': 'C', 'TGT': 'C', 'TGA': 'Opal', 'TGG': 'W',
    }
    # start translation
    for codon_start in range(0, len(dna_seq)-3, 3):
        codon = dna_seq[codon_start] + dna_seq[codon_start + 1] + dna_seq[codon_start + 2]
        if codon == transcription_initiator and not in_mrna:
            in_mrna = True
            num_start = codon_start
        if not in_mrna:
            pass
        else:
            amino_acid = trna[codon]
            # if we don't have a stop codeon add on to the chain
            if amino_acid not in ['Ochre', 'Opal', 'Amber']:
                poly_peptide_chain += amino_acid
                if make_dna:
                    edit_dna_seq += codon
            # if it is a stop add the polypeptide to the chane and print out name found
            else:
                chain_name = name_base + "_start_" + str(num_start) + "_endpos_" + \
                             str(codon_start) + "_stop-type_" + amino_acid
                ###print( "Found " + chain_name +"\n")
                if len(poly_peptide_chain) > min_length:
                    if make_dna:
                        poly_peptide_chain = edit_dna_seq
                    name_chain_llist.append([chain_name, poly_peptide_chain, False])
                #else:
                #    sys.stderr.write("Found a chain of " + str(len(poly_peptide_chain))
                #                     + " aminoacids and discarding it \n")
                # reset things
                poly_peptide_chain = ""
                edit_dna_seq = ""
                in_mrna = False
                if not pos:
                    num_start += 1
    # add sequence match bool to look at later
    for chain in name_chain_llist:
        for test_string in seqs_match:
            if test_string in chain[1]:
                chain[2] = True
                break
    return name_chain_llist


# dose math to calculate  desired slice inputs for sequences
def calc(x, alpha):
    return int(((alpha/2)*(x*x))+(((3*alpha + 2)/2)*x)+alpha)


# dose what is says on the tin...
def make_frame_edited_sequences(parent_seq,par_name):
    DNA_pol_III = {'A': 'T', 'T': 'A', 'U': 'A', 'C': 'G', 'G': 'C'}
    # find the num over a codon we are, save it as over_baces
    len_pseq = len(parent_seq)
    over_bases = len(parent_seq) % 3

    # get backwards complement
    backwards_seq = ''
    for bp in parent_seq[::-1]:
        backwards_seq += DNA_pol_III[bp]

    # get splice formatted end of seq
    s_key_0 = calc(-over_bases,len_pseq)
    s_key_1 = calc(-abs(1-over_bases), len_pseq)
    s_key_2 = calc(-abs(2-over_bases), len_pseq)

    # return sequences in the right ORFs
    return [[par_name + '_0+3n', parent_seq[0:s_key_0:]],
            [par_name + '_1+3n', parent_seq[1:s_key_1:]],
            [par_name + '_2+3n', parent_seq[2:s_key_2:]],
            [par_name + '_0-3n', backwards_seq[0:s_key_0:]],
            [par_name + '_1-3n', backwards_seq[1:s_key_1:]],
            [par_name + '_2-3n', backwards_seq[2:s_key_2:]]]


# as inputs are .fa parse them to usable formats.
def fasta_parse(path):
    fata_file = open(path)
    startup = True
    fasta_llist = []
    t_list = []

    for line in fata_file:
        print(line)
        if line[0] == '>' and startup:
            startup = False
            t_list = [line.rstrip("\n"),'']
        elif line[0] == ">":
            fasta_llist.append(t_list)
            t_list = [line.rstrip("\n"), '']
        elif not startup:
            t_list[1] += line.rstrip("\n")
    fasta_llist.append(t_list)
    fata_file.close()
    if debug:
        print(fasta_llist)
    return fasta_llist


# parse user arguments
def io():
    help_stm = '''
                Usage: Python3 exons-to-polypep.py -i <path to exons> -o <path to output fasta>
                ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n
                Optional flags:
                -h  \t --help       \t Print this usage statement
                -p  \t --position   \t add position information to print out
                -D  \t --DNA        \t output DNA instead of a polypeptide chain
                -c  \t --min-count  \t an int giving the minimum amino acids needed to consider the chain 
                -o2 \t --output-two \t put sequences that do not contain the matching seqs here instead of 
                    \t              \t making own name default == <path to output fasta>_not_matching.fa
                -m  \t --matching   \t a comma separated list of sequence that all true sequences should contain
                    \t              \t this will cause a second output tyo be created use the -o2 flag to name it
                -e  \t --ending     \t amino acid sequence to check for at the end to add in 
                -E  \t --p-log      \t only use '-e' on files named with a comma separated list 
                -L  \t --log        \t output log info to file instead of stderr
                '''
    input_f = ""
    output_f = ""
    output_2_f = ""
    DNA_f = False
    position_f = False
    min_length_f = 64
    output_2_f = ""
    reg_ex_f = "GTGYCKCPEGFLGEYCQHRD,PCEKNRCQNGGTCVAQAMLG,KATCRCASGFTGEDCQYST," + \
               "PCLNGGTCHMLSRDTYECTC,QVGFTGKECQWTDACLSH," + \
               "BGSTCTTVANQFSCKCLTGF,TGFTGQKCETDVNECDIPGHCQHGG"
    ends = "KRTR,KRTR,KEHDEN,KEHDEN,KEHDEN"
    ptyps = "NOTCH2,NOTCH2NL-D,NOTCH2NL-A,NOTCH2NL-B,NOTCH2NL-C"
    endings_f = []
    log_f = ""
    # parse the input
    for n in range(len(sys.argv)):
        if sys.argv[n] in ['-i', '--input']:
            input_f = sys.argv[n + 1]
        elif sys.argv[n] in ['-o', '--output']:
            output_f = sys.argv[n + 1]
        elif sys.argv[n] in ['-D', '--DNA']:
            DNA_f = True
        elif sys.argv[n] in ['-p', '--position']:
            position_f = True
        elif sys.argv[n] in ['-o2', '--output-two']:
            output_2_f = sys.argv[n + 1]
        elif sys.argv[n] in ['-c', '--min-count']:
            min_length_f = int(sys.argv[n + 1])
        elif sys.argv[n] in ['-m', '--matching']:
            if sys.argv[n + 1][0] == '-':
                reg_ex_f = ""
            else:
                reg_ex_f = sys.argv[n + 1]
        elif sys.argv[n] in ['-e', '--ending']:
            ends = sys.argv[n + 1]
        elif sys.argv[n] in ['-E', '--p-log']:
            ptyps = sys.argv[n + 1]
        elif sys.argv[n] in ['-L', '--log']:
            log_f = sys.argv[n + 1]
        elif sys.argv[n] in ['-h', '--help'] or len(sys.argv) < 3:
            sys.stderr.write(help_stm)
            sys.exit(0)
    if debug:
        sys.stderr.write(input_f)
    if output_2_f == "":
        output_2_f = output_f.rstrip(".fa") + "_not_matching.fa"

    if ends != "":
        ends_l = ends.split(",")
        for strings in ends_l:
            endings_f.append(["", len(strings), strings])
        if ptyps != "":
            ptyps_l = ptyps.split(",")
            if len(ends_l) == len(ptyps_l):
                for i in range(len(ptyps_l)):
                    endings_f[i][0] = ptyps_l[i]
            else:
                sys.stderr.write("Expected same length for '-e' and '-E' flags\n")
                sys.exit(-1)
    return input_f, output_f, output_2_f, DNA_f, position_f, min_length_f, reg_ex_f, endings_f , log_f


# use check the endings to see if we should include this sequence
def check_ends(pt_seq, end_ll):
    # loop through all endings
    for pt_len_ed in end_ll:
        # if there is a ptype matching ( null counts as matching )
        if pt_len_ed[0] in pt_seq[0]:
            sys.stdout.write("\n" + pt_seq[0] + "\t" + pt_seq[1][-pt_len_ed[1]:] + "\t" + pt_len_ed[2]
                             + "\t" + str(pt_seq[2]))
            # see if the last positions match
            if pt_seq[1][-pt_len_ed[1]:] == pt_len_ed[2]:
                sys.stdout.write("\tTrue")
                return True
    return pt_seq[2]


def main():

    input_path, output_path, alt_output, p_dna, pos_dna, min_len, regx, endings, log_path = io()
    out_fasta = open(output_path, 'w')
    # open a second flie if we are matching, so we can see what didn't match
    if regx != "":
        out_fasta_2 = open(alt_output, 'w')
    if log_path:
        log_file = open(log_path,"w")

    # loop through the sequences in the fasta
    for fasta_list in fasta_parse(input_path):
        foud_an_orf = False
        # loop through the 6 frames
        for pair in make_frame_edited_sequences(fasta_list[1], fasta_list[0]):
            # loop through all the polypeptide chains
            for name_seq in translate(pair[1], pair[0], p_dna, pos_dna, min_len, regx):
                # do ending checks if asked
                if len(endings) > 0:
                    name_seq[2] = check_ends(name_seq, endings)
                # write the polypeptide chain out in fasta format
                if name_seq[2]:
                    foud_an_orf = True
                    out_fasta.write(name_seq[0] + "\n")
                else:
                    out_fasta_2.write(name_seq[0] + "\n")
                for i in range(0, len(name_seq[1]), 78):
                    # see if we found one of the strings to match, if so put it in primary
                    # output, then the all else into second.
                    if name_seq[2]:
                        out_fasta.write(name_seq[1][i:i+78] + "\n")
                    else:
                        out_fasta_2.write(name_seq[1][i:i+78] + "\n")
        if not foud_an_orf:
            if log_path:
                log_file.write(fasta_list[0] + " has no polypeptides\n")
            else:
                sys.stderr.write(fasta_list[0] + " has no polypeptides\n")

    # only close the second file if we opened it
    if regx != "":
        out_fasta_2.close()
    if log_path:
        log_file.close()
    out_fasta.close()


if __name__ == '__main__':
    main()
