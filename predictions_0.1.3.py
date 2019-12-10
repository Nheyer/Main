def dict_maker(path, Depth_format_flag, reference_format_flag, alternitive_format_flag, contig):
    import pysam
    import numpy as np
    import sys


    reference_dictionary = {}
    variant_dictionary = {}
    # inport the vcf
    vcf = pysam.VariantFile(path)


    # get the samples from the header
    sample = (list(vcf.header.samples))
    # convert sample names to genotype
    genotype = []
    mult = []
    for name in sample:
        try:
            genotype.append(name.split('_')[1][0])
            if len(name.split('_')) >= 4:
                mult.append(int(name.split('_')[3][0]))
            else:
                mult.append(1)
        except:
            mult.append(0)
            genotype.append('')
    # calculate counts of the given paratypes
    a_type = 0
    b_type = 0
    c_type = 0
    d_type = 0
    n_type = 0
    gfloat = 0.0
    for paratype in range(len(genotype)):
        gfloat += mult[paratype]
        if genotype[paratype] == 'N':
            n_type = n_type + mult[paratype]
        elif genotype[paratype] == 'B':
            b_type = b_type + mult[paratype]
        elif genotype[paratype] == 'A':
            a_type = a_type + mult[paratype]
        elif genotype[paratype] == 'C':
            c_type = c_type + mult[paratype]
        elif genotype[paratype] in ['D', 'R']:
            d_type = d_type + mult[paratype]
    # calc a normilization value from the types
    n_norm = gfloat / (5 * n_type)
    d_norm = gfloat / (5 * d_type)
    a_norm = gfloat / (5 * a_type)
    b_norm = gfloat / (5 * b_type)
    c_norm = gfloat / (5 * c_type)
    print(n_norm, d_norm, a_norm, b_norm, c_norm, gfloat)
    print(n_norm*n_type, d_norm*d_type, a_norm*a_type, b_norm*b_type, c_norm*c_type)
    # loop through the vcf
    for var in vcf.fetch(contig=contig):
        if var.alts is not None:
            pos_key = str(var.pos) + "-".join(var.alts)
        elif var.alts is None and var.ref is not None:
            pos_key = str(var.pos) + "-".join(var.ref)
        else:
            sys.stderr.write("NO REF OR ALT ALLILE found at pos:\t%s\n exiting!!\n" % var.pos)
            sys.exit(-1)

        # initiolize number of reads of each genotypes seen across all samples
        ref_N = 0
        ref_A = 0
        ref_B = 0
        ref_C = 0
        ref_D = 0
        var_N = 0
        var_A = 0
        var_B = 0
        var_C = 0
        var_D = 0

        # for all samples in the variant line
        for i in range(len(var.samples)):
            if Depth_format_flag and reference_format_flag in var.samples[i]:  # Potential to blacklist things here
                # initiolize a varyable for the number of reads
                rr = 0.0
                vr = 0.0
                # for help with debuging try to do the next part if you can't print what you saw
                    # if this pos did not have a variant in the vcf assume it was reference here, if not use the full depth and number of ref observed to get proportions.
                if str(var.samples[i][reference_format_flag]) in (None , 'None'):
                    continue
                else:
                    rr = float(var.samples[i][reference_format_flag])/float(var.samples[i][Depth_format_flag])
                if alternitive_format_flag not in var.samples[i] or str(var.samples[i][alternitive_format_flag][0]) == 'None':
                    vr = 0.0
                else:
                    vr = float(var.samples[i][alternitive_format_flag][0])/float(var.samples[i][Depth_format_flag])
                # now for the current sample add the number of reads to the prober bin based on assigned genotype
                # and normolize to make each paratype 1/5th of dictionary
                if genotype[i] == 'N':
                    ref_N = ref_N + n_norm * rr * mult[i]
                    var_N = var_N + n_norm * vr * mult[i]
                elif genotype[i] == 'A':
                    ref_A = ref_A + a_norm * rr * mult[i]
                    var_A = var_A + a_norm * vr * mult[i]
                elif genotype[i] == 'B':
                    ref_B = ref_B + b_norm * rr * mult[i]
                    var_B = var_B + b_norm * vr * mult[i]
                elif genotype[i] == 'C':
                    ref_C = ref_C + c_norm * rr * mult[i]
                    var_C = var_C + c_norm * vr * mult[i]
                elif genotype[i] == 'D':
                    ref_D = ref_D + d_norm * rr * mult[i]
                    var_D = var_D + d_norm * vr * mult[i]
            # get a total sum for the reference totalthen it it is not all zero calculate proportions and put it in a dictionary
            ref_sum = ref_N + ref_A + ref_B + ref_C + ref_D
            if ref_sum > 0.0:
                reference_dictionary[pos_key] = np.array([ref_N, ref_A, ref_B, ref_C, ref_D]) / ref_sum
            # do the same for variant reads.
            var_sum = var_N + var_A + var_B + var_C + var_D
            if var_sum > 0.0:
                variant_dictionary[pos_key] = np.array([var_N, var_A, var_B, var_C, var_D]) / var_sum
    # close the vcf
    vcf.close()
    # return a list of two dictionarys that contain the proportion of genotypesin a position that is reference or variant
    return [variant_dictionary,reference_dictionary]


def find_total_genotype(path, dictionarys, dp, ref, alt, contig):
    import numpy as np
    import pysam

    vars = dictionarys[0]
    refs = dictionarys[1]
    # create running probubility score for the sample
    run_sample_prop_score = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
    vcf = pysam.VariantFile(path)

    # look at all the variants in the file on contig supplied
    for variant in vcf.fetch(contig):
            # get out current position
            if variant.alts is not None:
                pos_key = str(variant.pos) + "-".join(variant.alts)
            elif variant.alts is None and variant.ref is not None:
                pos_key = str(variant.pos) + "-".join(variant.ref)
            else:
                sys.stderr.write("NO REF OR ALT ALLILE found at pos:\t%s\n exiting!!\n" % variant.pos)
                sys.exit(-1)
            # find proportion of reads that were reference and proportion that were variant
            if variant.samples[0][dp] is not None and variant.samples[0][dp] > 0:
                if alt in variant.samples[0] and variant.samples[0][alt][0] is not None:
                    var_scaling_factor = float(variant.samples[0][alt][0])/float(variant.samples[0][dp])
                else:
                    var_scaling_factor = 0
                if variant.samples[0][ref] is not None:
                    ref_scaling_factor = float(variant.samples[0][ref])/float(variant.samples[0][dp])
                else:
                    ref_scaling_factor = 0
            else:
                continue
            # if we have this position in out varinat dictionary... add a scaled value of it to our running score
            if pos_key in vars:
                run_sample_prop_score = run_sample_prop_score + var_scaling_factor * vars[pos_key]
            if pos_key in refs:
                run_sample_prop_score = run_sample_prop_score + ref_scaling_factor * refs[pos_key]
    # return a vector of cumulative proportions based on all positions
    return run_sample_prop_score/sum(run_sample_prop_score)


def props_to_guess(props):
    # import packages
    import numpy as np
    # create a list of the true possabilities for our data
    allowable_NOTCH = [
        ["NNAAAA", np.array([0.333333333333, 0.666666666667, 0.0, 0.0, 0.0])],
        ["NNAAAAD", np.array([0.285714285714, 0.571428571429, 0.0, 0.0, 0.142857142857])],
        ["NNAAAADD", np.array([0.25, 0.5, 0.0, 0.0, 0.25])],
        ["NNAAAAC", np.array([0.285714285714, 0.571428571429, 0.0, 0.142857142857, 0.0])],
        ["NNAAAACD", np.array([0.25, 0.5, 0.0, 0.125, 0.125])],
        ["NNAAAACDD", np.array([0.222222222222, 0.444444444444, 0.0, 0.111111111111, 0.222222222222])],
        ["NNAAAACC", np.array([0.25, 0.5, 0.0, 0.25, 0.0])],
        ["NNAAAACCD", np.array([0.222222222222, 0.444444444444, 0.0, 0.222222222222, 0.111111111111])],
        ["NNAAAACCDD", np.array([0.2, 0.4, 0.0, 0.2, 0.2])],
        ["NNAAAB", np.array([0.333333333333, 0.5, 0.166666666667, 0.0, 0.0])],
        ["NNAAABD", np.array([0.285714285714, 0.428571428571, 0.142857142857, 0.0, 0.142857142857])],
        ["NNAAABDD", np.array([0.25, 0.375, 0.125, 0.0, 0.25])],
        ["NNAAABC", np.array([0.285714285714, 0.428571428571, 0.142857142857, 0.142857142857, 0.0])],
        ["NNAAABCD", np.array([0.25, 0.375, 0.125, 0.125, 0.125])],
        ["NNAAABCDD", np.array([0.222222222222, 0.333333333333, 0.111111111111, 0.111111111111, 0.222222222222])],
        ["NNAAABCC", np.array([0.25, 0.375, 0.125, 0.25, 0.0])],
        ["NNAAABCCD", np.array([0.222222222222, 0.333333333333, 0.111111111111, 0.222222222222, 0.111111111111])],
        ["NNAAABCCDD", np.array([0.2, 0.3, 0.1, 0.2, 0.2])],
        ["NNAABB", np.array([0.333333333333, 0.333333333333, 0.333333333333, 0.0, 0.0])],
        ["NNAABBD", np.array([0.285714285714, 0.285714285714, 0.285714285714, 0.0, 0.142857142857])],
        ["NNAABBDD", np.array([0.25, 0.25, 0.25, 0.0, 0.25])],
        ["NNAABBC", np.array([0.285714285714, 0.285714285714, 0.285714285714, 0.142857142857, 0.0])],
        ["NNAABBCD", np.array([0.25, 0.25, 0.25, 0.125, 0.125])],
        ["NNAABBCDD", np.array([0.222222222222, 0.222222222222, 0.222222222222, 0.111111111111, 0.222222222222])],
        ["NNAABBCC", np.array([0.25, 0.25, 0.25, 0.25, 0.0])],
        ["NNAABBCCD", np.array([0.222222222222, 0.222222222222, 0.222222222222, 0.222222222222, 0.111111111111])],
        ["NNAABBCCDD", np.array([0.2, 0.2, 0.2, 0.2, 0.2])],
        ["NNABBB", np.array([0.333333333333, 0.166666666667, 0.5, 0.0, 0.0])],
        ["NNABBBD", np.array([0.285714285714, 0.142857142857, 0.428571428571, 0.0, 0.142857142857])],
        ["NNABBBDD", np.array([0.25, 0.125, 0.375, 0.0, 0.25])],
        ["NNABBBC", np.array([0.285714285714, 0.142857142857, 0.428571428571, 0.142857142857, 0.0])],
        ["NNABBBCD", np.array([0.25, 0.125, 0.375, 0.125, 0.125])],
        ["NNABBBCDD", np.array([0.222222222222, 0.111111111111, 0.333333333333, 0.111111111111, 0.222222222222])],
        ["NNABBBCC", np.array([0.25, 0.125, 0.375, 0.25, 0.0])],
        ["NNABBBCCD", np.array([0.222222222222, 0.111111111111, 0.333333333333, 0.222222222222, 0.111111111111])],
        ["NNABBBCCDD", np.array([0.2, 0.1, 0.3, 0.2, 0.2])],
        ["NNBBBB", np.array([0.333333333333, 0.0, 0.666666666667, 0.0, 0.0])],
        ["NNBBBBD", np.array([0.285714285714, 0.0, 0.571428571429, 0.0, 0.142857142857])],
        ["NNBBBBDD", np.array([0.25, 0.0, 0.5, 0.0, 0.25])],
        ["NNBBBBC", np.array([0.285714285714, 0.0, 0.571428571429, 0.142857142857, 0.0])],
        ["NNBBBBCD", np.array([0.25, 0.0, 0.5, 0.125, 0.125])],
        ["NNBBBBCDD", np.array([0.222222222222, 0.0, 0.444444444444, 0.111111111111, 0.222222222222])],
        ["NNBBBBCC", np.array([0.25, 0.0, 0.5, 0.25, 0.0])],
        ["NNBBBBCCD", np.array([0.222222222222, 0.0, 0.444444444444, 0.222222222222, 0.111111111111])],
        ["NNBBBBCCDD", np.array([0.2, 0.0, 0.4, 0.2, 0.2])]
    ]
    # create five gating varyables, these will also act as our a=maximum difference threshold
    gate = [0.5,0.6,0.7,0.8,0.9]
    # make a place to store the  values
    return_string = ['','','','','']
    # loop through the list
    for posability in allowable_NOTCH:
        if np.sum(np.absolute(posability[1] - props)) < gate[0]:
            # this line is a new top value shift things down
            gate[4] = gate[3]
            return_string[4] = return_string[3]
            gate[3] = gate[2]
            return_string[3]=return_string[2]
            gate[2] = gate[1]
            return_string[2] = return_string[1]
            gate[1] = gate[0]
            return_string[1] = return_string[0]
            # now we can set out top values
            gate[0] = np.sum(np.absolute(posability[1] - props))
            return_string[0] = posability[0]
        elif np.sum(np.absolute(posability[1] - props)) < gate[1]:
            # we have a new value for second shift things down
            gate[4] = gate[3]
            return_string[4] = return_string[3]
            gate[3] = gate[2]
            return_string[3] = return_string[2]
            gate[2] = gate[1]
            return_string[2] = return_string[1]
            # now we have space we can set the value
            gate[1] = np.sum(np.absolute(posability[1] - props))
            return_string[1] = posability[0]
        elif np.sum(np.absolute(posability[1] - props)) < gate[2]:
            # we have a new third value shif things down
            gate[4] = gate[3]
            return_string[4] = return_string[3]
            gate[3] = gate[2]
            return_string[3] = return_string[2]
            # now we have space we can set the value
            gate[2] = np.sum(np.absolute(posability[1] - props))
            return_string[2] = posability[0]
        elif np.sum(np.absolute(posability[1] - props)) < gate[3]:
            # we hav a new fourth value shift things down
            gate[4] = gate[3]
            return_string[4] = return_string[3]
            # now we have space we can set the value
            gate[3] = np.sum(np.absolute(posability[1] - props))
            return_string[3] = posability[0]
        elif np.sum(np.absolute(posability[1] - props)) < gate[4]:
            gate[4] = np.sum(np.absolute(posability[1] - props))
            return_string[4] = posability[0]
    # return a list of lists with the decided values and the difference scores
    return [return_string,gate]


def import_info():
    import sys
    # loop through all the inputs
    unk_input_path = ''
    ref_input_path = ''
    ful_flag = 'DP'
    ref_flag = 'RO'
    ctig = "NOTCH2NL-consensus"
    alt_flag = "AO"
    # loop through the arguments till we hind a flag then assign the next term to the appropriate varyable.
    for arg_index in range(len(sys.argv)):
        if sys.argv[arg_index] in ['-h','--help']:
            sys.stderr.write('''
            Usage python predictions.py [options] -i/--input -r/--reference \n
            ________________________________________________________________________
            Options:\n
            -i\t--input     \t defines the single sample vcf that the program will make inferences about 
            -r\t--reference \t defines the multisample vcf used to make known data, the sample names should have the pariotype leter after the first '_'
            -A\t--alt-depth \t defines the format flag for the depth of reads supporting the alternitive in the dictionary  ; default == 'AO'
            -R\t--ref-depth \t defines the format flag for the depth of reads supporting the reference in the dictionary ; default == 'RO'
            -D\t--full-depth\t defines the format flag for the depth of reads in the dictionary  ; default == 'DP'
            -c\t--contig    \t define contig to grab vcf records from default == NOTCH2NL-consensus
            -h\t--help      \t print this information
            ''')
            sys.exit()
        elif sys.argv[arg_index] in ['-i','--input']:
            unk_input_path = sys.argv[arg_index+1]
        elif sys.argv[arg_index] in ['-r','--reference']:
            ref_input_path = sys.argv[arg_index+1]
        elif sys.argv[arg_index] in ['-A','--alt-depth']:
            alt_flag = sys.argv[arg_index+1]
            sys.stderr.write('setting alternative depth format flag to %s\n' % (sys.argv[arg_index+1]))
        elif sys.argv[arg_index] in ['-c','--contig']:
            ctig = sys.argv[arg_index + 1]
        elif sys.argv[arg_index] in ['-R','--ref-depth']:
            ref_flag = sys.argv[arg_index+1]
            sys.stderr.write('setting referenve depth format flag to %s\n' % (sys.argv[arg_index+1]))
        elif sys.argv[arg_index] in ['-D', '--full-depth']:
            ful_flag = sys.argv[arg_index + 1]
            sys.stderr.write('setting referenve depth format flag to %s\n' % (sys.argv[arg_index + 1]))


    return ref_input_path, ful_flag, ref_flag, alt_flag, unk_input_path, ctig


def main():
    # allow useer to pass in needed info
    ref_vcf, dp_flag, ref_flag, alt_flag, unknown, contig = import_info()
    print("ref vcf:\t%s\ndepth-ref-alt:\t%s-%s-%s\nUnknown vcf:\t%s\nContig to look at:\t%s" %
          (ref_vcf, dp_flag, ref_flag, alt_flag, unknown, contig)
          )
    # calculate the reference dictionaries
    dict_list = dict_maker(ref_vcf, dp_flag, ref_flag, alt_flag, contig)
    sample_prop = find_total_genotype(unknown, dict_list, dp_flag, ref_flag, alt_flag, contig)
    print(sample_prop)
    datamtx = props_to_guess(sample_prop)
    # generate our own man
    decode = ["N","A","B","C","D"]
    guess = ""
    for i in range(len(sample_prop)):
        while True:
            if sample_prop[i] > 0.1:
                sample_prop[i] -= 0.1
                guess += decode[i]
            else:
                break
        print("Current guess:\t", guess)
        print(sample_prop)
    print("I think it is:\t ", guess, "\n but off by:\t", sum(sample_prop))

    outstring = '''
    %s \t with a score of : %.3f\n
    %s \t with a score of : %.3f\n
    %s \t with a score of : %.3f\n
    %s \t with a score of : %.3f\n
    %s \t with a score of : %.3f\n
    ''' % (datamtx[0][0],datamtx[1][0],
           datamtx[0][1],datamtx[1][1],
           datamtx[0][2],datamtx[1][2],
           datamtx[0][3],datamtx[1][3],
           datamtx[0][4],datamtx[1][4]
           )
    print(outstring)


if __name__ == '__main__':
    main()