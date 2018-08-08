def dict_maker(path,variant_format_flag,reference_format_flag):
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
    for name in sample:
        try:
            genotype.append(name.split('_')[1][0])
        except:
            genotype.append('')
    # loop through the vcf
    for var in vcf:
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
            # initiolize a varyable for the number of reads
            rr = 0.0
            vr = 0.0
            # for help with debuging try to do the next part if you cant prit what you saw
            try:
                # see if ether the variant reads or the reeference reads are none if so set it to zero instead of casting.
                if str(var.samples[i][reference_format_flag]) == 'None':
                    rr = 0.0
                else:
                    rr = float(var.samples[i][reference_format_flag])
                if str(var.samples[i][variant_format_flag]) == 'None':
                    vr = 0.0
                else:
                    try:
                        vr = float(var.samples[i][variant_format_flag])
                    except:
                        vr = 0.0
                        for val in var.samples[i][variant_format_flag]:
                            if val == None:
                                vr = vr + 0.0
                            else:
                                vr = vr + float(val)

            except:
                print(var.samples[i][reference_format_flag])
                print(var.samples[i][variant_format_flag])
            # now for the current sample add the number of reads to the prober bin based on assigned genotype
            if genotype[i] == 'N':
                ref_N = ref_N + rr
                var_N = var_N + vr
            elif genotype[i] == 'A':
                ref_A = ref_A + rr
                var_A = var_A + vr
            elif genotype[i] == 'B':
                ref_B = ref_B + rr
                var_B = var_B + vr
            elif genotype[i] == 'C':
                ref_C = ref_C + rr
                var_C = var_C + vr
            elif genotype[i] == 'D':
                ref_D = ref_D + rr
                var_D = var_D + vr
        # get a total sum for the reference totalthen it it is not all zero calculate proportions and put it in a dictionary
        ref_sum = ref_N + ref_A + ref_B + ref_C + ref_D
        if ref_sum > 0.0:
            reference_dictionary[var.pos] = np.array([ref_N , ref_A , ref_B , ref_C , ref_D]) / ref_sum
        # do the same for variant reads.
        var_sum = var_N + var_A + var_B + var_C + var_D
        if var_sum > 0.0:
            variant_dictionary[var.pos] = np.array([var_N, var_A, var_B, var_C, var_D]) / var_sum
    # close the vcf
    vcf.close()
    # return a list of two dictionarys that contain the proportion of genotypesin a position that is reference or variant
    return [variant_dictionary,reference_dictionary]


def find_total_genotype(path,dictionarys):
    import numpy as np
    import pysam
    vars = dictionarys[0]
    refs = dictionarys[1]
    # create running probubility score for the sample
    run_sample_prop_score = np.array([0.0,0.0,0.0,0.0,0.0])
    vcf = pysam.VariantFile(path)

    # look at all the variants in the file...
    for variant in vcf:
        # actualy only the ones that are from chromosome 1
        if True:
            # get out current position
            pos = variant.pos
            # find proportion of reads that were reference and proportion that were variant
            var_scaling_factor = sum(list(map(float,variant.info['AO'])))/float(variant.info['DP'])
            ref_scaling_factor = float(variant.info['RO'])/float(variant.info['DP'])
            # if we have this position in out varinat dictionary... add a scaled value of it to our running score
            if vars.has_key(pos):
                run_sample_prop_score = run_sample_prop_score + var_scaling_factor * vars[pos]
            if refs.has_key(pos):
                run_sample_prop_score = run_sample_prop_score + ref_scaling_factor * refs[pos]
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
    gate = [0.5,0.5,0.5,0.5,0.5]
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
    alt_flag = 'AD'
    ref_flag = 'RD'
    # loop through the arguments till we hind a flag then assign the next term to the appropriate varyable.
    for arg_index in range(len(sys.argv)):
        if sys.argv[arg_index] in ['-h','--help']:
            sys.stderr.write('''
            Usage python predictions.py [options] -i/--input -r/--reference \n
            ________________________________________________________________________
            Options:\n
            -i\t--input\t defines the single sample vcf that the program will make inferences about 
            -r\t--reference\t defines the multisample vcf used to make known data, the sample names should have the pariotype leter after the first '_'
            -A\t--alt-depth\t defines the format flag for the depth of reads supporting the alternte in the dictionary (assumes it is a collated depth of all alternatives) ; default == 'AD'
            -R\t--ref-depth\t defines the format flag for the depth of reads supporting the reference in the dictionary ; default == 'RD'
            -h\t--help\tprint this information
            ''')
            sys.exit()
        elif sys.argv[arg_index] in ['-i','--input']:
            unk_input_path = sys.argv[arg_index+1]
        elif sys.argv[arg_index] in ['-r','--reference']:
            ref_input_path = sys.argv[arg_index+1]
        elif sys.argv[arg_index] in ['-A','--alt-depth']:
            alt_flag = sys.argv[arg_index+1]
            sys.stderr.write('setting alternative depth format flag to %s\n' % (sys.argv[arg_index+1]))
        elif sys.argv[arg_index] in ['-R','--ref-depth']:
            ref_flag = sys.argv[arg_index+1]
            sys.stderr.write('setting referenve depth format flag to %s\n' % (sys.argv[arg_index+1]))

    return [ref_input_path,alt_flag,ref_flag,unk_input_path]


def main():
    # allow useer to pass in needed info
    IOinfo = import_info()
    print (IOinfo)
    # calculate the reference dictionaries
    dict_list = dict_maker(IOinfo[0],IOinfo[1],IOinfo[2])
    sample_prop = find_total_genotype(IOinfo[3],dict_list)
    print(sample_prop)
    datamtx = props_to_guess(sample_prop)

    outstring = '''
    %s \t deviation of : %.3f\n
    %s \t deviation of : %.3f\n
    %s \t deviation of : %.3f\n
    %s \t deviation of : %.3f\n
    %s \t deviation of : %.3f\n
    ''' % (datamtx[0][0],datamtx[1][0],
           datamtx[0][1],datamtx[1][1],
           datamtx[0][2],datamtx[1][2],
           datamtx[0][3],datamtx[1][3],
           datamtx[0][4],datamtx[1][4]
           )
    print(outstring)


if __name__ == '__main__':
    main()