def main():
    # import some packages
    import subprocess
    import numpy
    import sys
    #  define some vars
    in_dir, t, out_file , check_unsplit, name , field = inputs()
    indiv_nanme = ''
    # open the file to write to
    outfile = open(out_file,'w')
    # write a header line
    outfile.write("Sample,A+B,A+B (A154T),A+B (T197I),A+B (M>I),A+B (R113*),C,C ( 2 ntdelletion),N,D\n")
    # see the sub directories thatg will become our individual names
    files = subprocess.check_output(['ls',in_dir]).rstrip('\n').split('\n')
    # use a list comprehension to fix the paths
    I_paths = [in_dir + x for x in files]
    # look through the directories seen
    for individual in I_paths:
        # initiolize a np aray to store counts of A/B paratype, A/B with A154T, A/B with T197I, A/B I>M, A/B R113*, C paratype, C dell , N para , D para
        indiv_type = numpy.array([0,0,0,0,0,0,0,0,0])
        # get a list of files in the directory ( assuming they are bgziped and indexed using tabix)
        vcfs_and_index = subprocess.check_output(['ls',individual+'/']).rstrip('\n').split('\n')
        vcfs = []
        # remove the tabix files from the list of vcfs to parse
        for nm in vcfs_and_index:
            if nm.split('.')[len(nm.split('.'))-1] == 'gz':
                vcfs.append(nm)
        # loop thorugh the vcf files
        for vcf_path in vcfs:
            mult = 1
            real_path = individual + '/' + vcf_path
            if check_unsplit and vcf_path.split("_")[field + 1][0].isdigit():
                mult = int(vcf_path.split("_")[field + 1][0])
            elif len(vcf_path.split("_")) > (field + 2):
                 if check_unsplit and vcf_path.split('_')[field + 2][0].isdigit():
                     mult = int(vcf_path.split('_')[field + 2][0])
            print (real_path, "\t" ,mult)
            # add this vcf's paratype information to the individuals a\overall vector
            indiv_type = indiv_type + mult * clasify_muti_type(real_path,field)
        if name == 'long':
            indiv_nanme = individual
        elif name == 'normal':
            indiv_nanme = individual.split('/')[len(individual.split('/'))-1]
        elif name == 'short':
            indiv_nanme = individual.split('/')[len(individual.split('/'))-1].split('_')[0]
        # print out the final csv line for the individual
        outfile.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (indiv_nanme,indiv_type[0],indiv_type[1],indiv_type[2],indiv_type[3],indiv_type[4],indiv_type[5],indiv_type[6],indiv_type[7],indiv_type[8]))




def clasify_muti_type(vcf_path , field):
    # import some packages
    import pysam
    import numpy
    ptype_array = numpy.array([0,0,0,0,0,0,0,0,0])
    # get the paratype as defined by the first character of the vcf file name
    path_list = vcf_path.split('/')

    paratype = path_list[len(path_list)-1].split('_')[field][0]
    # print the paratype for user conformation
    print (paratype)
    # open the vcf
    vcf = pysam.VariantFile(vcf_path)
    # if it is N then just set the vector accordingly
    if paratype == 'N':
        ptype_array = numpy.array([0, 0, 0, 0, 0, 0, 0, 1, 0])
    # if it is D then just set the vector accordingly
    elif paratype == 'D':
        ptype_array = numpy.array([0, 0, 0, 0, 0, 0, 0, 0, 1])
    # if it is A or B then cheack to see if any of the b/A typ mutations are present, and change the vector accordingly
    elif paratype in ['A', 'B']:
        ptype_array = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0])
        for val in (vcf.fetch('NOTCH2NL-consensus', 127050, 200127)):
            if val.pos == 199997:
                if 'A' in val.alts:
                    ptype_array = ptype_array + numpy.array([0, 1, 0, 0, 0, 0, 0, 0, 0])
            elif val.pos == 200127:
                if 'T' in val.alts:
                    ptype_array = ptype_array + numpy.array([0, 0, 1, 0, 0, 0, 0, 0, 0])
            elif val.pos == 127052:
                if 'A' in val.alts:
                    ptype_array = ptype_array + numpy.array([0, 0, 0, 1, 0, 0, 0, 0, 0])
            elif val.pos == 191863:
                if 'T' in val.alts:
                    ptype_array = ptype_array + numpy.array([0, 0, 0, 0, 1, 0, 0, 0, 0])
    # checks to see if it is a C paratype...
    elif paratype == 'C':
        ptype_array = numpy.array([0, 0, 0, 0, 0, 1, 0, 0, 0])
        for val in (vcf.fetch('NOTCH2NL-consensus',127060,127066)):
            if val.pos == 127063:
                # cheack to see if the vcf denotes a deletion....
                if len(val.ref) - numpy.average(numpy.array([len(x) for x in val.alts])) == 2.0:
                    ptype_array = ptype_array + numpy.array([0, 0, 0, 0, 0, 0, 1, 0, 0])

    return ptype_array




def inputs():
    import sys
    # set defaults
    input_dir = '../VCfs'
    out_file = 'OUTS.csv'
    cores = 30
    unsplit_in = False
    samp_nm = 'normal'
    f = 0
    # define my usage statement
    usage = '''
    python gen_csv.py -i [input_dir] \n  
    Manditory Flags:\n
    -i\t--input_dir\t the path to the directory contaning directories for each individual and in those directory have vcfs with the first character being the paratype.\n
    Optional Flags:\n
    -o\t--out_file\t Print the output here instead of at OUTS.csv\n
    -@\t--cores    \t Set the maximum number of cores to use ( not currently threaded) \n
    -u\t--UNSPLIT\t use when there is a VCF that represents 2 or more reads and should be multiplied file should be named "partype_#of groups it should represent_other info for user.vcf.gz"    
    -l\t--longpath\tuse the full path to the sample directory instead of just the directory name (use -sh to make it only the charicteres before the first '_'
    
    '''
    # parse user input
    for arg_index in range(len(sys.argv)):
        if sys.argv[arg_index] in ['--help','-h']:
            sys.stderr.write(usage)
            sys.exit()
        elif sys.argv[arg_index] in ['-i','--input_dir']:
            input_dir = sys.argv[arg_index + 1]
        elif sys.argv[arg_index] in ['-@','--cores']:
            cores = int(sys.argv[arg_index +1])
        elif sys.argv[arg_index] in ['-o','--out_file']:
            out_file = sys.argv[arg_index+1]
        elif sys.argv[arg_index] in ['--UNSPLIT', '-u'] :
            unsplit_in = True
        elif sys.argv[arg_index] in ['-l','--longpath']:
            samp_nm = 'long'
        elif sys.argv[arg_index] == '-sh':
            samp_nm = 'short'
        elif sys.argv[arg_index] in ['-f', '--field']:
            f = int(sys.argv[arg_index + 1])

    # tell the user what we saw
    sys.stderr.write("Taking inputs from %s \n Using %s threads (no current threading so this dosn't mean anything...) \n Outputing to %s\n"    % (input_dir,cores,out_file))

    return input_dir , cores ,out_file , unsplit_in, samp_nm, f



if __name__ == '__main__':
    # if you are not importing the function run it
    main()
