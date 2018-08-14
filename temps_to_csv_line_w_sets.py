
# define function to do the things it says in the help use -h to view
def inpath():
    import sys
    inputa=''
    inputb=''
    inputs='unexpected line in info line please check manually!!!'
    insert_file=''
    for f in range(len(sys.argv)):
        if sys.argv[f] in ['-h','--help']:
            sys.stderr.write('Usage:\ttemps_to_csv_line [flags] >> [output.csv]\n'
                       '-h\t--h\tPrints help\n'
                       '-ia\t--inputa\t should be followed with the input file\n'
                       '-ib\t--inputb\t should be followed with the input file\n'
                       '-iq\t--inputfastq\tshould be followed with an input fastq \n'
                       '-ii\t --input_insert \t should be followed by the file with the insert information')
        elif sys.argv[f] in ['-ia', '--inputa']:
            inputa = sys.argv[f + 1]
        elif sys.argv[f] in ['-ib', '--inputb']:
            inputb = sys.argv[f + 1]
        elif sys.argv[f] in ['-iq','--inputfastq']:
            inputs = sys.argv[f+1]
        elif sys.argv[f] in ['-ii','--input_insert']:
            insert_file = sys.argv[f+1]
    return [inputa,inputb,inputs,insert_file]

def main():
    # inport paths to files from cmand line
    inpaths = inpath()

    ## inpaths = ['/home/sheepless/Desktop/BD2K/SandBox/ziped_fastq_files/SI-GA-C3-1_S25_L001_R1_001.infoa',
    ##           '/home/sheepless/Desktop/BD2K/SandBox/ziped_fastq_files/SI-GA-C3-1_S25_L001_R1_001.infob',
    ##           '/home/sheepless/Desktop/BD2K/SandBox/ziped_fastq_files/SI-GA-C3-1_S25_L001_R1_001_sorted.sam']
    # open the files
    infilea = open(inpaths[0], 'r')
    infileb = open(inpaths[1], 'r')
    infiles = open(inpaths[2], 'r')
    infilei = open(inpaths[3], 'r')
    # initiolize varyables
    smpl = ''
    numr = 0.0
    numd = 0.0
    numq = 0.0
    numm = 0.0
    prem = 0.0
    pred = 0.0
    preq = 0.0
    b_codes = set()
    inserts = ''
    # get sample name from path
    path_list = inpaths[0].split('/')
    smpl = path_list[len(path_list)-1].split('.')[0]
    histpath = inpaths[0].rstrip('.infoa') + '_hist.pdf'

    #deal with flagstat file with all the reads using the first word in the line to decide what num is there
    for line in infilea:
        linelist = line.rstrip('\n').split('+')
        ##print (linelist[1].split(' ')[2])
        if linelist[1].split(' ')[2] == 'in':
            # assign value of totle number of reads
            numr = float(linelist[0])
        elif linelist[1].split(' ')[2] == 'duplicates':
            # assign value of number of duplicate reads
            numd = float(linelist[0])
        elif linelist[1].split(' ')[2] == 'mapped':
            # assign number of reads that map to the reference sequence
            numm = float(linelist[0])

    # deal with file with only 1q21.1 reads
    for line in infileb:
        linelist = line.rstrip('\n').split('+')
        if linelist[1].split(' ')[2] == 'in':
            #fig number of reads that map to the regins defined in the shell script
            numq = float(linelist[0])

    # deal with the input fastq file (this is the longest step
    seq_line = False
    for line in infiles:
        # see if the flag that tells us if the last line was the line before the sequence
        if seq_line:
            # if it is check the barcode agenst seen barcodes
            b_codes.add(line[0:15])
            # then wether or not you saw it the next line will not be sequence so set that flag to false
            seq_line = False
        # if it isn't a sequence check if the next line is annd set the flag if you ned to
        elif line[0] == '@':
            seq_line =True
    # at the end see how many keys there are in our lst
    unkr = float(len(b_codes))
    # look at the insert info file and set an indexing varyable to 0
    i_index = 0
    for line in infilei:
        # loop through the file , and if it is line 8
        i_index = i_index + 1
        if i_index == 8:
            # see if the file is as we expect if so xet the varyable
            if line.rstrip('\n').split('\t')[8] == 'FR':
                inserts = str((line.rstrip('\n').split('\t')[0])) + u'\u00B1' + str((line.rstrip('\n').split('\t')[2]))



    #calculate precentages from counts obtained
    prem = 100 * (numm/numr)
    pred = 100 * (numd/numr)
    preq = 100 * (numq/numr)
    read_per_unk = numr/unkr

    print ('%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%s,%s' %(smpl,numr,numd,numq,preq,prem,pred,unkr,read_per_unk,inserts,histpath))
    # close the files we opened
    infilea.close()
    infileb.close()
    infiles.close()
    infilei.close()



if __name__ == '__main__':
    main()