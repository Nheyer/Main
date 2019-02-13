#!usr/bin/env python3



# this assumes you are runing it in a directory with the reference genome index with prefix ref.fasta

def paths():
    # initialise and set defaults
    import sys
    # the gate variables are boolions that let the script know if the current eliment in the sys.argv is part of a given list
    input_gate = False
    memopt_gate = False
    genes_gate = False
    # these are just variable that need to be defined a their respective nulls
    inputs = []
    outputfilenm = ''
    memopts = ''
    sam_threds = 4
    # these are default values more info on them in the main loop  of the method
    qk = False
    qual_flag = 10
    genes_list = ['119911553,120069703,N,chr1:',
                  '146151908,146229032,A,chr1:',
                  '148602850,148679774,B,chr1:',
                  '149390591,149468389,C,chr1:',
                  '120723917,120798050,D,chr1:']
    splkey = 14


    for n in range(len(sys.argv)):
        # if they use these flags print help
        if sys.argv[n] in ['-h','--help']:
            sys.stderr.write('''
            ~~~The script assumes the reference index is in the same directory you are runing this from and has the prefix ref.fasta ~~~\n
            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~updated for samtools1.9.0~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n
            'Usage: \t bin_clasifier -i/--Inputfas [file paths] -o/--OutputFile [filename.csv] \n
            -----------------------------------------------------------------------------------------\n
            Manditory flags:\n
            -i \t --Inputfas \t a space seperated list of paths to the fasta or fastq files in the bins\n
            -o \t --Outputfile \t the path including the filename of the csv you would like the clasifications in\n
            ------------------------------------------------------------------------------------------\n
            Optional flags:\n
            -h \t --help \t this will print the usage statement to the standard error \n
            -q \t --QuickKill \t if this flag is used it will remove all sam files as quickly as possable\n
            -b \t --memOptions \t a space seperated list of flags to pass to mem\n
            -Q \t --Quality \t minimum read quality score to consider default == 10\n
            -s \t --splits \t how much to split the genes [1 to 30] default == 14\n
            -g \t --genes \t a space seperated list of start,end,name,chr[int]: default == NOTCH2NL\n
            
            ''')

            # exit if you print help
            sys.exit()
        # if you see the -i flag the next elimints till the next flag are input files so open the input gate
        # make sure to close all the other gates as well becouse this is where the lists stop and a new one starts
        elif sys.argv[n] in ['-i','--Inputfas']:
            input_gate =True
            memopt_gate = False
            genes_gate = False
        # set output if you see the flag
        elif sys.argv[n] in ['-o','--Outputfile']:
            outputfilenm = sys.argv[n+1]
            input_gate = False
            memopt_gate = False
            genes_gate = False
        # input a quality threshold from the user if not default
        elif sys.argv[n] in ['-Q', ' --Quality']:
            qual_flag = sys.argv[n+1]
            input_gate = False
            memopt_gate = False
            genes_gate = False
        # import the number of times the gene should be partitioned if not default
        elif sys.argv[n] in ['-s','--splits']:
            splkey = sys.argv[n+1]
            input_gate = False
            memopt_gate = False
            genes_gate = False
        # if you see the flag the user wants you to kill the sam files as fast as possible so set qk to true
        elif sys.argv[n] in ['-q','--QuickKill']:
            qk = True
            input_gate = False
            memopt_gate = False
            genes_gate = False
        # if you see this flag the next eliments till a defined flag are flaags for mem (there is no overlap in flags)
        #open the mem option gate
        elif sys.argv[n] in ['-b','--memOptions']:
            input_gate = False
            memopt_gate = True
            genes_gate = False
        # see if we need to edit the gene list  if so cleart the list then start up the workhorse below
        elif sys.argv[n] in ['-g','-genes']:
            genes_list = []
            genes_gate = True
            input_gate = False
            memopt_gate = False
	elif sys.argv[n] in ['-R','--reference']:
            path_to_fa = sys.argv[n+1]
            genes_gate = False
            input_gate = False
            genes_gate = False
        elif sys.argv[n] == '-t' and not memopt_gate:
            sam_threds = sys.argv[n+1]

        # workhorses for the flags that input lists (essentially adds stuff to python lists or strings as needed)
        elif input_gate :
            inputs.append(sys.argv[n])
        elif memopt_gate :
            memopts = memopts + ' ' + sys.argv[n]
        elif genes_gate:
            genes_list.append(sys.argv[n])
    # we have now delt with the io return all the things we need
    return [inputs,outputfilenm,qk,memopts,str(qual_flag),splkey,genes_list,path_to_fa,str(sam_threds)]

def inputconf(list, split):
    # set the spliting factor defined as the defult or user input
    SplitIntoParts = split
    # create the dictionary to be added to, with a none key already inside to deal with situations that
    # result from no reeds in any region being tested
    diction = {'none':'**'}


    # a bunch of empty lists to be filled based off the splits needed
    genes = []
    llist =[]
    for dummy in range(int(SplitIntoParts)):
        llist.append([])
    # loop through the input list defined by the user or default
    for eli in list:
        genes.append(eli.split(','))
    # loop through all the gene regions we found in the list
    for gene in genes:
        # see if the user made an error I can easaly see happening and fix it...
        if SplitIntoParts == 0:
            dif = float(int(gene[1]) - int(gene[0]))
        # set the difference between regions as the diefference between the beginning and end devided by the number of parts
        else:
            dif = float(int(gene[1]) - int(gene[0]))/float(SplitIntoParts)
        # set the running varyable to the start of the region
        runing = int(gene[0])
        # use the outer loop to construct parts needed to add a key to the dictionary then do that
        for i in range(int(SplitIntoParts)):
            # combine the chromosome input with the running start (same as the end of the last iteration)
            strn1 = gene[3] + str(int(round(runing))) + "-"
            # update the running var to the end of the current location using the calculated dif
            runing = runing + dif
            # create the second part of the key term
            strn1 = strn1 + str(int(round(runing)))
            # find the key value needed
            strn2 = "|" + gene[2] + str(i)
            # combine the strings so they can be split properly to be added to the dictionary
            strn = strn1 + strn2
            # finally add them to the dictionary
            diction[strn.split('|')[0]] = strn.split('|')[1]
            # inner loop adds the key to the correct sub list if it is add the key to the list
            for j in range(len (llist)):
                if i == j:
                    llist[j].append(strn1)
    # return a dictionary, and a list of the sublists
    return [diction,llist]


def main():

    # import needed packages
    import os
    import subprocess
    import sys




    # get needed paths from the shell and open the output file
    io_list = paths()
    input_list = io_list[0]
    ouput_path = io_list[1]
    # create the gene lists and dictionary we need based on defaults or user flags
    dictionAndlists = inputconf(io_list[6],io_list[5])
    outfile = open(ouput_path, 'w')
    # see how many lists are in the list of lists of gene regions and set our first inner loop to loop that may times
    strloop = len(dictionAndlists[1])
    reference_file_fa = ' '+ io_list[7] + ' '


    # start the main loop to loop through input files (bins)
    for loopfile in input_list:
        # write the file name to the line then add a comma
        outfile.write("%s," % (loopfile))
        # rum bwa with the flags the user defined  referencing the genmome called ref.fasta and stores the sam to be used
        # the sam file is marked for the cleanup later , then the sam file is redused to only quality 60 reads and
        # converted to a bam file that is also marked for cleanup
        print(io_list[3])
        os.system('bwa mem '+io_list[3]+ reference_file_fa +loopfile+' > tobekilledfoo_quick_working_sequence.sam')
        os.system('samtools view -@ '+io_list[8]+' -Sbq'+io_list[4]+' tobekilledfoo_quick_working_sequence.sam > tobekilledfoo_working_sequence.bam')
        # if the user flaged quick cleanup of sam files remove the sam file now, if not we can do it at the end
        if io_list[2]:
            os.system('rm tobekilledfoo_quick_working_sequence.sam')
        # sort and index the bam file
        os.system('samtools sort tobekilledfoo_working_sequence.bam > tobekilledfoo_sort_seq.bam')
        os.system('samtools index tobekilledfoo_sort_seq.bam')
        # see explonation for this above
        for r in range(strloop):
            # reseting things to zero
            read_gate = 0
            best_gene = 'none'
            # if the current sublist had no eliments ie is an empty list break the loop so the code dose not and none keys till 30.
            if dictionAndlists[1][r] == []:
                break
            # loop through all the gene regions ( NOTCH2NL by default) and create temperary bams including only reads from that region and give a cout as an int
            for gene in dictionAndlists[1][r]:
                os.system('samtools view -b tobekilledfoo_sort_seq.bam '+gene+'> tobekilledfoo_gene_limed.bam')
                # get the flagstat info from samtools  and use string manipulation to find the number of reads in the bam
                flgst = subprocess.check_output(['samtools' , 'flagstat' , 'tobekilledfoo_gene_limed.bam'])
                cur_read = (int(str(flgst).split('+')[0].lstrip("b'")))
                # if that count is greater then the previous largest set this region to be the best gene, and update greatist read
                if cur_read > read_gate:
                    best_gene = gene
                    read_gate = cur_read
            # write the best matched region to the user defined output file
            outfile.write('%s' % (dictionAndlists[0][best_gene]))
        outfile.write("\n")
    # cleanup marked files we will also tell the user here
    sys.stderr.write("Cleaning up ...\n")
    os.system("rm -f tobekilledfoo_*")
    outfile.close()











if __name__ == '__main__':
    main()
