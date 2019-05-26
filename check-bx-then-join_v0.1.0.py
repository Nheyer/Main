import sys
import pysam
import os

def inputs():
    bam_2 = ''
    bam_1 = ''
    threads = 10
    tag_temp = 'BX'
    txt_out = ''
    file_out = ''
    usage_stmnt='''
    Usage: pyhton check_BX_for_dups.py -i/--inputs <the two bams you are checking> -@/--threads <max number of threads>\n
    ------------------------------------------------------------------------------------------------------------------
    Mandatory flags:
    -i\t --inputs\t The two .bam files from longranger you want to merge.
    ------------------------------------------------------------------------------------------------------------------
    Optional flags:
    -@\t --threads\t The maximum number of threads to use default [10]
    -L\t --logs\t Put important print outs here default [/dev/stdout]
    -h\t --help\t Print this information 
    -o\t --out\t print to givin file instead of modified input file
    -T\t --tag\t What tag to check for collisions [BX]\n
    '''
    for id in range(len(sys.argv)):
        if sys.argv[id] in ['-h','--help']:
            sys.stderr.write(usage_stmnt)
            sys.exit()
        elif sys.argv[id] in ['-i', '--inputs']:
            bam_1 = sys.argv[id + 1]
            bam_2 = sys.argv[id + 2]
        elif sys.argv[id] in ['-@','--threads']:
            threads = sys.argv[id + 1]
        elif sys.argv[id] in ['-L','--logs']:
            txt_out = sys.argv[id + 1]
        elif sys.argv[id] in ['-o', '--outs']:
            file_out = sys.argv[id+1]
        elif sys.argv[id] in ['-T', '--tags']:
            tag_temp = sys.argv[id+1]
    return (tag_temp,bam_1, bam_2, threads, txt_out, file_out)
def make_BX_set(bam, TAG):
    barcodes = set()
    #open the bam file
    var_bam = pysam.AlignmentFile(bam)
    # loop through all the lines in the file and add the barcode to a set
    for contig in var_bam:
        barcodes.add(contig.get_tag(TAG))
    # return a set containing all of the barcodes in the bam file
    var_bam.close()
    return barcodes

def see_dupes(bam_1,bam_2,out_s,TAG):
    dupes = False
    outline = 'If this results all is fine in %s and %s ' % (bam_1,bam_2)
    pr_sd = False
    if out_s == '':
        pr_sd = True
    # create sets for each bam
    bx1_set = make_BX_set(bam_1,TAG)
    bx2_set = make_BX_set(bam_2,TAG)
    # make a set of all barcodes
    unique_bx = bx1_set.union(bx2_set)
    # find total barcodes in each file and add the numbers
    total_num = len(bx1_set) + len(bx2_set)
    # fine the number of unique barcodes in each file
    unique_bx_num = len(unique_bx)
    # see if they are different if so tell us
    if total_num != unique_bx_num:
        outline = ("There are %s duplicates, trying to fix" % (total_num - unique_bx_num))
        dupes = True

    if pr_sd:
        print(outline)
    else:
        out_file = open(out_s, 'w')
        out_file.write(outline)
        out_file.close()
    return dupes

def mk_unque(bam_file,unique_str,TAG):
    new_bam = bam_file.rstrip('.bam') + 'unique-'+TAG+'.bam'
    # open the bam that needs bx changes, and one to write to
    bam = pysam.AlignmentFile(bam_file, 'rb')
    bam_n = pysam.AlignmentFile(new_bam, 'wb', template=bam)
    # loop through all the lines
    for line in bam:
        n = 0
        curr_tags = line.get_tags()
        # update the tags to include the new unique barcode
        for tag in curr_tags:
            if tag[0] == TAG:
                curr_tags[n] = (TAG, curr_tags[n][1]+unique_str)
            n = n + 1
        line.set_tags(curr_tags)
        # write this line to a new file
        bam_n.write(line)
    bam.close()
    return new_bam



def main():
    # get user inputs
    tg,b1,b2,t,out,npath = inputs()

    if see_dupes(b1,b2,out,tg):
        new_path_1 = mk_unque(b1, '-U-1', tg)
        new_path_2 = mk_unque(b2, '-U-2', tg)
        if see_dupes(new_path_1,new_path_2,out):
            sys.exit()
        else:
            b1 = new_path_1
            b2 = new_path_2
            print ('here')
    if npath == "" :  # if the user didn't give us an output path construct one
        npath = b1.split('/')[-1].split('.')[0] +'-con-combined.bam'
    join = 'samtools merge -@' + str(t) + ' -f ' + npath + ' ' + b1 + ' ' + b2
    print (join)
    os.system(join)

if __name__ == '__main__':
    main()


