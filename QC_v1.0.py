import os
import sys
import glob
import subprocess
DEBUG = True


# define function to do the things it says in the help use -h to view
def inpath():
    inputs = []
    ref_path = ""
    no_buf =  'chr1:119911553-120069703 chr1:146151908-146229032 chr1:148602850-148679774 chr1:149390591-149468389 chr1:120723917-120798050'
    mid_buf = 'chr1:119989248-120190000 chr1:149328818-149471561 chr1:148600079-148801427 chr1:146149145-146328264 chr1:120705669-120801220'
    big_buf = 'chr1:143200001-150600000 chr1:120400001-121700000'
    io_gate = False
    region = no_buf
    threads = 1
    picard = "~/bin/picard.jar"
    out = "QC.csv"
    help_str = "Usage:\n"
    help_str += "python QC_vX.X.py -i [paths] -R [path] -o [path]\n"
    help_str += "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    help_str += "-h\t--help\t\t Print this help message\n"
    help_str += "-i\t--input\t[paths]\t A space separated list of directeries containing the fastqs returned from sequencing\n"
    help_str += "-o\t--output\t[path]\t csv to write to\n"
    help_str += "-R\t--reference\t[path]\t A path to the indexed reference genome\n"
    help_str += "-p\t--picard-jar\t[path]\t A path to the installed copy of picard default assumes it is in ~/bin/picard.jar\n"
    help_str += "-@\t--threads\t[int]\t Maximum number of threads to use\n"
    help_str += "-b\t--buffer\t[mid/big/none]\t what size of buffer around N2NL should we use default in no,\n\t\t\t you may also use samtools like regions in a comma seperated list\n "

    for i in range(len(sys.argv)):
        if sys.argv[i] in ['-h','--help']:
            sys.stderr.write(help_str)
            sys.exit(-1)
        elif sys.argv[i] in ['-i', '--input']:
            io_gate = True
        elif sys.argv[i] in ['-R', '--reference']:
            ref_path = sys.argv[i+1]
            io_gate = False
        elif sys.argv[i] in ['-b', '--buffer']:
            if sys.argv[i + 1] == "mid":
                region = mid_buf
            elif sys.argv[i + 1] == "big":
                region = big_buf
            elif sys.argv[i + 1] == "none":
                region = no_buf
            else:
                sys.stderr.write("Treating -b/--buffer input as a comma seperated list of samtools like gene regions\n")
                region = sys.argv[i + 1].replace(",", " ")
            sys.stderr.write("Using " + region + " as the region to use for calculating enrichment.\n")
            io_gate = False
        elif sys.argv[i] in ['-@', '--threads']:
            threads = sys.argv[i + 1]
            io_gate = False
        elif sys.argv[i] in ['-p', '--picard-jar']:
            picard = sys.argv[i + 1]
            io_gate = False
        elif sys.argv[i] in ['-o', '--output']:
            out = sys.argv[i + 1]
            io_gate = False
        elif io_gate:
            inputs.append(sys.argv[i])
    return inputs, ref_path, region, threads, picard, out


def get_bams(r1ANDr2, threads, dir, ref, region):
    print(r1ANDr2, threads, dir, ref, region)
    # allign to ref, then sort and index
    subprocess.call("bgzip " + "-d -@ " + threads + " " + r1ANDr2[0], shell=True)
    subprocess.call("bgzip " + "-d -@ " + threads + " " + r1ANDr2[1], shell=True)
    cmd = "bwa mem -t " + threads + " -o " + dir + "/alligned.sam " + ref + " " + r1ANDr2[0].rstrip(".gz") + " " + r1ANDr2[1].rstrip(".gz")
    subprocess.call(cmd, shell=True)
    subprocess.call("samtools view -@ " + threads + " -Sb -o " + dir + "/alligned.bam " + dir + "/alligned.sam ", shell=True)
    subprocess.call("samtools sort -o " + dir + "/sorted.bam " + dir + "/alligned.bam ", shell=True)
    subprocess.call("samtools index " + dir + "/sorted.bam ", shell=True)
    # limit to enriched region, and re-sort/ re-index
    subprocess.call("samtools view -b -@ " + threads + " -o" + dir + "/limited.bam " + dir + "/sorted.bam " + region, shell=True)
    subprocess.call("samtools sort -o " + dir + "/lim-sorted.bam " + dir + "/limited.bam ", shell=True)
    subprocess.call("samtools index " + dir + "/lim-sorted.bam ", shell=True)
    raw_check = glob.glob(dir + "/*.bai")
    if len(raw_check) == 2:
        return raw_check[0].rstrip(".bai"), raw_check[1].rstrip(".bai")
    else:
        sys.stderr.write("failed to process files in dir:\t" + dir + "\n")
        sys.exit(-1)
def flagstat_list_to_vals(Flagstat_list):
    if(DEBUG):
        for i in range(len(Flagstat_list)):
            print(Flagstat_list[i], i)
    num_reads  = int(Flagstat_list[0].lstrip("'")) + int(Flagstat_list[2])
    num_maped  = int(Flagstat_list[18].lstrip("duplicates\\n")) + int(Flagstat_list[20])
    num_dups   = int(Flagstat_list[15].lstrip("supplementary\\n")) + int(Flagstat_list[17])
    prop_maped = Flagstat_list[22].lstrip("(")
    if (DEBUG):
        print(float(num_maped)/float(num_reads), " should equal ", prop_maped)
    return num_reads, num_maped, num_dups, prop_maped


def get_insert_stats(path):
    insert_file = open(path, 'r')
    mean = 0.0
    std = 0.0
    found_header = False
    for line in insert_file:
        if(line[0] != "#" and len(line) > 10):
            line_list = line.rstrip("\n").split("\t")
            if(DEBUG):
                print(line_list)
            if(line_list[5] == "MEAN_INSERT_SIZE"):
                found_header = True
            elif(found_header):
                mean = float(line_list[5])
                std = float(line_list[6])
                break
    insert_file.close()
    return str(mean) + u'\u00B1' + str(std)


def get_num_BX(path):
    import pysam
    bxs = set([])
    bam = pysam.AlignmentFile(path)
    for read in bam:
        try:
            bxs.add(read.get_tag('BX'))
        except:
            continue
    return len(bxs)


def main():
    # inport paths to files from cmand line
    inpaths, ref, reg, t, picard, out= inpath()
    # open file to write to
    outfile = open(out, 'w')
    outfile.write("Sample,reads,duplicates,1q21.1,1q21.1 Enrichment,% mapped,%duplicates,unique BCs,reads/unique BC,unique_1q12_BX,1q21 reads/unique BX,average insert size,path to histogram")
    for dir in inpaths:

        fqs = glob.glob(dir + "/*R1*")
        fqs.append(glob.glob(dir + "/*R2*")[0])
        full_bam, limited_bam = get_bams(fqs, str(t) , dir, ref, reg)
        full_flagstat = str(subprocess.check_output(["samtools", "flagstat", full_bam])).lstrip('b').split(" ")
        lim_flagstat  = str(subprocess.check_output(["samtools", "flagstat", limited_bam])).lstrip('b').split(" ")
        subprocess.call("java -jar ~/bin/picard.jar CollectInsertSizeMetrics I=" + full_bam +
                        " O="+full_bam.rstrip(".bam") + "inserts.txt H=" + full_bam.rstrip(".bam") + ".hist.pdf", shell=True)
        t_num_reads, t_num_maped, t_num_dup, t_prop_maped = flagstat_list_to_vals(full_flagstat)
        t_num_bx = get_num_BX(full_bam)
        l_num_reads, l_num_maped, l_num_dup, l_prop_maped = flagstat_list_to_vals(lim_flagstat)
        l_num_bx = get_num_BX(limited_bam)
        ins_size = get_insert_stats(full_bam.rstrip('.bam') + "inserts.txt")
        if(t_num_bx > 0):
            l_prop_bx = float(l_num_reads)/float(l_num_bx)
        else:
            l_prop_bx = -1
        if(t_num_bx > 0):
            t_prop_bx = float(t_num_reads)/float(t_num_bx)
        else:
            t_prop_bx = -1

        print(dir.split("/")[-1], t_num_reads, t_num_dup, l_num_reads, t_prop_maped,
              float(t_num_dup)/float(t_num_reads), float(l_num_reads)/float(t_num_reads),
              t_num_bx, t_prop_bx, l_num_bx, l_prop_bx, ins_size)



    outfile.close()




if __name__ == '__main__':
    main()