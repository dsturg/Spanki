
# The fasta file downloaded from GENCODE is not compatible with Spanki,
# because the pyfasta package used in Spanki require ``>chr1`` rather than 
# ``>chr1 extra information``. Thus, the fasta file needs a bit change, which 
# can be done with this program.


from optparse import OptionParser

if __name__ == "__main__":
    parser = OptionParser()

    parser.add_option("--ref_file", "-f", dest="ref_file", default=None,
        help="The reference file in fasta format.")

    (options, args) = parser.parse_args()
    fid = open(options.ref_file, 'r')
    fasta = fid.readlines()
    fid.close()
    
    fasta_out = ".".join(options.ref_file.split(".")[:-1])
    fid = open(fasta_out + ".spanki.fa", "w")

    for i in range(len(fasta)):
        if fasta[i][0] != ">":
            fid.writelines(fasta[i])
        else:
            fid.writelines(fasta[i].split()[0] + "\n")
    fid.close()