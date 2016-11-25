"""
LRGparser.py

created: November 2016

Tested on python versions 3.5.2 and 

@authors: Laura Carreto, Rosie Coates-Brown

usage: python LRGparser.py -g [LRG file name] -d [True/False] -i [True/False]

required parameters:
-g, --gene	[name of LRG file without .xml suffix]

optional parameters:
-h, --help shows this message and quits
-d, --diff, --difference [True/False] triggers or supresses the output of [LRG]_diffs.csv
-a, --annotations = [True/False] triggers or suppresses [LRG]_annotation.csv

output:

[LRG]_t1.bed: a tab separated bed file containing the chromosome number, exon start position, exon end position  
[LRG]_diffs.csv: a csv file containing the differences between 37 and 38, False will suppress this
[LRG]_annotation.csv: a csv file of gene information including synonyms, lsdb, long gene name

"""

try:
    import xml.etree.ElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import sys, os, csv, getopt


def read_file(genein):
    """
    read in the LRG.xml file
    test: does the file exist

    parameters: genein(string returned from command line arguments)
    returns: root(ElementTree root node), gene(user input string variable)

    """
    gene = genein
    
    file_name = gene+'.xml'
    #file_path = '/home/swc/Desktop/LRGParser/'
    file_path = '/Users/rosiecoates/Documents/Clinical_bioinformatics_MSc/programming/assignment/'
    full_path = file_path+file_name
    
    # check if required xml file exists in the directory defined in file_path
    try:
        tree = ET.parse(full_path)
    except:
        print("couldn't open file... check you have supplied an LRG file name without extension, and file is an XML")
        usage()
        sys.exit(2)
        
    # parse xml and create root object    
    tree = ET.ElementTree(file=full_path)
    root = tree.getroot()

    return root, gene
    
def write_csv(mylist, myfilename, mode):
    """
    Output to CSV file from a list
    
    Parameters: mylist (list), myfilename(string), mode (string); the mode options used in the code are: 'a'= append, 'w'=write
    
    """
    #
    out = csv.writer(open(myfilename, mode), quoting=csv.QUOTE_ALL)
    out.writerow(mylist)
        
    return


def bed_file(root, gene):
    """
    produces a bed file containing the chromosome, genomic start position and genomic end position

    test: is start position bigger than end position?
    test: is the calculate chromosome position between the given chromosome start and end?

    parameters: root(ElementTree root node), gene(user input string variable)

    returns: exon_ranges(dict)
    """
    # initialise dictionary of exon ranges; format: exon = list(exon_start, exon_end)    
    exon_ranges={}

    for mapping in root.findall("./updatable_annotation/annotation_set[@type='lrg']/mapping[@type='main_assembly']"):
        # for each branch, get start and end coordinates in the reference genome and chromosome number 
    
        chr = mapping.attrib['other_name'] # chromosome
        #print (chr)
        ref_start = int(mapping.attrib['other_start']) # start (converted to integer)
        #print (ref_start)
        ref_end = int(mapping.attrib['other_end']) # end (converted to integer)
        #print(ref_end)
        
        # get strand to account for + or - strand when converting exon lrg coordinates to genomic coordinates
        for mapping_span in mapping:  
            strand = mapping_span.attrib['strand']
            print (gene+' is coded in strand '+ strand)
        
        if strand == '1':
            offset = (ref_start) - 1 # offset start to compensate for string count
            offset_int = int(offset) # convert offset start to integer
            #print (offset)
            
        elif strand == '-1':#WORK IN PROGRESS; no difference at the moment for + or - strand
            offset = (ref_start) - 1 # modify 
            offset_int = int(offset) # modify
            #print (offset)


    for id in root.findall("./fixed_annotation/id"):
        id_tag = id.text

    for transcript in root.findall("./fixed_annotation/transcript"):
        # for each transcript, get exon coordinates into bed file
        trans_num =  (transcript.attrib['name']) # transcript number
        ref_name = gene + "_" + trans_num + ".bed" # define bed file name
        
        bedfile = open(ref_name, 'w') # open bed file
        
        for exon in transcript:
            #for each exon, get start and end coordinates as integers; offset lrg coordinates to get genomic coordinates 
            if exon.tag == 'exon':
                exon_number =  (exon.attrib['label'])

                for coord_sys in exon:
                    if (coord_sys.attrib['coord_system']) == id_tag:
                        start=(coord_sys.attrib['start'])
                        start_int = int(start)
                        end = coord_sys.attrib['end']
                        end_int = int(end)

                        gen_st = str(start_int + offset_int)
                        gen_st_int = int(gen_st)
                        gen_end = str(end_int + offset_int)
                        gen_end_int = int(gen_end)
                        
                        #assert that exon genomic coordinates are within transcript genomic start and end coordinates
                        assert (gen_st_int > ref_start or gen_st_int < ref_end), "calculated genomic position of exon start is not within given genomic coordinates"
                        assert (gen_end_int > ref_start or gen_end_int < ref_end), "calculated genomic position of exon end is not within given genomic coordinates"

                        if gen_end > gen_st:
                            # extra caution: check that genomic coordinates were correctly offset before writing to bed file
                            bed_list = [chr, gen_st, gen_end]
                            bedfile.write("\t".join(bed_list))
                            bedfile.write("\n")
                            
                        elif gen_end <= gen_st:
                            print ("error: end coord is not greater than start")
                            
                exon_ranges[exon_number]=[start_int,end_int]

    return exon_ranges

def get_diffs(exon_ranges, gene, root):
    """
    produces a csv file of the differences in the gene between build 37 and build 38

    parameters: exon_ranges(dict), gene(string), root(elementTree root node)


    """
    # initialise data structures
    diffexons={} # dictionary linking differences to exon number; value defaults to 'intronic' if position do not map within exon coordinates
    lrgstartlist=[] # list of differences' start position
    
    # define file name for differences file and create list of column headers to write to file    
    diff_file = gene + "_diffs.csv"
    diff_headers = ["position", "type", "lrg_start", "lrg_end", "other_start", "other_end", "LRG_seq", "other_seq"]
    
    # write headers into csv file
    write_csv(diff_headers, diff_file, 'w') # mode 'w' truncates any file with same name in the directory

    for mapping in root.findall("./updatable_annotation/annotation_set[@type='lrg']/mapping[@type='main_assembly']"):
        # get reference assembly id and print to console
        ref_assem = (mapping.attrib['coord_system'])
        print ("Reference assembly= "+ref_assem +" so differences are with respect to this bulid")
        
        # for each difference in xml, get attributes into list
        for span in mapping:
            for diff in span:
                lrg_start = int(diff.attrib['lrg_start'])
                lrgstartlist.append(lrg_start) #for each difference, append to list of start positions
                
                typeattrib = diff.attrib['type']
                lrg_start_str = (diff.attrib['lrg_start'])
                lrg_end = (diff.attrib['lrg_end'])
                other_start = (diff.attrib['other_start'])
                other_end = (diff.attrib['other_end'])
                LRG_seq = diff.attrib['lrg_sequence']
                other_seq = diff.attrib['other_sequence']

                for pos in lrgstartlist:
                    # check if start position of difference is within exon start and end coordinates
                    # using exon = (exon_start, exon_end) dictionary; values are a list;
                    for key, value in exon_ranges.items(): # for key, value pair in exon_ranges dictionary                     
                        # for exon in exon_ranges dictionary
                        # compare pos to exon_start and exon_end coordinates in all exons                       
                        if pos >= value[0] and pos <= value[1]: # value [0] refers to exon_start, value[1] refers to exon_end
                        # if pos within exon range,
                        # add pos = exon to diffexons dictionary
                            diffexons[pos] = 'exon '+str(key)
                            
                        else:
                            next

                    # if pos not in exons, it is not in diffexons yet
                    # must be added to diffexons as 'intronic'
                    if pos not in diffexons.keys():
                        diffexons[pos] = 'intronic'

                    for position, location in diffexons.items():
                        # for position in diffexons, location is either the exon number or 'intronic'
                        if position == lrg_start:
                            #pos_str = v
                            # set up list with attributes for k difference
                            diff_list = [location, typeattrib, lrg_start_str, lrg_end, other_start, other_end, LRG_seq, other_seq]
                            # write diff_list content to line in csv file                          
                            write_csv(diff_list, diff_file, 'a') #mode 'a' to append to existing file diff_file


def get_annotations(gene, root):    
   """
    Outputs gene annotations, including overlapping genes and respective gene name synonyms, to CSV file
    
    parameters: gene(string), root(elementTree root node)
    
   """
   #initialise file with headers
   annot_file = gene+"_annotation.csv"
   annot_headers = ['NCBI_ID','HGNC_symbol', 'LRG_start','LRG_end','Strand','Description','Synonyms' ]
   write_csv (annot_headers, annot_file, 'w')#write mode ('w') truncates file with same name in directory to avoid appending to an old file

   # loop through LRG file to get annotations into annotation_list
   # one list for each overlapping gene (if present)
   for gene in root.findall("./updatable_annotation/annotation_set[@type='ncbi']/features/gene"):
        annotation_list=[]
        
        if gene.attrib.get('source')=='NCBI-Gene':
            NCBI = gene.attrib['accession']
            #print (NCBI)
            annotation_list.append(NCBI)
        else: #in case no NCBI accession
            annotation_list.append('')
            
        for branch in gene:
            if branch.tag == 'symbol':
                    
                if branch.attrib.get('source')=='HGNC':
                    HGNC = branch.attrib['name']
                    annotation_list.append(HGNC)
                    #print (HGNC)
                    
                    #create a synonyms list
                    synonym_list=[]
                    for synonym in branch:
                        synonym_list.append(synonym.text)
                        
                    synonym_string = ', '.join(synonym_list) #convert synonym_list into string to append later into annotation_list as one element
                    #print ('synonyms= ', synonym_string)
                   
                    
                    
                else:#in case no HGNC name
                    annotation_list.append('')
                    synonym_string = ''

                    
            if branch.tag == 'coordinates':
                LRG_start = branch.attrib['start']
                LRG_end = branch.attrib['end']
                Strand = branch.attrib['strand']
                                
                annotation_list.append(LRG_start)
                annotation_list.append(LRG_end)
                annotation_list.append(Strand)
                #print (LRG_start+'-'+LRG_end+'; strand='+Strand)
                                
            if branch.tag == 'long_name':
                ln = branch.text
                annotation_list.append(ln)
                                    
                annotation_list.append(synonym_string)
                #print (annotation_list)
                                    
                write_csv (annotation_list, annot_file, 'a') #mode 'a' to append to annot_file


def main():
	"""
	Parses and handles command line arguments. Calls the other functions in the script.
	Functions produce a bedfile, an optional annotation file and an optional information file
	"""
	
	#parses the command line arguments to check that all flags passed are valid, exits if not
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'hg:d:i:', ['help', 'gene=', 'difference=', 'info='])
	except getopt.GetoptError as err:
		print (err)
		usage()
		sys.exit(2)
	#defines all parameters to allow the possibility of optional arguments
	genein = ''
	diff = "false"
	info = "false"
	for opt, arg in opts:
		if opt == '-h':
			print (__doc__)
			sys.exit(2)
		elif opt in ('-g', '--gene'):
			genein = arg
		elif opt in ('-d', '--difference', '-diff'):
			diff = arg
		elif opt in ('-i', '--information', '--info'):
			info = arg
		else:
			usage()
			sys.exit(2)
	#checks if an LRG file name has been passed, exits if not
	if genein == '':
		print ('Please supply LRG file name without extension')
		usage()
		sys.exit(2)
		
	# read xml; function returns root object and variable with gene name
	root, gene = read_file(genein)
	# create bed file and return dictionary of exon ranges
	exon_ranges = bed_file(root, gene)

	# create csv file with sequence differences
	if diff == "True":
		get_diffs(exon_ranges, gene, root)
	else:
		print ("-d = ", diff, " therefore no annotation file produced")

	# create csv file with annotations for overlaping genes and respective synonyms	
	if info == "True":
		get_annotations(gene, root)
	else:
		print ("-i is not True therefore no annotation file produced")
	
def usage():
	"""
	helpful hints about the usage of the script
	"""

	print ("usage:")
	print ("python LRGparser.py -g [LRG file name] -d [True/False] -i [True/False]")
		
#runs script if it is run as a script from the command line as opposed to as a function
if __name__ == "__main__":
    main()





