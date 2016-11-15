"""
created: November 2016

@authors: Laura Carreto, Rosie Coates-Brown

usage: python LRGparser.py --gene --difference --info

--gene = name of LRG file without .xml suffix
--difference = (y/n) y will trigger the output of a csv file containing the differences between 37 and 38, n will suppress this
--annotations = (y/n) y will cause a text file of gene information to be produced, n will suppress this file. This file contains information on synonyms, lsdb, long gene name

"""

try:
    import xml.etree.ElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import sys, os, csv


def read_file():
    """
    read in the LRG.xml file
    test: does the file exist

    input: sys.argv[1]
    returns: root(ElementTree root node), gene(user input string variable)

    """
    gene = 'LRG_7'#sys.argv[1]
    file_name = gene+'.xml'
    file_path = '/home/swc/Desktop/LRGParser/'#'/Users/rosiecoates/Documents/Clinical_bioinformatics_MSc/programming/assignment/'
    full_path = file_path+file_name
    try:
        tree = ET.parse(full_path)
    except:
        print("couldn't open file...")
    tree = ET.ElementTree(file=full_path)
    root = tree.getroot()

    return root, gene
    
def write_csv(mylist, myfilename, mode):
    """
    Output to CSV file from a list
    
    Parameters: mylist (list), myfilename(string), mode (string); the mode options used in the code are: 'a'= append, 'w'=write
    
    """
    
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
    exon_ranges={}


    for mapping in root.findall("./updatable_annotation/annotation_set[@type='lrg']/mapping[@type='main_assembly']"):

        ref_start = int(mapping.attrib['other_start'])
        #print (ref_start)
        ref_end = int(mapping.attrib['other_end'])
        #print(ref_end)
        offset = (ref_start) - 1
        offset_int = int(offset)
        #print (offset)
        chr = mapping.attrib['other_name']

    for id in root.findall("./fixed_annotation/id"):
        id_tag = id.text

    for transcript in root.findall("./fixed_annotation/transcript"):
        trans_num =  (transcript.attrib['name'])
        ref_name = gene + "_" + trans_num + ".bed"
        build37bed = open(ref_name, 'w')
        for exon in transcript:
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

                        assert (gen_st_int > ref_start or gen_st_int < ref_end), "claculated genomic position of exon start is not within given  genomic coordinates"
                        assert (gen_end_int > ref_start or gen_end_int < ref_end), "claculated genomic position of exon end is not within given  genomic coordinates"

                        if gen_end > gen_st:

                            bed_list = [chr, gen_st, gen_end]
                            build37bed.write("\t".join(bed_list))
                            build37bed.write("\n")
                        elif gen_end <= gen_st:
                            print ("error: end coord is not greater than start")
                exon_ranges[exon_number]=[start_int,end_int]


    return exon_ranges

def get_diffs(exon_ranges):
    """
    produces a csv file of the differences in the gene between build 37 and build 38

    parameters: exon_ranges(dict)


    """
    diffexons={}
    lrgstartlist=[]
    
    diff_file = gene + "_diffs.csv"
    diff_headers = ["position", "type", "lrg_start", "lrg_end", "other_start", "other_end", "LRG_seq", "other_seq"]
    write_csv(diff_headers, diff_file, 'w') #mode 'w' to write to diff_file; truncates any file with same name in the directory

    for mapping in root.findall("./updatable_annotation/annotation_set[@type='lrg']/mapping[@type='main_assembly']"):

        ref_assem = ("reference assembly="+(mapping.attrib['coord_system']))
        print (ref_assem)
        for span in mapping:
            for diff in span:
                lrg_start = int(diff.attrib['lrg_start'])
                lrgstartlist.append(lrg_start)
                typeattrib = diff.attrib['type']
                lrg_start_str = (diff.attrib['lrg_start'])
                lrg_end = (diff.attrib['lrg_end'])
                other_start = (diff.attrib['other_start'])
                other_end = (diff.attrib['other_end'])
                LRG_seq = diff.attrib['lrg_sequence']
                other_seq = diff.attrib['other_sequence']

                for pos in lrgstartlist:
                    for key, value in exon_ranges.items():
                        if pos >= value[0] and pos <= value[1]:
                            diffexons[pos] = [key]
                        else:
                            diffexons[pos] = ['intronic']

                    for k, v in diffexons.items():
                        if k == lrg_start:
                            pos_str = v[0]
                            diff_list = [pos_str, typeattrib, lrg_start_str, lrg_end, other_start, other_end, LRG_seq, other_seq]
                            write_csv(diff_list, diff_file, 'a') #mode 'a' to append to existing file diff_file


def get_annotations(gene):    
   """
    Outputs gene annotations, including overlapping genes and respective gene name synonyms, to CSV file
    
   """
   #initialise file with headers
   annot_file = gene+"_annotation.csv"
   annot_headers = ['NCBI_ID','HGNC', 'LRG_start','LRG_end','Strand','Description','Synonyms' ]
   write_csv (annot_headers, annot_file, 'w')#write mode ('w') truncates file with same name in directory to avoid appending to an old file

#loop through LRG file to get annotations into annotation_list; one list for each overlapping gene (if existing)
   for gene in root.findall("./updatable_annotation/annotation_set[@type='ncbi']/features/gene"):
        annotation_list=[]
        
        if gene.attrib.get('source')=='NCBI-Gene':
            NCBI = gene.attrib['accession']
            #print (NCBI)
            annotation_list.append(NCBI)
        else: #in case no NCBI accession
            annotation_list.append('')
            
        for leaf in gene:
            if leaf.tag == 'symbol':
                    
                if leaf.attrib.get('source')=='HGNC':
                    HGNC = leaf.attrib['name']
                    annotation_list.append(HGNC)
                        #print (HGNC)
                else:#in case no HGNC name
                    annotation_list.append('')
                        
                    #create a synonym list
                synonym_list=[]
                for synonym in leaf:
                    synonym_list.append(synonym.text)
                    synonym_string = ', '.join(synonym_list) #convert synonym_list into string to append later into annotation_list as one element
                #print ('synonyms= ', synonym_string)
                    
            if leaf.tag == 'coordinates':
                LRG_start = leaf.attrib['start']
                LRG_end = leaf.attrib['end']
                Strand = leaf.attrib['strand']
                                
                annotation_list.append(LRG_start)
                annotation_list.append(LRG_end)
                annotation_list.append(Strand)
                #print (LRG_start+'-'+LRG_end+'; strand='+Strand)
                                
            if leaf.tag == 'long_name':
                ln = leaf.text
                annotation_list.append(ln)
                                    
                annotation_list.append(synonym_string)
                #print (annotation_list)
                                    
                write_csv (annotation_list, annot_file, 'a') #mode 'a' to append to annot_file



root, gene = read_file()
exon_ranges = bed_file(root, gene)
get_diffs(exon_ranges)
get_annotations(gene)