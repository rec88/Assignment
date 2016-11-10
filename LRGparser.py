try:
    import xml.etree.ElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import sys, os


def read_file():
    """
    read in the LRG.xml file
    test: does the file exist ***dev note how do we handle the exceptions error correctly?***
    test: does fixed annotation ID match the user selected lrg file?
    """
    gene = sys.argv[1]
    file_name = gene+'.xml'
    file_path = '/Users/rosiecoates/Documents/Clinical_bioinformatics_MSc/programming/assignment/'
    full_path = file_path+file_name
    try:
        tree = ET.parse(full_path)
    except:
        print("couldn't open file...")
    tree = ET.ElementTree(file=full_path)
    root = tree.getroot()

    return root, gene


def bed_file(root, gene):
    """
    grabs the chromosome number, the exon number, the exon start and stop coordinates
    creates a dict with the exon ranges as keys, the exon_number as the value
    """
    exon_ranges={}
    ref_name = gene+"_37.bed"
    build37bed = open(ref_name, 'w')

    for annot_set in root.findall("./updatable_annotation/annotation_set"):
        if annot_set.attrib.get('type')=='lrg':
            for mapping in annot_set:
                if mapping.tag == 'mapping':
                    #if mapping.attrib['coord_system'] == "GRCh37.p13":
                    if mapping.attrib['type'] == "main_assembly":
                        ref_start = mapping.attrib['other_start']
                        offset = int(ref_start) - 1
                        offset_int = int(offset)
                        #print (offset)
                        chr = mapping.attrib['other_name']

    for id in root.findall("./fixed_annotation/id"):
        id_tag = id.text

    for exon in root.findall("./fixed_annotation/transcript/exon"):
        exon_number = exon.attrib['label']
        #print(exon_number)
        for coord_sys in exon:
            if (coord_sys.attrib['coord_system']) == "LRG_7":
                start=(coord_sys.attrib['start'])
                start_int = int(start)
                end = coord_sys.attrib['end']
                end_int = int(end)
                #exon_range = range(start_int,end_int)
                #print (exon_range)
                gen_st = str(start_int + offset_int)
                gen_end = str(end_int + offset_int)
                if gen_end > gen_st:
                    #coords = [start, end]
                    bed_list = [chr, gen_st, gen_end]
                    build37bed.write("\t".join(bed_list))
                    build37bed.write("\n")
                elif gen_end <= gen_st:
                    print ("error: end coord is not greater than start")
        exon_ranges[exon_number]=[start_int,end_int]

    #for keys, values in exon_ranges.items():
       #print(keys)
       #print(values)
    return exon_ranges

def get_diffs(exon_ranges):
    ref_name = gene + "_diffs.csv"
    diff_file = open(ref_name, 'w')
    diff_headers = ["type", "lrg_start", "lrg_end", "other_start", "other_end", "LRG_seq", "other_seq"]
    diff_file.write(",".join(diff_headers))
    diff_file.write("\n")

    for annot_set in root.findall("./updatable_annotation/annotation_set"):
        if annot_set.attrib.get('type')=='lrg':
            for mapping in annot_set:
                if mapping.tag == 'mapping':
                    if mapping.attrib['type'] == "main_assembly":
                       for span in mapping:
                           for diff in span:
                               type = diff.attrib['type']
                               lrg_start = int(diff.attrib['lrg_start'])
                               lrg_end =  int(diff.attrib['lrg_end'])
                               other_start = int(diff.attrib['other_start'])
                               other_end = int(diff.attrib['other_end'])
                               LRG_seq = diff.attrib['lrg_sequence']
                               other_seq = diff.attrib['other_sequence']
                               for key, value in exon_ranges.items():
                                   print(key)
                                   print("\t")
                                   print(value[0])
                                   print("\t")
                                   print(value[1])
                                   #if lrg_start >= value[0] and lrg_start <= value[1]:
                                        #print (key)
                                #csv file doesn't work as all variables are integers- remove int and it will be fine
                               #diff_list = [type, lrg_start, lrg_end, other_start, other_end, LRG_seq, other_seq]
                               #diff_file.write(",".join(diff_list))
                               #diff_file.write("\n")



root, gene = read_file()
exon_ranges = bed_file(root, gene)
get_diffs(exon_ranges)