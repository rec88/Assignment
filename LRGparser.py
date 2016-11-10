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
                    if mapping.attrib['coord_system'] == "GRCh37.p13":
                        chr = mapping.attrib['other_name']

    for id in root.findall("./fixed_annotation/id"):
        id_tag = id.text

    for exon in root.findall("./fixed_annotation/transcript/exon"):
        exon_number = exon.attrib['label']
        #print(exon_number)
        for coord_sys in exon.iterfind('coordinates[@coord_system="LRC_7"]'):
            start = coord_sys.attrib['start']
            end = coord_sys.attrib['end']
            if end > start:
                #coords = [start, end]
                bed_list = [chr, start, end]
                build37bed.write("\t".join(bed_list))
                build37bed.write("\n")
            elif end <= start:
                print ("error: end coord is not greater than start")
        exon_ranges[range(int(start), int(end))] = exon_number

    #for keys, values in exon_ranges.items():
    #    print(keys)
    #    print(values)
    return exon_ranges


root, gene = read_file()
exon_ranges = bed_file(root, gene)
