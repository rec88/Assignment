import xml.etree.ElementTree as ET
import sys, os

tree = ET.parse("/media/sf_Programming/Desktop/LRG_292_BRCA1.xml")
source = tree.find('sequence_source').text
print 'source: ', source

