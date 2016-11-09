import elementtree.ElementTree as ET
import sys, os

tree = ET.parse("../LRG_7.xml")
doc = tree.getroot()
thingy = doc.find('timeSeries')

print thingy.attrib