"""
LRGparser_testsuite.py

created: December 2016

This unit test script assumes you are using python version 3.5.2 
It tests the function of LRGparser.py
LRGparser.py is tested on python version 3.5.2 and python 2.7.8


If you are using python 2.7.8 please change lines 36 and 93, and the boolean at line 56 to 1 to reflect this.
 
@authors: Laura Carreto, Rosie Coates-Brown


usage: python LRGtesting.py
"""

import unittest, sys, os
from LRGparser import read_file
from LRGparser import versiontest
from LRGparser import bed_file
from LRGparser import get_diffs
from LRGparser import get_annotations
from LRGparser import check_status

def version():
    version =".".join(map(str, sys.version_info[:3]))
    if version in {"3.5.2"}:
        print ("goodnews! you are using a compatible version of python: ", version,)
    else:
        print (__doc__)
        sys.exit(2)
 

class testLRG(unittest.TestCase):
 
    def setUp(self):
        pass

    #when unittest is imported all functions with test as the first 4 letters are run as a test
    def testread_fileLRG7(self):
        """tests the read file function by checking the file name generated. 
        this function is required for other tests, which I don't think is ideal.""" 
        #calls the read_rile function from LRGparser.py
        actual1, actual2, actual3 = read_file("LRG_7", "url", "3.5.2")
        expected = "LRG_7.xml"
        #asserts if the actual and expected values are not equal
        self.assertEqual(actual3, expected)
        
        
    def testbedfunction(self):
        """tests the bed_file function is working correctly by checking the strand. 
        This is vital for the correct offset calculation when generating genomic coordinates.
        This function is required for other tests which is not ideal"""
        testroot, testgene, testname = read_file("LRG_7", "url", "3.5.2")
        #calls the bedfile function in LRGparser.py
        actual1, actual2 = bed_file(testroot, testgene)
        expected_str = "-1"
        self.assertEqual(actual2, expected_str)
        
        
    def testversion(self):
        """tests the testversion function is correctly generating the python version number"""
        actual_v, version = versiontest()
        expected_v = 0
        self.assertEqual(actual_v, expected_v)
        
    def testdiffs(self):
        """tests the get_diffs function by checking the expected diff list is produced for LRG_7"""
        testroot, testgene, testname = read_file("LRG_7", "url", "3.5.2")
        test_exonranges, actual2 = bed_file(testroot, testgene)

        expected_list = ['intronic', 'mismatch', '145881', '145881', '13365580', '13365580', 'A', 'G']
        actual_list = get_diffs(test_exonranges, testgene, testroot) 

        self.assertEqual(actual_list, expected_list)
        

    def testannotation(self):
        """tests the get_annotations function by checking the expected last long name is produced for LRG_7"""
        testroot, testgene, testname = read_file("LRG_7", "url", "3.5.2")
        

        expected_name = "calcium voltage-gated channel subunit alpha1 A"
        actual_name =  get_annotations(testgene, testroot) 

        self.assertEqual(actual_name, expected_name)    
        
    def teststatus(self):
        """tests the check_status function by checking the expected status of LRG_7"""
        
        actual_status = check_status("url", "3.5.2", "LRG_7")
        expected_status = "public"
        
        self.assertEqual(actual_status, expected_status)

if __name__ == '__main__':
    version()
    unittest.main()