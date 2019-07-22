# Takes a reference and a target loom file, and based on a shared col attribute
# appends a column attribute from the reference to the target file
# 
# Positional ars
# key id for refence loom file 
# value id for reference loom file -- this columns will be join to target value
# key id for target loom file, should have shared values with key id for reference
# path to loom reference file
# path to loom target file
# path to loom output file


import sys
import loompy as lp
from shutil import copy2

def main():
    
    # Read args
    key_ref = sys.argv[1]
    value_ref = sys.argv[2]
    key_target = sys.argv[3]
    loom_ref_file = sys.argv[4]
    loom_target_file = sys.argv[5]
    loom_out_file = sys.argv[6]
    
    # Getting files ready 
    copy2(loom_target_file, loom_out_file)
    loom_ref = lp.connect(loom_ref_file)
    loom_out = lp.connect(loom_out_file)
    
    
    # Reading key-value pairs from the refrence file
    attr = dict()
    for key,value in zip(loom_ref.ca[key_ref], loom_ref.ca[value_ref]):
        attr[key] = value
    
    # Order and write values in target file
    new_col = [ attr[key] if key in attr else None for key in loom_out.ca[key_target] ]
    loom_out.ca[value_ref] = new_col

if __name__ == "__main__":
    main()
    
