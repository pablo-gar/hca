# Gets a scRNA-seq matrix from DCP using matrix service by project title

import hca
import re
import sys
import numpy as np
import loompy as lp
from shutil import copy2
#from dcplib.etl import DSSExtractor

# Globs
SWAGGER_URL = "https://dss.data.humancellatlas.org/v1/swagger.json"
REPLICA = "aws"
LOOM_UUID = "dss_bundle_fqid"

def main():
    
    loom_in_file = sys.argv[1]
    loom_out_file = sys.argv[2]
    
    copy2(loom_in_file, loom_out_file)
    
    loom_out = lp.connect(loom_out_file)
    
    # get sra ids
    new_features = get_sra(loom_out)
    
    # adds them to file
    loom_out.ca.insdc_run_accessions = new_features["insdc_run_accessions"]
    loom_out.ca.file_name = new_features["file_name"]
    
    loom_out.close()
    

def get_sra(loom):
    
    client = hca.dss.DSSClient(swagger_url = SWAGGER_URL)

    # Fields of interest from the sequence_file_0.json
    fields =  ["insdc_run_accessions", "file_name"]
    
    # Initialize new cellular features 
    features = dict()
    for i in fields:
        features[i] = [""] * loom.shape[1]
    
    for i in range(loom.shape[1]):
        
        print("Working with cell", i, end = "\r")
        
        # Get bunddle id
        uuid = loom.ca[LOOM_UUID][i].split(".")[0]
        
        # Get json uuid
        bundle = client.get_bundle(uuid = uuid, replica = REPLICA)
        sequence_file_uuid = get_file_uuid(bundle, "sequence_file_0.json")
        
        # Get the file
        seq_file = client.get_file(uuid = sequence_file_uuid, replica = REPLICA)
        
        
        # Getting sra ids and seq names
        features["file_name"][i] = seq_file["file_core"]["file_name"]
        if "insdc_run_accessions" in seq_file:
            features["insdc_run_accessions"][i] = seq_file["insdc_run_accessions"][0]
        else:
            features["insdc_run_accessions"][i] = re.sub(r'(SRR[a-zA-Z0-9]+).*', r'\1', seq_file["file_core"]["file_name"])
            
    
    print("")
    return features
        

def get_file_uuid(bundle, name):
    
    ''' 
    given a dict bundle from hca.dss.get_bundle() gets the uuid for the sequence_file_0 file
    '''
    
    uuid = ""
    for file_json in bundle['bundle']['files']:
        if file_json["name"] == name:
            uuid = file_json["uuid"]
            break
    
    if uuid == "":
        raise Exception("No file found for " + name)
    
    return uuid

def test():
    
    '''
    testing an elastic query search for faster results
    '''
    
    extractor = DSSExtractor(staging_directory=staging_directory,
                             content_type_patterns=content_type_patterns,
                             filename_patterns=filename_patterns,
                             dss_client=get_dss_client(deployment_stage))
    extractor.extract(
        query=query, #elastic search query
        transformer=transformer_cb,
        finalizer=finalizer_cb,
        max_workers=max_workers,
        max_dispatchers=max_dispatchers,
        dispatch_executor_class=dispatcher_executor_class
    )
        
if __name__ == "__main__":
    main()
