# Gets a scRNA-seq matrix from DCP using matrix service by project title
# Relies on REST API https://matrix.data.humancellatlas.org/v1, makes post request on project.project_core.project_title filter
#
# Usage
# python3 project_title outfile.loom
#
# Example
# python3 get_HCA_matrix_by_project_title.py "Ischaemic sensitivity of human tissue by single cell RNA seq" Ischaemic_sensitivity.loom

import requests, sys, time

# Globs
MATRIX_URL = "https://matrix.data.humancellatlas.org/v1"
FIELD = "project.project_core.project_title"

def main():
    
    project_title = sys.argv[1]
    out_file = sys.argv[2]
    
    print("Getting matrix:", project_title, sep = "\n")
    get_matrix(project_title, out_file)

def get_matrix(project_title, out_file):
    
    '''
    returns the filename of the downloaded matrix, throws error if it fails
    '''
    
    
    # POST
    resp = requests.post(
            MATRIX_URL + "/matrix",
            json={"filter":
                    {"op" : "=", 
                     "value" : project_title, 
                     "field" : FIELD
                    }
                 }
            )
                    
    
    # Wait for response
    resp = wait_response_matrix(resp)
    
    # Save matrix
    save_matrix(resp, out_file)
    

def wait_response_matrix(resp):
    
    '''
    waits for POST request to finish 
    '''
    
    while True:
        resp = update_response(resp)
        message = resp.json()["message"]
        status = resp.json()["status"]
        
        if status == "Complete":
            break
        if status == "Failed":
            raise Exception("Failed to get matrix with error: {}".format(message))
        
        print("Status:", status, "Waiting...")
        
        time.sleep(10)
    
    return resp

def save_matrix(resp, out_file):
    resp = update_response(resp)
    url = resp.json()["matrix_url"]
    
    with open(out_file, "wb") as out:
        out.write(requests.get(url).content)

    
def update_response(resp):
    return requests.get(MATRIX_URL + "/matrix/" + resp.json()["request_id"])
        
        
if __name__ == "__main__":
    main()
