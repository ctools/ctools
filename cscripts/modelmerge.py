#!/usr/bin/env python
import gammalib
import sys

# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Initialise flags
    need_help = False

    # Test for command line arguments
    if (len(sys.argv) > 1):
        if sys.argv[1] == "-h":
            need_help = True
    else:
        need_help = True

    # Print help if needed and exit
    if need_help:
        print("This script merges all model XML files provided as command line arguments. The usage is:\n")
        print("$ modelmerge.py modelfile1.xml modelfile2.xml ...\n")
        print("The resulting model will be written into \"merged_models.xml\" by default")
        print("This can be modified by specifying the keyword \"outmodel\", e.g.\n")
        print("$ modelmerge.py modelfile1.xml modelfile2.xml outmodel=my_merged_models.xml")
        sys.exit()

    # Get arguments
    args = sys.argv[1:]
    
    # Initialise empty model container
    models = gammalib.GModels()
    
    # Initialise outmodel file
    outmodel = "merged_models.xml"
    
    # Loop over input parameters
    for arg in args:
        if "outmodel=" in arg:
            outmodel = arg.replace("outmodel=","")
            continue    
        
        # Load model file
        add_model = gammalib.GModels(arg)
        
        # Append each model to output container
        for model in add_model:
            models.append(model)
    
    # Save models    
    models.save(outmodel)