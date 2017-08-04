#!/usr/bin/python
#
# Author: Mike Gloudemans
# 
# Tabix lists of SNPs for easy access.
# 

# Function tabix_file
# 
# Arguments: Config object
# Action: tabix all files in Config.
# Returns: None
def tabix_all(config):

    # Call tabix_file for each file in config.
    for file in eqtl_files:
        tabix_file(file)


# Arguments: name of file, whether to reindex if file exists already
# Action: tabix the file
# Returns: None
def tabix_file(filename, redo=FALSE):

    # Check whether file has already been tabixed. If so,
    # quit.
    if file.exists("{0}.tbi".format(filename))
        return

    # Otherwise, tabix the file.
    print("Indexing files with tabix. This may take several hours, but only is necessary \
            on the first run with a given file.")

    # Run tabix
    # TODO: Look back at other code to see the right format.
    # Note that tabix will require inputs formatted in a specific order. 
    
