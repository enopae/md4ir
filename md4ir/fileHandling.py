"""
Generally used file handling functions.
This code is part of md4ir, available at:
https://gitlab.ethz.ch/paenurke/md4ir
"""

__all__ = ['read_file','read_cols']


# Functions

def read_file(datafile):
    """
    Reads the content of a file and splits the line content into a list
    """
    content = []
    with open(datafile, encoding='utf8') as f:
        for line in f:
                content.append(line.split())
        f.close()
    return content

def read_cols(datafile):
    """
    Returns the first 3 columns of the data file if more than 1 columns
    Else returns the first column
    No other number of columns are currently supported
    """
    data = read_file(datafile)
    if len(data[0]) > 1:
        data = [[i[0],i[1],i[2]] for i in data]
    else:
        data=[[i[0]] for i in data]
    return data