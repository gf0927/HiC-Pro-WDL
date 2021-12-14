#!/usr/bin/env python

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant, Eric Viara
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

"""
Script to merge any files with the same template
"""

import getopt
import sys
import glob
import os
from collections import OrderedDict

def usage():
    """Usage function"""
    print("Usage : merge_statfiles.py")
    print("-d/--dir <files directory>")
    print("-p/--pattern <files pattern>")
    print("[-v/--verbose] <Verbose>")
    print("[-h/--help] <Help>")
    return


def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

def get_args():
    """Get argument"""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "d:p:vh",
            ["dir=",
             "pattern=",
             "verbose", "help"])
    except getopt.GetoptError:
        usage()
        sys.exit(-1)
    return opts

if __name__ == "__main__":
    ## Read command line arguments
    opts = get_args()
    verbose = False
    path = None
    pattern = None
    output = "-"

    if len(opts) == 0:
        usage()
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-d", "--dir"):
            path = arg
        elif opt in ("-p", "--pattern"):
            pattern = arg
        elif opt in ("-v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    ## Verbose mode
    if verbose:
        print("## merge_statfiles.py")
        print("## dir=", path)
        print("## pattern=", pattern)
 
    infiles = [name for name in glob.glob(os.path.join(os.path.abspath(path), pattern)) if os.path.isfile(os.path.join(path,name))]
    li = len(infiles)

    if li > 0:
        if verbose:
            print("## Merging "+ str(li)+" files")
 
        ## Reading first file to get the template
        template = OrderedDict()
        if verbose:
            print("## Use "+infiles[0]+" as template")
        with open(infiles[0]) as f:
            for line in f:
                if not line.startswith("#"):
                    lsp = line.strip().split("\t")
                    data = map(num, lsp[1:len(lsp)])
                    template[str(lsp[0])] = list(data)
        #取得total_reads
        for k,v in template.items():
            if k=="Total_pairs_processed":
                TotalReads = v[0]
        #print(template)
        if len(template) == 0:
            print("Cannot find template files !")
            sys.exit(1)

        ## Int are counts / Float are percentage
        for fidx in range(0, li):
            #把 template 打開
            with open(infiles[fidx]) as f:
                #針對每一行的template做處理
                for line in f:

                    #只處理非#開頭的行數
                    if not line.startswith("#"):
                        #把該行以'/t'做分隔，並儲存成list
                        lsp = line.strip().split("\t")
                        if lsp[0] in template:
                            for i in range(1, len(lsp)):
                                #針對整數的部分，通常是reads數量
                                if isinstance(num(lsp[i]), int):
                                    #直接把template給複製保存(不動)
                                    template[lsp[0]][i-1] = num(lsp[i])
                                else:
                                    template[lsp[0]][i-1] = round((num(lsp[i-1])/TotalReads)*100 ,3)
                        else:
                            sys.stderr.write("Warning : '"+lsp[0]+"' not found in template ["+infiles[fidx]+"]\n")


        
        ## Print template
        for x in template:
            sys.stdout.write(x)
            for y in template[x]:
                sys.stdout.write("\t"+str(y))
            sys.stdout.write("\n")

    else:
        print("No files to merge - stop")
        sys.exit(1)

