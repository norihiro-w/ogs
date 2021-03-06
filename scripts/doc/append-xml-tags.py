#!/usr/bin/python

# This script augments the parameter documentation pages by information
# such as if they are required/optional, their data typ and in which
# end-to-end tests they are used.
# It uses the cache file generated by normalize-param-cache.py

# prevent broken pipe error
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

import os
import sys
import xml.etree.cElementTree as ET

github_src_url  = "https://github.com/ufz/ogs/tree/master"
github_data_url = "https://github.com/ufz/ogs-data/tree/master"

if len(sys.argv) != 4:
    print("Usage:")
    print("{0} EXT DATADIR DOCAUXDIR".format(sys.argv[0]))
    sys.exit(1)

ext       = sys.argv[1]
datadir   = sys.argv[2]
docauxdir = sys.argv[3]

extension = '.' + ext
datadir   = os.path.abspath(datadir)
docauxdir = os.path.abspath(docauxdir)
docdir    = os.path.join(docauxdir, "dox", "ProjectFile")

# used to expand documentation entry points to full xml tag paths
# that are used in the prj file.
tag_path_expansion_table = {
    "initial_condition":  "process_variables.process_variable.initial_condition",
    "boundary_condition": "process_variables.process_variable.boundary_conditions.boundary_condition",
    "linear_solver":      "linear_solvers.linear_solver",
    "process":            "processes.process",
    "parameter":          "parameters.parameter",
    "prj": "",
}

# maps tags to the set of xml files they appear in
dict_tag_files = dict()

# maps tags to additional parameter info obtained prior to this script
dict_tag_info = dict()

def dict_of_set_append(dict_, key, value):
    if key in dict_:
        dict_[key].add(value)
    else:
        dict_[key] = set((value,))

def dict_of_list_append(dict_, key, value):
    if key in dict_:
        dict_[key].append(value)
    else:
        dict_[key] = [value]


def print_tags(node, path, level, filepath):
    global dict_tag_files

    tag = node.tag
    if level>1: # skip root node
        tagpath = path + "." + tag
    else:
        tagpath = tag

    if level>0: # skip root node
        dict_of_set_append(dict_tag_files, (True, tagpath), filepath)
        for k in node.attrib:
            dict_of_set_append(dict_tag_files, (False, tagpath + "." + k), filepath)

    for child in node:
        print_tags(child, tagpath, level + 1, filepath)

# gather info from xml files
for (dirpath, _, filenames) in os.walk(datadir):
    for f in filenames:
        if not f.endswith(extension): continue

        filepath = os.path.join(dirpath, f)
        xmlroot = ET.parse(filepath).getroot()
        print_tags(xmlroot, "", 0, filepath[len(datadir)+1:])

if False:
    first = True
    for (tag, files) in sorted(dict_tag_files.items()):
        if first:
            first = False
        else:
            print()

        print("T |" if tag[0] else "A |", tag[1])
        for f in sorted(files):
            print("   ", f)

# read parameter cache (generated by normalize-param-cache.py)
with open(os.path.join(docauxdir, "documented-parameters-cache.txt")) as fh:
    for line in fh:
        line = line.strip().split("@@@")
        if line[0] == "OK":
            tagpath = line[3]
            dict_of_list_append(dict_tag_info, tagpath, line)

# traverse dox file hierarchy
for (dirpath, _, filenames) in os.walk(docdir):
    reldirpath = dirpath[len(docdir)+1:]
    istag = True

    for f in filenames:
        if not f.endswith(".dox"): continue

        if f.startswith("i_") or f.startswith("c_"):
            tagpath = reldirpath
        elif f.startswith("t_"):
            tagpath = os.path.join(reldirpath, f[2:-len(".dox")])
            istag = True
        elif f.startswith("a_"):
            tagpath = os.path.join(reldirpath, f[2:-len(".dox")])
            istag = False

        tagpath = tagpath.replace(os.sep, ".")

        path = os.path.join(dirpath, f)
        with open(path, "a") as fh:
            # TODO this can currently only expand the top level
            tagpathparts = tagpath.split(".")
            if tagpathparts[0] in tag_path_expansion_table:
                tagpathparts[0] = tag_path_expansion_table[tagpathparts[0]]
            else:
                tagpathparts[0] = "NONEXISTENT"
            tagpath_expanded = ".".join(tagpathparts).lstrip(".")

            if tagpath:
                fh.write("\n\n# Additional info\n")
                if tagpath in dict_tag_info:
                    for info in dict_tag_info[tagpath]:
                        path = info[1]; line = info[2]
                        fh.write(("\n## From {0} line {1}\n\n")
                                .format(path, line))

                        method = info[6]
                        if method.endswith("Optional"):
                            fh.write("- This is an optional parameter.\n")
                        elif method.endswith("List"):
                            fh.write("- This parameter can be given arbitrarily many times.\n")
                        elif method: # method not empty
                            fh.write("- This is a required parameter.\n")

                        datatype = info[5]
                        if datatype: fh.write("- Data type: <tt>{0}</tt>\n".format(datatype))

                        fh.write("- Expanded tag path: {0}\n".format(tagpath_expanded))

                        fh.write("- Go to source code: [&rarr; ufz/ogs/master]({2}/{0}#L{1})\n"
                                .format(path, line, github_src_url))
                else:
                    fh.write("\nNo additional info.\n")

            if tagpath_expanded:
                fh.write("\n\n# Used in the following test data files\n\n")
                try:
                    datafiles = dict_tag_files[(istag, tagpath_expanded)]

                    for df in sorted(datafiles):
                        fh.write("- \\[[&rarr; ogs-data/master]({1}/{0})\\]&emsp;{0}\n"
                                .format(df, github_data_url))
                except KeyError:
                    fh.write("Used in no end-to-end test cases.\n")
            else:
                # no additional output for the main doc page
                pass

            fh.write("\n*/\n")
