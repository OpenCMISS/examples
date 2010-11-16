#! /usr/bin/env python
#
# Script to aid in setting up a 42 Testing points example.
#
# Bojan Blazevic
# Nov 2010
#
#
# USAGE:
# ./make42tree.py EXECUTABLE ARGS ROOT_DIR
#
#   EXECUTABLE :  Path to the compiled example binary
#                 This must be with respect to OPENCMISS_ROOT,
#                 which will be prepended to this string.
#   ARGS       :  Default arguments that the example requires [= '']
#   ROOT_DIR   :  Directory to build the folder structure [= ./]
#
# This file can also be loaded as a module, so that more parameters can be
# changed, individual functions run etc.
#
#
# BEHAVIOUR:
# 1) Creates the directory structure as defined by the 'spec' argument of the 
#    'traverseDirTree' function.
# 2) Creates the 'input' 'output' and 'expected_results' folders in each leaf
#    folder. Creates a 'run42.sh' that points to the executable passed in on
#    the command line.
# 3) Runs every 'run42.sh' in the directory tree
#
#
# USE CASES:
# 'I have an example file that takes care of all 42 points in one go. It is
#  located in $OPENCMISS_ROOT/cm/examples/.../jimmysexample/':
#  ./make42tree.py /cm/examples/.../jimmysexample/bin/x86_64-linux/mpich2/gnu/jimmysexample
#
# 'I have made some changes and I need to regenerate my expected output':
# python
# import make42tree
# run42FirstTime()
#
# 'But I don't want to re-run _everything_!':
# Read the 'TODO' below
#
#
# ARGUMENT FORMAT:
# -DIM = 2D | 3D
# -ELEM = TRI | TET | QUAD | HEX | HERMITE
# -BASIS = linear | quadratic | cubic | hermite
# -LEVEL = 1 | 2 | 3
#
#
# TODO:
# - Add docstrings
# - Add 'exceptions' (folders to be skipped when traversing the tree). We kinda
#   need this for steps 2/3 above.
# - Add parameter checking/automatic replacement of path with $OPENCMISS_ROOT


import sys
import os
import subprocess
import stat


# This function traverses a directory tree, calling 'nodeFunc' in every folder
# that has subfolders and calling 'leafFunc' in the remaining folders. By
# changing the functions passed in 'traverseDirTree' can be used for each step
# of setting up a '42 Testing Points' example and the structure specification
# format allows flexibility in the exact structure of the tree. (Eg. 'cubic'
# was changed to 'CubicVelocityLinearPressure' in the Navier Stokes tree)
#
# The tree structure specification is defined below
#
#     tree   ::= [level, level, ...]
#     level  ::= [option, option, ...]
#              | (tree, tree, ...)
#
# Here 'option' must be a string - this denotes an actual folder name.
#
# INTERPRETATION:
# This structure in essence just defines a tree where each folder has a sub-
# folders that are defined by the next 'level' spec. in the list. 
#
# Eg. [['left', 'right'], ['up', 'down'], ['forward', 'reverse']] would
#     define a tree where 'left' and 'right' had as subfolders 'up' and 'down'
#     etc.
#
# We can define more complex trees by throwing in tuples of trees where lists
# of folder names once were. A tuple of trees effectively results in each tree
# being substituted into the whole spec and executed as above.
#
# Eg. [['L', 'R'], ([['U'], ['F', 'R']], [['D'],['N', 'S']]), ['foo', 'bar']]
#     results in the following two trees:
#     [['L', 'R'], ['U'], ['F', 'R'], ['foo', 'bar']]
#     and
#     [['L', 'R'], ['D'], ['N', 'S'], ['foo', 'bar']]
#
# Notice how difficult this is to read and write. One should probably proceed
# by carefully modifying the specification defined later in the document.

def traverseDirTree(spec, funcs, sofar = []):
    if len(spec) == 0:              # Leaf node
        funcs.leafFunc(sofar)
    else:
        cur = spec[0]
        if type(cur) == tuple:      # Tuple trees
            for t in cur:
                traverseDirTree(t + spec[1:], funcs, sofar)
        else:                       # An actual tree
            for n in cur:
                funcs.nodeFunc(sofar + [n])
                try:
                    os.chdir(n)
                    traverseDirTree(spec[1:], funcs, sofar + [n])
                    os.chdir('..')
                except OSError:
                    print('%s does not exist, continuing...' % n)



class makeDirFuncs(object):
    def nodeFunc(self, sofar):
        d = sofar[-1]
        print('Making directory %s' % d)
        try:
            os.mkdir(d)
        except OSError:
            print('Cannot create %s (or it exists), continuing...' % d)
            
    def leafFunc(self, sofar): pass


    
class makeScriptFuncs(object):
    def __init__(self, defaultargs = '', execdir = '/'):
        self.defaultargs = defaultargs
        self.execdir = execdir
    
    @staticmethod
    def makeArgs(sofar):
        # The easiest way to do this
        args = ''
        args += '-DIM=' + sofar[0] + ' '
        args += '-ELEM=' + sofar[1] + ' '
        if (sofar[1] == 'HERMITE'):
            args += '-BASIS=hermite '
        else:
            args += '-BASIS=' + sofar[2] + ' '
        args += '-LEVEL=' + sofar[-1][-1]
        return args
    
    @staticmethod
    def makeLocalDirs(sofar):
        try:
            os.mkdir('input')
            os.mkdir('output')
            os.mkdir('expected_results')
        except OSError:
            print('Error making leaf folders, continuing...')
    
    def nodeFunc(self, sofar): pass
        
    def leafFunc(self, sofar):
        print('Making leaf %s' % '/'.join(sofar))
        makeScriptFuncs.makeLocalDirs(sofar)
        options = makeScriptFuncs.makeArgs(sofar)
    
        # Make the run script
        try:
            run42 = open('run42.sh', 'w')
        except IOError:
            print('Cannot make run42.sh in leaf, continuing...!?')
        else:
            run42.write('#!/bin/bash\n')
            # Should look at making following line more robust, people will most likely make mistakes
            run42.write('$OPENCMISS_ROOT/' + self.execdir + ' ' + self.defaultargs + ' ' + options + '\n')
            run42.write('mv *.exnode *.exelem output/\n')
            run42.close()
            
            # This method did not work, need to read more help
            # os.chmod('run42.sh', stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            os.chmod('run42.sh', 0755)
            # Neither does this
            # subprocess.call('chmod +x run42.sh')
        
        
        
class firstTimeRunFuncs(object):
    def nodeFunc(self, sofar): pass
    
    def leafFunc(self, sofar):
        print('Running %s' % '/'.join(sofar) + '/run42.sh')
        try:
            ret = subprocess.call('/bin/bash run42.sh', shell=True)
        except Exception:
            print('Could not run subprocess')
        else:
            if not ret == 0:
                print('Bad return code %d, continuing...' % ret)
            else:
                subprocess.call('mv *.exnode *.exelem expected_results')
     
        

# The specification for the default 42 Testing Points tree
spec42 = [(
           [['2D'], (
                     [['TRI', 'QUAD'], ['linear', 'quadratic', 'cubic']],
                     [['HERMITE']]
                    )
           ],
           [['3D'], (
                     [['TET', 'HEX'], ['linear', 'quadratic', 'cubic']],
                     [['HERMITE']]
                    )
           ]
          ), 
           ['LEVEL_1', 'LEVEL_2', 'LEVEL_3']
         ]
 


# Nice wrappers for the 3 main stages of setup
def make42Structure(spec = spec42):
    traverseDirTree(spec, makeDirFuncs())
    
def make42Scripts(defaultargs, execdir, spec = spec42):
    traverseDirTree(spec, makeScriptFuncs(defaultargs, execdir))
    
def run42FirstTime(spec = spec42):
    traverseDirTree(spec, firstTimeRunFuncs())


# When used as a script
if __name__ == '__main__':

    def printUsageDie():
        print('make42tree.py')
        print('')
        print('Usage:')
        print('  ./make42tree.py EXECUTABLE ARGS ROOT_DIR')
        print('')
        print('   EXECUTABLE :  Path to the compiled example binary')
        print('                 This must be with respect to OPENCMISS_ROOT,')
        print('                 which will be prepended to this string.') 
        print('   ARGS       :  Default arguments that the example requires '
                                  '[= \'\']')
        print('   ROOT_DIR   :  Directory to build the folder structure '
                                  '[= ./]')
        sys.exit('Incorrect usage!')
    
    if len(sys.argv) < 2:
        printUsageDie()
        
    execdir = sys.argv[1]
    
    if len(sys.argv) < 3:
        rootdir = './'
    else:
        rootdir = sys.argv[2]
        
    if len(sys.argv) < 4:
        defaultargs = ''
    else:
        defaultargs = sys.argv[3]
    
    try:
        os.chdir(rootdir)
    except OSError:
        print('Cannot go to ROOT_DIR=%s' % rootdir)
        printUsageDie()
    
    make42Structure()
    make42Scripts(defaultargs, execdir)
    run42FirstTime()


