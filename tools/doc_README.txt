#'documentation.py' is from cms software github(https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/test)
#the original file, 'inspectNanoFile.py' in "test" directory is modified only in lines, 291-292, because of 'dict_values' issue)
#only ROOT and PYTHON(especially, written for Python version 3) required
#'documentation.py' runs with output root file of NanoAODplus 

#usage:

#documentation of branch titles
>> python3 documentation.py --doc=DOC test.root

#two additional documentations with 'documentation.py' if being useful
>> python3 documentation.py --size=SIZE test.root
>> python3 documentation.py --json=JSON test.root

(DOC, SIZE, JSON can be any output name)

------------------------------------------------------------------------

##########################
###PYTHON MODULE ISSUES###
##########################

#if 'six' module is not found in current Python version 3, please copy the file, 'six.py' and then run it
>> python3 six.py

(in naf-cms systems, with adding six module in the system, 'documentation.py' works)

