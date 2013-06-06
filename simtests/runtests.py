#!/usr/bin/python
"""
    TestRunner and TestDefinition for simple full program testing
    rv july 2007

    1. every test setup is given in a separate directory
    2. the directory content is given by
          input/ run/ ref/ misc/ testdef
         - in input/ is the unmodified input data
         - run/ is the directory where input data is copied to and executed
         - ref/ is the reference output (so a copy of an older run/ directory)
         - misc/ is a directory that may contain helper scripts to analyze any errors etc.
         - testdef is a file defining the files and colums that should be compared
"""

import os;
import re;
import sys;
import time;
import socket;

class TestDefinition:
    """
        TestDefinition class
        holds data file name, information on colums that should be compare etc.
    """
    def  __init__(self, dataFileName, tagName, colums, rows, relativeError, absoluteError):
        self.dataFileName  = dataFileName
        self.tagName       = tagName
        self.colums        = colums
        self.rows          = rows
        self.relativeError = float(relativeError)
        self.absoluteError = float(absoluteError)

        if self.relativeError <= 0 and self.absoluteError < 0:
            print "the error definitions for %s are bullshit! " % tagName
            sys.exit(2)

    def compare(self, rootDirectory):
        """
            compare test values with reference
        """
        reasons = []
        dataRun = self.readXYDataFile(rootDirectory + "/run/" + self.dataFileName)
        if dataRun == -1:
            reasons.append("could not open or read data file in directory run/");
        else:
            dataRef = self.readXYDataFile(rootDirectory + "/ref/" + self.dataFileName)
            if dataRef == -1:
                reasons.append("could not open or read reference data file");
            else:
                for cidx in self.colums.split(':'):
                    ret = self.compareColumns(dataRun,dataRef,int(cidx))
                    if not ret == -1:
                        reasons.append(ret)

        if len(reasons) == 0:
            print " PASSED - %s in %s" % (self.tagName, rootDirectory)
            return 0
        else:
            print " FAILED - %s in %s" % (self.tagName, rootDirectory)
            for r in reasons:
                print "          %s" % r
            return 1

    def compareColumns(self, dataRun, dataRef, cidx):
        """
            compare column cidx data's
        """
        if len(dataRun) != len(dataRef):
            return 'length of datasets is unequal (run: %d, ref: %d)' % (len(dataRun),len(dataRef))
        maxErrorAbsolute = [-1, -1];
        maxErrorRelative = [-1, -1];

        # restrict yourself to some rows (needed e.g. for DOS, where calculation is numerically unstable)
        minRow = 0
        maxRow = len(dataRun)
        if(self.rows != -1):
            rowRestriction = self.rows.split(":");
            if len(rowRestriction) != 2:
                print("invalid row restriction string: " + str(self.rows));
                sys.exit(1)
            minRow = max(minRow,int(rowRestriction[0]))
            maxRow = min(maxRow,int(rowRestriction[1]))

            if(minRow >= maxRow):
                print("invalid row range: " + str(self.rows))
                sys.exit(1)
        ret = ""
        num_errors = 0

        for ii in range(len(dataRun)):
            if minRow <= ii and maxRow > ii:

                a = dataRun[ii][cidx]
                b = dataRef[ii][cidx]
                f = 1;
                if a != 0:
                    f = a
                elif b != 0:
                    f = b
                actualErrorRelative = abs(a - b) / abs(f);
                actualErrorAbsolute = abs(a - b);

                if self.absoluteError != -1 and actualErrorAbsolute > self.absoluteError:
                    if maxErrorAbsolute[0] == -1 or maxErrorAbsolute[1] < actualErrorAbsolute:
                        maxErrorAbsolute[0] = ii
                        maxErrorAbsolute[1] = actualErrorAbsolute
                if self.relativeError != -1 and actualErrorRelative > self.relativeError:
                    if maxErrorRelative[0] == -1 or maxErrorRelative[1] < actualErrorRelative:
                        # skip it if absolute values are smaller than the absolute error threshold
                        if abs(a) > self.absoluteError and abs(b) > self.absoluteError:
                            maxErrorRelative[0] = ii
                            maxErrorRelative[1] = actualErrorRelative

                if maxErrorRelative[0] != -1 or maxErrorAbsolute[0] != -1:
                    if num_errors == 10:
                        ret = ret + "\n          too many errors! stopping output!\n"
                    elif num_errors < 10:
                        if ret != "":
                            ret = ret + "\n          "
                        if maxErrorRelative[0] != -1:
                            ret = ret + 'column: %d, row: %d, max. relative err. of %g > %g (abs value: %g)' % (cidx,maxErrorRelative[0], maxErrorRelative[1],self.relativeError, abs(a))
                        if maxErrorAbsolute[0] != -1:
                            if len(ret) > 0 and maxErrorRelative[0] != -1:
                                ret = ret + ', '
                            ret = ret + 'column: %d, row: %d, max. absolute err. of %g > %g ' % (cidx, maxErrorAbsolute[0], maxErrorAbsolute[1],self.absoluteError)
                    num_errors += 1
                    maxErrorAbsolute = [-1, -1];
                    maxErrorRelative = [-1, -1];

        if ret != "":
            ret = ret + "\n          %d out of %d data sets failed" % (num_errors, len(dataRun));
            return ret
        else:
            return -1



    def readXYDataFile(self, fileName):
        """
            read and parse XY Data file consisting either of doubles or even complex numbers
        """
        commentProg = re.compile('^\s*#');
        doubleProg  = re.compile('([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)')
        complexProg = re.compile('\(([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?),([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\)')
        ret         = []
        try:
            f = open(fileName, 'r')
        except IOError, (errno, strerror):
            return -1

        for line in f.readlines():
            if not commentProg.search(line):
                dataLine = [];
                items = line.split(' ');
                for item in items:
                    if complexProg.search(item):
                        mobj = complexProg.match(item)
                        real = float(mobj.expand('\\1'))
                        imag = float(mobj.expand('\\5'))
                        dataLine.append(complex(real,imag))
                    elif doubleProg.search(item):
                        dataLine.append(float(item))
                ret.append(dataLine)
        if len(ret) == 0:
            return -1
        return ret;





    def __str__(self):
        return 'TestDefinition with\n  file = %s\n  colums = %s  \n  rel = %s\n  abs = %s' % (self.dataFileName,self.colums, self.relativeError, self.absoluteError)


class TestRunner:
    """
        TestRunner Class Executing and Evalulating Tests
    """
    tests = []
    timings = []
    def __init__(self,tdkpExecutable):
        print " starting TestRunner"
        mypath = os.popen('pwd').read();
        mypath = mypath.strip('\n');
        self.tdkpexec = "%s/%s" % (mypath, tdkpExecutable)

        if not os.access(self.tdkpexec, os.R_OK | os.X_OK):
            print " *error*: executable " + self.tdkpexec + " does not exist or can not be executed"
            sys.exit(1)


    def findTests(self, directory):
        """
            searches current directory tree for test definitions
        """
        self.tests = []
        out = os.popen("find " + directory + " -name \"testdef\"").read();
        for ss in out.split():
            self.tests.append(ss[0:-8]);
        print (" found %d test(s)" % (len(self.tests)))

    def memcheck(self):
        """
            use valgrind to perform memchecks
        """
        for testDirectory in self.tests:
            print " -------------------------------------------------------- "
            print "  doing a valgrind in %s " % testDirectory
            print " -------------------------------------------------------- "
            myPipe = os.popen("cd " + testDirectory + "/run;valgrind --leak-check=full --leak-resolution=high --log-file-exactly=valgrind.out " + self.tdkpexec + " run");
            out = myPipe.read()
            ret = myPipe.close()
            if ret != None:
                print " ERROR  - calling the script run in %s failed with error-code: %s" % (testDirectory, ret)
                print " program output is: "
                print out

            f = open(testDirectory + "/run/valgrind.out", 'r');
            self.analyzeValgrindOutput(f.read());
            f.close();

    def analyzeValgrindOutput(self, string):
        """
            analyze valgrind file output
        """
        print ""

    def run(self):
        """
            run tests
        """
        errors    = 0
        testcount = 0
        self.timings = []
        for testDirectory in self.tests:
            # check if test already passed
            if os.path.exists(testDirectory + '/run/PASSED'):
                print " test " + testDirectory + " already passed last time"
                self.timings.append(0)
            else:
                # read file
                f = open(testDirectory + '/testdef', 'r');
                definitionString = f.read();
                f.close();
                testDef = self.readTestDefinitions(definitionString)
                self.setupSimulation(testDirectory)
                self.timings.append(self.runSimulation(testDirectory))
                errors    += self.compareData(testDirectory,testDef)
                testcount += len(testDef)
        self.showStats(errors,testcount)
        self.storeBenchmarks()

    def storeBenchmarks(self):
        """
			store benchmark timings in benchmark.dat
        """
        f = open("benchmarks.dat", 'w');
        f.write("# -----------------------------------------------\n");
        f.write("# tdkp benchmarks evaluated on " + time.strftime("%d.%h.%Y %H:%M:%S", time.localtime()) + "\n");
        f.write("# on host " + socket.gethostname() + "\n");
        f.write("# -----------------------------------------------\n");

        for ii in range(len(self.tests)):
            testDirectory = self.tests[ii];
            tm = self.timings[ii];
            f.write(testDirectory + "  " + str(tm) + "\n");
        f.close();



    def evaluate_only(self):
        """
            don't calculate tests, but just evaluate older runs
        """
        errors    = 0
        testcount = 0
        for testDirectory in self.tests:
            # read file
            f = open(testDirectory + '/testdef', 'r');
            definitionString = f.read();
            f.close();
            testDef = self.readTestDefinitions(definitionString)
            errors    += self.compareData(testDirectory,testDef)
            testcount += len(testDef)
        self.showStats(errors,testcount)

    def showStats(self,errors,testscount):
        if testscount > 0:
            ratio = float(errors) / float(testscount);
        else:
            ratio = 0
        if ratio == 0:
            msg = "a perfect day!"
        elif ratio < 0.05:
            msg = "almost o.k., i will manage that!"
        elif ratio < 0.1:
            msg = "well ... that will probably need you some time to fix it"
        elif ratio < 0.2:
            msg = "ouuuhhu ... that does not look good!"
        elif ratio < 0.5:
            msg = "hmm ... maybe you will have cancel your further plans ..."
        elif ratio < 0.8:
            msg = "did you ever thought about getting an easier job?"
        else:
            msg = "boahh ... you broke it completely!"


        print "\n\n\n ------------------------------------------------------ "
        print " succedded in %d of %d tests: %s" % (testscount - errors, testscount, msg)
        print " ------------------------------------------------------ "



    def cleanUp(self):
        """
            delete contents in directory run
        """
        for testDirectory in self.tests:
            os.system('rm -f ' + testDirectory + '/run/*');


    def readTestDefinitions(self,definitionString):
        """
            read test definitions and create TestDefinition objects
        """
        tests    = [];
        outer    = definitionString.lstrip(' \n').split('}');
        tagFind  = re.compile('([a-z]+)\s*=\s*(.+)');
        for inner in outer:
            arguments = inner.split(",")
            if len(arguments) > 1:
                myTags = {"rel" : -1, "abs" : -1, "rows" : -1}
                for ii in range(len(arguments)):
                    arguments[ii] = arguments[ii].lstrip("{ \n")
                    mobj = tagFind.match(arguments[ii])
                    if mobj != None:
                        tag = mobj.expand('\\1');
                        arg = mobj.expand('\\2');
                        myTags[tag] = arg;
                tests.append(TestDefinition(myTags["file"], myTags["tag"], myTags["col"], myTags["rows"], myTags["rel"], myTags["abs"]));

        return tests;



    def setupSimulation(self, testDirectory):
        """
            delete contents in run directory and copy files from input to run
        """
        if len(os.listdir(testDirectory + '/run')) > 0:
            os.system('rm  ' + testDirectory + '/run/*');
        os.system('cp ' + testDirectory + '/input/* ' + testDirectory + '/run');

    def runSimulation(self, testDirectory):
        """
            run given simulation in given directory
        """
        print " -------------------------------------------------------- "
        print "  calculating simulation in %s " % testDirectory
        print " -------------------------------------------------------- "
        startingTime = time.time();
        myPipe = os.popen("cd " + testDirectory + "/run;" + self.tdkpexec + " run");
        out = myPipe.read()
        ret = myPipe.close()
        endTime = time.time();
        # store output into file
        f = open(testDirectory + "/run/run.output", 'w')
        f.write(out)
        f.close()


        if ret != None:
            print " ERROR  - calling the script run in %s failed with error-code: %s" % (testDirectory, ret)
            print " program output is: "
            print out

        return endTime - startingTime;


    def compareData(self, testDirectory, testDefinitions):
        """
            compare run and reference data
        """
        errors = 0
        for testDef in testDefinitions:
            errors += testDef.compare(testDirectory)
        # create file passed to indicate that results are o.k.
        if errors == 0:
            os.system('touch ' + testDirectory + '/run/PASSED');
        return errors












# -----------------------------------------------------------
# main routine
# -----------------------------------------------------------
start_time = time.time()
command = ""

if len(sys.argv) > 1:
    command = sys.argv[1]


runner  = TestRunner("../tdkpshell.bin")
if command == 'all':
    runner.findTests('.');
    runner.run();
elif command == 'single':
    if len(sys.argv) > 2:
        for ii in range(2,len(sys.argv)):
            directory = sys.argv[ii];
            runner.findTests(directory);
            runner.run();
    else:
        print " *error*: no directories defined!"
        sys.exit(1);
elif command == 'clean':
    directory = "."
    if len(sys.argv) > 2:
        directory = sys.argv[2]
    runner.findTests(directory)
    runner.cleanUp()
elif command == 'evaluate':
    directory = "."
    if len(sys.argv) > 2:
      directory = sys.argv[2]
    runner.findTests(directory)
    runner.evaluate_only()
elif command == 'update':
    if len(sys.argv) == 4:
        source = sys.argv[2]
        target = sys.argv[3] + "/" + source + "/ref"
        if not os.path.exists(source):
            print " *error*: source directory " + source + " does not exist"
            sys.exit(2)
        if not os.path.exists(target):
            print " *error*: target directory " + target + " does not exist"
            sys.exit(2)

        existingFiles = os.listdir(target)
        if len(existingFiles) > 0:
            os.system("rm " + target + "/*")
            print " deleted old reference files"
        os.system("cp " + source + "/run/* " + target)
	# delete any cache file in target dir
        os.system("rm -f " + target + "/cache_*");
        # delete any binary file in target dir (we don't save binary data
        # because that would take too much space!)
        os.system("rm -f " + target + "/*.bin");
        existingFiles = os.listdir(target)
        print " updated " + target + ". directory has now the following files"
        for file in existingFiles:
            print "   " + file

    else:
        print " Usage: runtests.py update <source> <target>"
elif command == 'memcheck':
    # -----------------------------------------
    # unfunctional at the moment
    # ----------------------------------------
    if sys.argv[2] == 'all':
        runner.findTests('.');
        runner.memcheck();
    else:
        for ii in range(2,len(sys.argv)):
            directory = sys.argv[ii];
            runner.findTests(directory);
            runner.memcheck();
else:
    print " Usage: runtests.py <command> [COMMAND OPTIONS]\n"
    print "Available commands:"
    print "  all                          run all available tests"
    print "  single <dirname> [<dir>..]   run tests in directories named <dirname>"
    print "  clean [<dirname>]            clean up evaluation data [in directory <dir> if defined]"
    print "  evaluate [<dirname>]         only evaluate available results against reference"
    print "  update <source> <target>     store run data in test dir <source> to target "
    print "                               root directory <target>"
    print "\nAnother piece of software    provided by r.v. ;-)"
    sys.exit(0)


end_time = time.time()

print("runtests.py was running for %d seconds" % int(end_time - start_time))
