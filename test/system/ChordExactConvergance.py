from computeConvergance import *

def findFloatsInString(string):
    f = []
    for t in string.split():
        try:
            f.append(float(t))
        except ValueError:
            pass
    return f
                
def extractErrors_ChordExact(filename):
    with open(filename) as file:
        errorsReported=False
        componentNames = []
        errorNames = []
        errors = []
        for line in file:
            if "-> Error norms at time" in line:
                print("Found start of solution error at time " +
                      str(findFloatsInString(line)[0]) +
                      " from " + filename)
                errorsReported=True
            if "Component" in line and errorsReported:
                componentNames.append(line.split()[-1])
                errors.append([])
            if "(error)" in line and errorsReported:
                errorName = line.split()[0]
                compNum = len(componentNames)-1
                if errorName in errorNames:
                    errors[compNum].insert(errorNames.index(errorName),
                                           findFloatsInString(line)[0])
                else:
                    errorNames.append(errorName)
                    errors[compNum].append(findFloatsInString(line)[0])
        return [errors, componentNames, errorNames]

def convergenceFromChordOutput(files, gridSize, dim, output):
    errors=[]
    componentNames=''
    errorNames=''
    for f in files:
        [er, cn, en] = extractErrors_ChordExact(f)
        errors.append(er)
        if componentNames == '':
            componentNames = cn
            errorNames = en
        assert componentNames == cn
        assert errorNames == en
    output(errors, gridSize, dim, componentNames, errorNames)

files = ['OoA/pout_16', 'OoA/pout_32', 'OoA/pout_64', 'OoA/pout_128']

dim = 2
gridSize = [16**dim, 32**dim, 64**dim, 128**dim]

# output type needs to be either convergeReport or convergeReportFiles
output_type = convergeReport
convergenceFromChordOutput(files, gridSize, dim, output_type)
