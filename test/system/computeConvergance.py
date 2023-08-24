from math import *

import sys
if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3!")

def convergeRateRatio(error, refRatio):
    assert len(error) == len(refRatio)+1
    assert len(error) >= 2
    assert all(r>0 for r in refRatio)
    return [log(error[i+1]/error[i])/log(1/refRatio[i]) if error[i] != 0 else 0 for i in range(len(error)-1)]

def convergeRateDOF(error, dof, dim):
    assert dim > 0
    assert len(dof) >= 2
    assert all(dof[i] < dof[i+1] for i in range(len(dof)-1))
    return dim*convergeRateRatio(error, [(dof[i+1]/dof[i])**(1./dim) for i in range(len(dof)-1)])

def convergeReport(errors, dofs, dim, componentNames, errorNames):
    # process each component separate
    for c in range(len(componentNames)):
        print("Errors for " + componentNames[c] + ": ")
        printTitles = ["DoF"]
        printArr = [dofs]
        for e in range(len(errorNames)):
            er = [errors[r][c][e] for r in range(len(dofs))]
            rate = convergeRateDOF(er, dofs, dim)
            printTitles.append(errorNames[e])
            printArr.append(er)
            printTitles.append("Rate(" + errorNames[e] + ")")
            printArr.append([nan] + rate)
        # Print the titles 
        s=''
        for row in range(len(printTitles)):
           s+=format(printTitles[row], '<11s')
        print(s)
        # Print the table
        for row in range(len(printArr[0])):
            s=''
            for col in range(len(printArr)):
                if isfinite(printArr[col][row]):
                    s+=format(printArr[col][row], '<11.3e')
                else:
                    s+=format('---', '^11s')
            print(s)

def convergeReportFiles(errors, dofs, dim, componentNames, errorNames):
    # process each component separate
    for c in range(len(componentNames)):
        printTitles = ["DoF"]
        root = lambda x, dim=dim: x^(1./dim)
        printArr = [list(map(sqrt,dofs))]
        for e in range(len(errorNames)):
            er = [errors[r][c][e] for r in range(len(dofs))]
            rate = convergeRateDOF(er, dofs, dim)
            printTitles.append(errorNames[e])
            printArr.append(er)
        # Print the titles
        filename = componentNames[c] + ".dat"
        output_file = open(filename,'w')
        s=''
        for row in range(len(printTitles)):
           s+=printTitles[row] + "\t"
        output_file.write(s + "\n")
        # Print the table
        for row in range(len(printArr[0])):
            s=''
            for col in range(len(printArr)):
                s+=format(printArr[col][row], 'e') + "\t"
            output_file.write(s + "\n")
        output_file.close()
