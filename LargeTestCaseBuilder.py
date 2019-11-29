import glob, os

# When ran this file converts *.graph files downloaded from https://www.cc.gatech.edu/dimacs10/archive/delaunay.shtml
# to *.txt files with DIMACS format
# the output *.txt files are stored in the tests/ folder

def convertGraphToMaxflowFormat(fileName):

    dir_path = 'tests/'
    filepath = fileName.split('/')
    graphName = filepath[1][:-6]
    maxFlowGraphFile = open(os.path.join(dir_path, graphName + ".txt"), "w")

    with open(fileName) as fp:
        line = fp.readline()
        nums = line.split()
        numVertices = nums[0]
        maxFlowGraphFile.write("p max "  + line)
        maxFlowGraphFile.write("n 1 s\n")
        maxFlowGraphFile.write("n "+numVertices+" t\n")
        line = fp.readline()
        currVertex = 1
        while line:
            neighbors = line.split()
            for neighbor in neighbors:
                maxFlowGraphFile.write("a " + str(currVertex) + " " + neighbor+" 1\n")
            currVertex += 1
            line = fp.readline()

    maxFlowGraphFile.close()



def main():
    for file in glob.glob("raw_test_files/*.graph"):
        print(file)
        convertGraphToMaxflowFormat(file)

if __name__ == '__main__':
    main()