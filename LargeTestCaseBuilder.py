import glob, os

def convertGraphToMaxflowFormat(fileName):

    graphName = fileName[:-6]
    maxFlowGraphFile = open(graphName + ".txt", "w")

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
    os.chdir("tests")
    for file in glob.glob("*.graph"):
        print(file)
        convertGraphToMaxflowFormat(file)

if __name__ == '__main__':
    main()