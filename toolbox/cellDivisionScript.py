from toolbox import tissueSample
from toolbox import cellDivision
import random

def main():
    sample = tissueSample.Sample(configDir = "samples/", simulationTime = 500, tissueType = "periodic")
    crossBoundary = True
    while crossBoundary:
        cellID = random.choice(list(sample.cells_.keys()))
        crossBoundary = sample.cells_[cellID].crossBoundary_
    sample.cells_[cellID].is_mother_ = True
    print("Mother cell ID: ", cellID)
    cellDivision.dumpCellVtk(sample, cellID)
    sample = cellDivision.evaluatePostDivisionTopology(sample, cellID)
    cellDivision.dumpSample(sample)
    return

if __name__ == "__main__":
    main()