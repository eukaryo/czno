
import re, random
from czno import getMinBarrierPathDijkstra, basePairList, energyOfStr
from previous_direct_methods import MorganHiggs1998GreedyDirect, Voss2004GreedyDirect

def HammingDistance(structure1, structure2):
    bp1 = set(basePairList(structure1))
    bp2 = set(basePairList(structure2))
    return len(bp1 ^ bp2)

def InterpretDatasetFile(filename):

    dataset = []
    seq = ""
    structures = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if re.fullmatch(r"^[ACGU]+$", line) is not None:
                if len(structures) >= 1:
                    dataset.append([seq, structures])
                seq = line
                structures = []
            elif re.fullmatch(r"^[(.)]+$", line) is not None:
                structures.append(line)
            else:
               assert False

    if len(structures) >= 1:
        dataset.append([seq, structures])
    return dataset

def PathwayToBarrier(sequence, pathway):
    return max([energyOfStr(sequence, s) for s in pathway])
    

def DirectPathSingleExperiment(sequence, structure1, structure2):

    result = []

    pathway1 = getMinBarrierPathDijkstra(sequence, structure1, structure2)
    result.append(PathwayToBarrier(sequence, pathway1))

    pathway2 = MorganHiggs1998GreedyDirect(sequence, structure1, structure2)
    result.append(PathwayToBarrier(sequence, pathway2))

    pathway3 = Voss2004GreedyDirect(sequence, structure1, structure2)
    result.append(PathwayToBarrier(sequence, pathway3))

    best_value4 = float("inf")
    for i in range(10):
        pathway4 = Voss2004GreedyDirect(sequence, structure1, structure2, 10, i)
        best_value4 = min(best_value4, PathwayToBarrier(sequence, pathway4))
    result.append(best_value4)

    return result


if __name__ == '__main__':
    dataset = InterpretDatasetFile("s151-localminima-dataset.txt")

    for data in dataset:
        starts = random.sample(data[1], min(1, len(data[1])))
        for i in range(len(starts)):
            g1 = list(set(data[1]) - set(starts)) + starts[i+1:]
            g2 = [HammingDistance(starts[i],g1[x]) for x in range(len(g1))]
            g3 = [g1[x] for x in range(len(g1)) if 5 <= g2[x] and g2[x] <= 15]
            goals = random.sample(g3, min(1, len(g3)))
            for j in range(len(goals)):
                hamdist = HammingDistance(starts[i], goals[j])
                if 5 <= hamdist and hamdist <= 15:
                    result = DirectPathSingleExperiment(data[0], starts[i], goals[j])
                    print(str(len(data[0]))+" "+str(hamdist)+" "+" ".join([str(x) for x in result]))

