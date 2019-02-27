
import re
from czno import getMinBarrierPathDijkstra, basePairList, energyOfStr

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
    

def SingleExperiment(sequence, structure1, structure2):

    result = []

    pathway = getMinBarrierPathDijkstra(sequence, structure1, structure2)
    result.append(max([energyOfStr(sequence, s) for s in pathway]))

    #TODO: ここで、先行手法をsubprocessで叩いて結果を得てエネルギーバリアを求める

    return result


if __name__ == '__main__':
    dataset = InterpretDatasetFile("s151-localminima-dataset.txt")
    for data in dataset:
        for i in range(len(data[1])):
            for j in range(i+1, len(data[1])):
                hamdist = HammingDistance(data[1][i], data[1][j])
                if 5 <= hamdist and hamdist <= 20:
                    result = SingleExperiment(data[0], data[1][i], data[1][j])
                    print(str(len(data[0]))+" "+str(hamdist)+" "+" ".join([str(x) for x in result]))

