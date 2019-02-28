
import random
from czno import energyOfStr, basePairList, isExclusive

def basePairSetToDotNotation(sequence, bp):
    answer = ["." for _ in sequence]
    for x in bp:
        assert x[0] < x[1]
        answer[x[0]]="("
        answer[x[1]]=")"
    return "".join(answer)



def MorganHiggs1998GreedyDirect(sequence, structure1, structure2):
    A = set(basePairList(structure1))
    B = set(basePairList(structure2))
    answer = [structure1]
    while A != B:
        NeedToClose = sorted(list(B - A))
        NeedToOpen = sorted(list(A - B))
        if len(NeedToClose) >= 1:
            NextClose = min([(sum([isExclusive(x,y) for y in NeedToOpen]), x) for x in NeedToClose])[1]
            NextOpens = [y for y in NeedToOpen if isExclusive(NextClose, y)]
            for x in NextOpens:
                A.remove(x)
                answer.append(basePairSetToDotNotation(sequence, A))
            A.add(NextClose)
            answer.append(basePairSetToDotNotation(sequence, A))
        else:
            for x in NeedToOpen:
                A.remove(x)
                answer.append(basePairSetToDotNotation(sequence, A))
    return answer

def Voss2004GreedyDirect(sequence, structure1, structure2, k = 1, seed = 12345):
    random.seed(seed)
    A = set(basePairList(structure1))
    B = set(basePairList(structure2))
    answer = [structure1]
    while A != B:
        NeedToClose = sorted(list(B - A))
        NeedToOpen = sorted(list(A - B))
        PossibleClose = [x for x in NeedToClose if sum([isExclusive(x,y) for y in NeedToOpen]) == 0]
        PossibleChange = NeedToOpen + PossibleClose
        PossibleNextState = [basePairSetToDotNotation(sequence, A^set([x])) for x in PossibleChange]
        PossibleNextStateAndEnergy = [(energyOfStr(sequence, x), x) for x in PossibleNextState]
        NextState = random.choice(sorted(PossibleNextStateAndEnergy)[0:k])[1]
        answer.append(NextState)
    return answer