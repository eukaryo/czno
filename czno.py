
import numpy as np
import heapq





#from clote_misc import energyOfStr, basePairList

import os,tempfile,subprocess


def energyOfStr(rna,secStr):
    with open(".tmp1",'w') as file:
        file.write(rna+"\n")
        file.write(secStr)
    subprocess.call('RNAeval < %s > %s' % (".tmp1", ".tmp2"), shell=True)
    with open(".tmp2",'r') as file:
        rna0 = file.readline().strip()
        assert( rna0 == rna )
        energy = file.readline().split()[-1]
        energy = energy.replace('(','')
        energy = energy.replace(')','')
        energy = float(energy)
        return energy

  #tmpfilename = '.tmp'
  #file        = open(tmpfilename,'w')
  #file.write(rna+"\n")
  #file.write(secStr)
  #file.close()
  #cmd  = 'RNAeval < %s' % tmpfilename
  #file = os.popen(cmd)
  #rna0 = file.readline().strip()
  #assert( rna0 == rna )
  #energy = file.readline().split()[-1]
  #energy = energy.replace('(','')
  #energy = energy.replace(')','')
  #energy = float(energy)
  #return energy



def basePairList(secStr):
    L = []; stack = []
    for i in range(len(secStr)):
        if secStr[i]=="(":
            stack.append(i)
        elif secStr[i]==")":
            x = stack.pop() #return and pop end of stack 
            L.append( (x,i) ) 
    return L






def isExclusive(bp1, bp2):
    #引数は相異なる2つの塩基対とする。
    #2つの塩基対が同時に存在し得ないならTrue

    assert(bp1 != bp2)
    #同じ塩基とくっつくならTrue
    if bp1[0] == bp2[0]: return True
    if bp1[0] == bp2[1]: return True
    if bp1[1] == bp2[0]: return True
    if bp1[1] == bp2[1]: return True

    #pseudo-knotの位置関係にあるならTrue
    if bp1[0] < bp2[0] < bp1[1] < bp2[1]: return True
    if bp2[0] < bp1[0] < bp2[1] < bp1[1]: return True

    return False

def makeDotBracket(structure, bp_set1, diff_list):
    for d in [s for s in diff_list if s in bp_set1]:
        structure = structure[:d[0]]+"."+structure[d[0]+1:d[1]]+"."+structure[d[1]+1:]
    for d in [s for s in diff_list if s not in bp_set1]:
        structure = structure[:d[0]]+"("+structure[d[0]+1:d[1]]+")"+structure[d[1]+1:]
    return structure


def getMinBarrierPathDijkstra(sequence, structure1, structure2):

    #引数の構造はDot-Bracket記法だが、これを塩基対のリストに変換する。
    bp_set1 = set(basePairList(structure1))
    bp_set2 = set(basePairList(structure2))

    #塩基対のリストの相異なるものを抽出する。これがMAX_DIST以内でなければ諦める
    diff_pos_set = bp_set1 ^ bp_set2
    diff_pos_list = list(diff_pos_set)
    MAX_DIST = 25
    H = len(diff_pos_set)
    if H > MAX_DIST: return None


    #prerequisite_flagというテーブルを作る。
    #prerequisite_flag[x]のy番目のビットが1になっている⇔
    #diff_pos_list[x]を変化させるとき、予めdiff_pos_list[y]を変化させていなければならない。
    #
    #(1)(i,j)が塩基対を組んでいるとき、(i,k)や(k,j)が同時に組むことは出来ない
	#(2)pseudoknot禁止
    #これらを満たすためにprerequisite_flagが要請される。
    prerequisite_flag = np.zeros(H, dtype='uint64')
    for i in range(H):
        if diff_pos_list[i] not in bp_set1: continue
        for j in range(H):
            if diff_pos_list[j] not in bp_set2: continue

            #iを踏んでからjを踏む必要がある場合、obstacle[j] |= 1 << i とする。
			#そのようなケースは、(max-loop制約を無視すると)
			#(1)スタート構造でiが組まれていてjは組まれていない
			#(2)ゴール構造ではiは組まれておらずjが組まれている
			#ようなケースに限られる。
			#(max-loop制約においては、先に外すことがinvalidになるケースがあるが、今回はRNAevalを使うので問題ない。)

            if isExclusive(diff_pos_list[i], diff_pos_list[j]):
                prerequisite_flag[j] = np.bitwise_or(prerequisite_flag[j], np.uint64(2**i))
                
    #result[x]はスタート構造から構造xまでの部分解
    result = np.full(2**H, np.inf)

    #energy[x]は構造xのエネルギーの値
    energy = np.full(2**H, np.inf)

    #Dijkstraは二分ヒープで、要素を小さい順に出す。
    #Dijkstraの要素は(部分解、-時刻、構造のインデックス)
    #-時刻を入れるのは、部分解が同じだとLIFOで出て欲しいから。
    Dijkstra = []
    count = 0

    energy[0] = energyOfStr(sequence, structure1)
    heapq.heappush(Dijkstra, (energy[0], count, np.uint64(0)))
    while len(Dijkstra) > 0:
        x = heapq.heappop(Dijkstra)
        index = np.uint64(x[2])
        if result[index] != np.inf: continue
        result[index] = x[0]
        if index == np.uint64(2 ** H - 1): break
        for i in range(H):
            if np.bitwise_and(index, np.uint64(2**i)) != np.uint64(0): continue
            if np.bitwise_and(index, prerequisite_flag[i]) != prerequisite_flag[i]: continue
            new_index = np.bitwise_or(index, np.uint64(2**i))
            if result[new_index] != np.inf: continue
            if energy[new_index] == np.inf:
                new_structure = makeDotBracket(structure1, bp_set1, [diff_pos_list[x] for x in range(H) if np.bitwise_and(new_index, np.uint64(2**x)) != np.uint64(0)])
                energy[new_index] = energyOfStr(sequence, new_structure)
            new_result = max(result[index], energy[new_index])
            count -= 1
            heapq.heappush(Dijkstra, (new_result, count, new_index))


    #トレースバック
    answer = [structure2]
    index = np.uint64(2 ** H - 1)
    while index > np.uint64(0):
        flag = False
        for i in range(H):
            if np.bitwise_and(index, np.uint64(2**i)) == np.uint64(0): continue
            back_index = np.bitwise_xor(index, np.uint64(2**i))
            if np.bitwise_and(back_index, prerequisite_flag[i]) != prerequisite_flag[i]: continue
            if result[index] != max(energy[index], result[back_index]): continue
            answer.append(makeDotBracket(structure1, bp_set1, [diff_pos_list[x] for x in range(H) if np.bitwise_and(back_index, np.uint64(2**x)) != np.uint64(0)]))
            index = back_index
            flag = True
            break
        if flag == False:
            answer.append("fail")
            break
    answer.reverse()
    return answer

if __name__ == '__main__':
    
    seq1 = "CCCCAAAACCCCAAAAGGGGAAAAGGGG"
    str1 = "((((............))))........"
    str2 = "........((((............))))"
    ans = getMinBarrierPathDijkstra(seq1, str1, str2)

    print(seq1)
    for s in ans:
        print(s)
