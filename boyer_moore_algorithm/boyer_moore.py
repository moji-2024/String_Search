import pandas as pd
class badChar_goodSufixRuls:
    def __init__(self, pattern: str, alphabet=['A', 'T', 'C', 'G']):
        self.pattern = pattern
        self.alphabet = alphabet
        self.lenPattern = len(pattern)
        # create dict of goodSuffixIndexes
        self.goodSuffixData = {}
        # first create all possible Prefix (big length to small length) for get skips
        allPrefix = [pattern[:i] for i in range(1, self.lenPattern)][::-1]
        for ind_N in range(self.lenPattern - 1, -1, -1):
            suffix = pattern[ind_N:]
            lenSuffix = len(suffix)
            errorInd = pattern.rfind(suffix, 0, ind_N)
            # The suffix does NOT appear anywhere else
            if errorInd == -1:
                # Check if part of the suffix matches the beginning (prefix) of the pattern.
                for prefix in allPrefix:
                    if prefix == suffix[lenSuffix - len(prefix):]:
                        self.goodSuffixData[suffix] = self.lenPattern - 1 - len(prefix)
                        break
                # If no overlap exists, shift the whole pattern past the mismatch.
                else:
                    self.goodSuffixData[suffix] = self.lenPattern - 1
            # The suffix appears somewhere else in the pattern
            else:
                self.goodSuffixData[suffix] = self.lenPattern - 1 - (errorInd + lenSuffix - 1) - 1
        # create table of badCharIndexes:
        # rows= characters of pattern with their indexes
        # column= characters of alphabets
        data = {}
        for ind_N in range(self.lenPattern - 1, -1, -1):
            Nucleotid = pattern[ind_N]
            for alpha_N in alphabet:
                if alpha_N == Nucleotid:
                    data.setdefault(f'{Nucleotid}_{ind_N}', []).append('-')
                    continue
                errorInd = pattern.rfind(alpha_N, 0, ind_N)
                if errorInd == -1:
                    data.setdefault(f'{Nucleotid}_{ind_N}', []).append(ind_N)
                else:
                    data.setdefault(f'{Nucleotid}_{ind_N}', []).append(ind_N - errorInd - 1)
        self.skipTable = pd.DataFrame(data, index=alphabet)

    def _skipIndexByBadCharRule(self, refSubStr: str, patternSubStr: str):
        return self.skipTable.loc[refSubStr, patternSubStr]

    def _skipIndexByGoodSuffixRule(self, GoodSuffix):
        return self.goodSuffixData[GoodSuffix]

    def findMatch(self, refGene: str):
        matches = []
        indexRef = 0
        while indexRef < len(refGene):
            refGenePortion = refGene[indexRef:indexRef + self.lenPattern]
            if len(refGenePortion) < self.lenPattern:
                break
            if refGenePortion == self.pattern:
                matches.append(indexRef)
                indexRef += self.lenPattern
                continue
            for indexPattern in range(self.lenPattern - 1, -1, -1):
                if refGenePortion[indexPattern] != self.pattern[indexPattern]:
                    # Check if current character of text is defined in pattern and alphabet
                    if refGenePortion[indexPattern] not in self.alphabet:
                        indexRef += indexPattern
                        break
                    else:
                        badCharInd = self._skipIndexByBadCharRule(refGenePortion[indexPattern],
                                                                  f'{self.pattern[indexPattern]}_{indexPattern}')
                        if indexPattern < self.lenPattern - 1:
                            goodSuffixInd = self._skipIndexByGoodSuffixRule(refGenePortion[indexPattern + 1:])
                        else:
                            goodSuffixInd = 0
                        indexRef += max(goodSuffixInd, badCharInd)
                        break
            indexRef += 1
        return matches
