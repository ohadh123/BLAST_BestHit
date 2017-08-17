# BLAST_BestHit Program
# Created by: Ohad Koronyo, Vidya Venkatraman
# Van Eyk Laboratory, August 2017

import sys
import csv
import getopt
import copy
import string
import math
import itertools

# Read input blast+ file
def readTabular(tabularFile):
    matchDict = {}
    keyOrder = []
    inputDatabits = []
    
    try:
        with open(tabularFile,"rU") as txtFile:
            reader = csv.reader(txtFile, delimiter = '\t')
            counter = 0
            insideCounter = 0
            largeCounter = 0
            for line in reader:
                counter = counter + 1
                if '#' not in line[0][0]:
                    key = line[0] #Each key represents a query
                    key = key + "~~^$%@^&@*`&<>" + str(largeCounter)
                    keyOrder.append(key)
                    if matchDict.has_key(key): #Each key maps to a 2-D array consisting subject matches
                        matchDict[key].append(line)
                    else:
                        matchDict[key] = [line]
                
                elif '#' in line[0][0] and counter == 4:
                    toParse = line[0][10:]
                    inputDatabits = toParse.split(", ")
                
                if '#' in line[0][0]:
                    insideCounter = insideCounter + 1
                    if insideCounter == 5:
                        insideCounter = 0
                        largeCounter = largeCounter + 1
    except:
        print >> sys.stderr, "Could not find/read file. Please re-check input txt file for meeting all preconditions"
        exit(1)
    
    checked = []
    for e in keyOrder:
        if e not in checked:
            checked.append(e)
        
    txtFile.close()

    return matchDict, checked, inputDatabits

# Apply Defaults filter/User filter
def applyUserFilters(arrayOfMatch, inputDatabits):
    array = arrayOfMatch
    if not array:
        return array
    filteredArray = []
    unwantedArr = []
    
    for match in array:
        try:    
            if float(match[inputDatabits.index("% identity")]) < userPreferences.get('pident',75):
                unwantedArr.append(match)
        except:
            print >> sys.stderr, "Error: Missing Percent Identity (% identity); ignoring user input for -p threshhold"

        try:
            if ((float(match[inputDatabits.index("alignment length")])-float(match[inputDatabits.index("gaps")])) / float(match[inputDatabits.index("query length")])) *100 < userPreferences.get('alength',75):
                unwantedArr.append(match)
        except:
            print >> sys.stderr, "Error: Missing Alignment length (alignment length) or Query Length (query length); ignoring user input for -a threshhold"

        try:
            if float(match[inputDatabits.index("evalue")]) > userPreferences.get('evalue',10):
                unwantedArr.append(match)
        except:
            print >> sys.stderr, "Error: Missing E-Value (evalue); ignoring user input for -e threshhold"
            
    
    for match in array:
        if match not in unwantedArr:
            filteredArray.append(match)
    return filteredArray

# Apply Identity and Coverage filter
def applyIdentityFilter(arrayOfMatch, inputDatabits, softOrHard):
    array = arrayOfMatch
    if not array:
        return array
    
    maxALQLI = 0
    maxReviewedIdentity = 0

    for match in array:
        try:
            alignmentLength = float(match[inputDatabits.index("alignment length")]) - float(match[inputDatabits.index("gaps")])
        except:
            print >> sys.stderr, ("Error: Missing Alignment length (alignment length) or Number of Gaps (gaps), skipping identity filter")
            return array
        try:
            queryLength = float(match[inputDatabits.index("query length")])  
        except:
            print >> sys.stderr, ("Error: Missing Query Length (query length), skipping identity filter")
            return array
        try:
            percentIdentity = float(match[inputDatabits.index("% identity")])  
        except:
            print >> sys.stderr, ("Error: Missing Percent Identity (% identity), skipping identity filter")
            return array
        
        try:
            calculatedALQLI = alignmentLength/queryLength * percentIdentity
            calculatedALQLI = round(calculatedALQLI,2)
        except:
            print >> sys.stderr, ("Error: in identity calculation; skipping identity filter")
            return array

        if calculatedALQLI > maxALQLI:
            maxALQLI = calculatedALQLI
        try:
            if match[inputDatabits.index("subject id")][0:match[inputDatabits.index("subject id")].find('|')] == 'sp' and calculatedALQLI > maxReviewedIdentity:
                maxReviewedIdentity = calculatedALQLI
        except:
            print >> sys.stderr, ("Error: Missing Subject ID (subject id), skipping identity filter")
            return array
        if softOrHard == 'soft':
            match.append(calculatedALQLI)
    
    filteredArray = []
    
    if softOrHard == 'soft':
        for match in array:
            threshhold = ((maxALQLI - match[-1])/(maxALQLI))*100
            if (float(match[-1]) == maxALQLI or threshhold < 20):
                if softOrHard == 'soft':
                    match.append(maxReviewedIdentity)
                filteredArray.append(match)
    else:
        for match in array:
            if (float(match[-1]) == maxALQLI): 
                if softOrHard == 'soft':
                    match.append(maxReviewedIdentity)
                filteredArray.append(match)

    return filteredArray

# Apply Review Status filter
def applyReviewFilter(arrayOfMatches, inputDatabits):
    array = arrayOfMatches
    if not array:
        return array
    filteredArray = []
    filteredURArray = []
    
    for match in arrayOfMatches:
        try:
            reviewedOrUnreviewed = match[inputDatabits.index("subject id")][0:match[inputDatabits.index("subject id")].find('|')]
        except:
            print >> sys.stderr, ("Error: Missing Subject ID (subject id), skipping review filter")
            return array

        try:
            mismatchCalculation = (float(userPreferences.get('mismatch',0))/float(match[inputDatabits.index("query length")]))*100
        except:
            print >> sys.stderr, ("Error: Missing Query Length (query length), skipping review filter")
            return array

        try:
            mismatchCalculation = match[-1] + (math.ceil(mismatchCalculation/100)*100)
            if reviewedOrUnreviewed == 'sp':
                match.pop()
                filteredArray.append(match)
            elif reviewedOrUnreviewed == 'tr' and (match[-2] > mismatchCalculation):
                match.pop()
                filteredURArray.append(match)
        except:
            print >> sys.stderr, ("Error: in identity comparison; skipping review filter")
            return array
    
    if filteredURArray:
        return filteredURArray
    elif filteredArray:
        return filteredArray
    else:
        return array

# Apply Bitscore filter
def applyBitscoreFilter(arrayOfMatches, inputDatabits, softOrHard):
    array = arrayOfMatches
    if not array:
        return array
        
    maxBitscore = -1;
    for match in array:
        try: 
            bitscore = float(match[inputDatabits.index("bit score")])  
        except:
            print >> sys.stderr, ("Error: Missing Bitscore (bit score), skipping bitscore filter")
            return array
            
        if bitscore > maxBitscore:
            maxBitscore = bitscore
            
    filteredArray = []
    if softOrHard == 'soft':
        for match in array:
            threshhold = ((maxBitscore - float(match[inputDatabits.index("bit score")])  )/(maxBitscore))*100
            if float(match[inputDatabits.index("bit score")]) == maxBitscore or threshhold < 20: 
                filteredArray.append(match)   
    else:
        for match in array:
            if float(match[inputDatabits.index("bit score")]) == maxBitscore: 
                filteredArray.append(match)   
    
    return filteredArray

# Apply Evalue filter
def applyEvalueFilter(arrayOfMatches, inputDatabits, softOrHard):
    array = arrayOfMatches
    if not array:
        return array
        
    minEvalue = float("inf");
    
    for match in array:
        try:
            evalue = float(match[inputDatabits.index("evalue")])
        except:
            print >> sys.stderr, ("Error: Missing E-Value (evalue), skipping evalue filter")
            return array
        
        if evalue < minEvalue: 
            minEvalue = evalue
    
    if minEvalue < 0.0000001:
        threshholdValue = 10000000000
    elif minEvalue > 0.0000001 and minEvalue < 0.1:
        threshholdValue = 10000000
    else:
        threshholdValue = 10000
        
    filteredArray = []
    
    if softOrHard == 'soft':
        for match in array:
            threshhold = abs(((minEvalue - float(match[inputDatabits.index("evalue")]))/(minEvalue))*100)
            if float(match[inputDatabits.index("evalue")]) == minEvalue or threshhold < threshholdValue:
                filteredArray.append(match)
    else:
        for match in array:
            if float(match[inputDatabits.index("evalue")]) == minEvalue:  
                filteredArray.append(match)   
    
    return filteredArray

# Apply Isoform filter
def applyIsoformFilter(arrayOfMatches, inputDatabits):
    array = arrayOfMatches
    if not array:
        return array
    filteredArray = []
    filteredArrayFinal = []
    dictionary = {}
    fullDictionary = {}
    isoformCount = {}

    for match in array:
        subjAccess = ''
        try:
            matchIsoform = match[inputDatabits.index("subject acc.")] #Retrieve subject accession
        except:
            print >> sys.stderr, ("Error: Missing Subject Accessifier (subject acc.), skipping isoform filter")
            return array
        
        if '-' in matchIsoform: #subjAccess is everything before the dash
            subjAccess = matchIsoform[0:matchIsoform.find('-')]  
        else:
            subjAccess = matchIsoform
        try:
            if not isoformCount[subjAccess]:
                isoformCount[subjAccess] = 0
        except:
            isoformCount[subjAccess] = 0 #isoformCount of the subjAccess is initialized to 0

        if '-' in matchIsoform: #If isoform, isoformCount for the subjAccess is incremented by 1
            isoformCount[subjAccess] += 1
            
        if dictionary.has_key(subjAccess): #Add matchIsoform under subjAccess in the dict and match under subjAccess in fullDictionary
            dictionary[subjAccess].append(matchIsoform)
            fullDictionary[subjAccess].append(match)
        else:
            dictionary[subjAccess] = [matchIsoform]
            fullDictionary[subjAccess] = [match]

    for key in dictionary:
        if '-' in key:
            sbjA = key[0:key.find('-')]
        else:
            sbjA = key
        if sbjA in dictionary[key]:
            dictionary[key] = [sbjA]
        else:
            fullDictionary[key] = applyEvalueFilter(fullDictionary[key], inputDatabits,'hard')
            helperArray = []
            for match in fullDictionary[key]:
                helperArray.append(match[inputDatabits.index("subject acc.")])
            dictionary[key] = helperArray
            dictionary[key] = sorted(dictionary[key])
            while len(dictionary[key]) != 1:
                dictionary[key].pop()


    try:
        for match in array:
            if [match[inputDatabits.index("subject acc.")]] in dictionary.values():
            
                matchIsoform = match[inputDatabits.index("subject acc.")]

                if '-' in matchIsoform and (isoformCount[matchIsoform[0:matchIsoform.find('-')]] > 1 or match[-1] < 100):
                    match[inputDatabits.index("subject acc.")] = matchIsoform[0:matchIsoform.find('-')]
                    title = match[inputDatabits.index("subject title")]
                    match[inputDatabits.index("subject title")] = title[title.find("of",5)+3:]
                    subID = match[inputDatabits.index("subject id")]
                    match[inputDatabits.index("subject id")] = subID.replace(matchIsoform,match[inputDatabits.index("subject acc.")],1)
            
                filteredArray.append(match)
    except:
        print >> sys.stderr, ("Error: Missing Subject Accessifier (subject acc.), Subject Title (subject title), or Subject ID (subject id); skipping isoform filter")
        return array

    if len(filteredArray) > 1:
        for match in filteredArray:
            mIsoform = match[inputDatabits.index("subject acc.")]
            if '-' in mIsoform:
                match[inputDatabits.index("subject acc.")] = mIsoform[0:mIsoform.find('-')]
                titler = match[inputDatabits.index("subject title")]
                match[inputDatabits.index("subject title")] = titler[titler.find("of",5)+3:]
                subjID = match[inputDatabits.index("subject id")]
                match[inputDatabits.index("subject id")] = subjID.replace(mIsoform,match[inputDatabits.index("subject acc.")],1)
            filteredArrayFinal.append(match)



    if filteredArrayFinal:
        return filteredArrayFinal
    elif filteredArray:
        return filteredArray
    else:
        return array

# Apply Species filter
def applySpeciesFilter(arrayOfMatches, inputDatabits,spec):
    array = arrayOfMatches
    if not array:
        return array
    
    filteredArray=[]
    for match in array:
        speciesOfMatch = ''
        try:
            speciesOfMatch = match[inputDatabits.index("subject id")][match[inputDatabits.index("subject id")].find('_')+1:]
        except:
            print >> sys.stderr, ("Error: Missing Subject ID (subject id), skipping species filter")
            return array
        if speciesOfMatch == spec:
            filteredArray.append(match)
            
    return filteredArray

# Apply Annotation filter
def applyNameFilter(arrayOfMatches,inputDatabits):
    array = arrayOfMatches
    if not array:
        return array
    
    characterizedPOQArray = []
    characterizedABCArray = []
    uncharacterizedPOQArray = []
    uncharacterizedABCArray = []

    for match in array:
        try:
            nameOfMatch = match[inputDatabits.index("subject title")]
        except:
            print >> sys.stderr, ("Error: Missing Subject Title (subject title), skipping name filter")
            return array
        try:
            charOfSubjectAcc = match[inputDatabits.index("subject acc.")][0]
        except:
            print >> sys.stderr, ("Error: Missing Subject Accessifier (subject acc.), skipping the second name filter")
            return array

        characterizedFlag = True
        for string in ["uncharacterized","Uncharacterized","putative","Putative","hypothetical","Hypothetical"]:
            if string in nameOfMatch:
                characterizedFlag = False
        POQFlag = False
        for char in ['P','O','Q']:
            if charOfSubjectAcc == char:
                POQFlag = True

        if characterizedFlag and POQFlag:
            characterizedPOQArray.append(match)
        elif characterizedFlag and (not POQFlag):
            characterizedABCArray.append(match)
        elif (not characterizedFlag) and POQFlag:
            uncharacterizedPOQArray.append(match)
        else:
            uncharacterizedABCArray.append(match)
    
    
    if characterizedPOQArray:
        return characterizedPOQArray
    elif characterizedABCArray:
        return characterizedABCArray
    elif uncharacterizedPOQArray:
        return uncharacterizedPOQArray
    elif uncharacterizedABCArray:
        return uncharacterizedABCArray
    else:
        return array

# Print method for debug session
def printMethod(string, array):
    if userPreferences.get('debugMode'):
        print(string)
        if not array:
            print("0")
            print('\n')
            return
        for match in array:
            try:
                print match[inputDatabits.index("subject id")]
            except:
                print >> sys.stderr, ("Error: Missing Subject ID (subject id), skipping debug log")
        print("Length: "+str(len(array)))
        print('\n')
    else:
        pass

# Apply Species filter for top pick selection across all species
def applyOtherSpeciesFilter(arrayOfMatches, inputDatabits):
    array = arrayOfMatches
    if not array:
        return array
    arrayOfSpecies = userPreferences.get('species',[])

    arrayDict = {}

    for match in array:
        try:
            speciesOfMatch = match[inputDatabits.index("subject id")][match[inputDatabits.index("subject id")].find('_')+1:]
        except:
            print >> sys.stderr, ("Error: Missing Subject ID (subject id), skipping species filter")
            return array
        if len(speciesOfMatch) == 0:
            speciesOfMatch = 'Unspecified'
        if len(arrayOfSpecies) == 0:
            print("Error: No species preferences detected")
            return array

        for i in range(0,20):
            if len(arrayOfSpecies) > i and speciesOfMatch == arrayOfSpecies[i]:
                if arrayDict.has_key(i):
                    arrayDict[i].append(match)
                else:
                    arrayDict[i] = [match]

    for i in range(0,20):
        if arrayDict.get(i) and len(arrayDict[i]) != 0:
            return arrayDict[i]
    return array

# Select top pick across all species
def finalPickOverSpecies(array, inputDatabits):

    arrayOfMatches = array

    printMethod('For Top Pick: Length Of Initial Array', arrayOfMatches)

    if not userPreferences.get('speciesOverScore', False):
        arrayOfMatches = applyIdentityFilter(arrayOfMatches, inputDatabits,'hard')
        printMethod('For Top Pick: Length Of Array After Ranking Filter', arrayOfMatches)
    
    arrayOfMatches = applyOtherSpeciesFilter(arrayOfMatches, inputDatabits)
    printMethod('For Top Pick: Length Of Array After Species Filter', arrayOfMatches)
    
    if userPreferences.get('debugMode'):
        print('\n-------------------------------------------------------------------')
    return arrayOfMatches

# Select top pick for every species of every query
def bestMatchPerEntry(matchDictionary, key, inputDatabits, spec):
    
    arrayOfMatches = matchDictionary[key]
    
    if userPreferences.get('debugMode'):
        print("Starting algorithm on query: " + key[:key.find("~~^$%@^&@*`&<>")] + "\n") 
    printMethod('Length Of Initial Array', arrayOfMatches)

    arrayOfMatches = applyUserFilters(arrayOfMatches,inputDatabits)
    printMethod('1) Length Of Array After User Filter', arrayOfMatches)

    arrayOfMatches = applySpeciesFilter(arrayOfMatches, inputDatabits, spec)
    printMethod('2) Length Of Array After Species Filter', arrayOfMatches)

    arrayOfMatches = applyIdentityFilter(arrayOfMatches, inputDatabits,'soft')
    printMethod('3) Length Of Array After Identity Filter', arrayOfMatches)
    
    arrayOfMatches = applyBitscoreFilter(arrayOfMatches, inputDatabits,'soft')
    printMethod('4) Length Of Array After Bitscore Filter', arrayOfMatches)
            
    arrayOfMatches = applyEvalueFilter(arrayOfMatches, inputDatabits,'soft')
    printMethod('5) Length Of Array After E-Value Filter', arrayOfMatches)

    if not userPreferences.get('skipReviewed',False):
        arrayOfMatches = applyReviewFilter(arrayOfMatches, inputDatabits)
        printMethod('6) Length Of Array After Review Filter',arrayOfMatches)

    arrayOfMatches = applyIdentityFilter(arrayOfMatches, inputDatabits,'hard')
    printMethod('7) Length Of Array After Ranking Filter', arrayOfMatches)
    
    arrayOfMatches = applyNameFilter(arrayOfMatches,inputDatabits)
    printMethod('8) Length Of Array After Annotation Filter', arrayOfMatches)

    arrayOfMatches = applyIsoformFilter(arrayOfMatches, inputDatabits)
    printMethod('9) Length Of Array After Isoform Filter', arrayOfMatches)

    arrayOfMatches = applyEvalueFilter(arrayOfMatches, inputDatabits,'hard')
    printMethod('10) Length Of Array After Top-Pick Filter', arrayOfMatches)

    if userPreferences.get('debugMode'):
        print('\n-------------------------------------------------------------------')
    return arrayOfMatches

# Cluster similar queries to provide uniform decisions
def groupingFunction(gDict,keyO):
    groupDictionary = gDict
    assignDict = {}
    kOrder = keyO
    groupCounter = 1
    
    while len(groupDictionary) != 0:
        assignDict[groupCounter] = []
        pre_proteinIDs = groupDictionary[kOrder[0]]
        newIDs = pre_proteinIDs
        assignDict[groupCounter].append(kOrder[0])
        groupDictionary.pop(kOrder[0],None)
        kOrder.remove(kOrder[0])

        while len(newIDs) != 0:
            
            proteinIDs, rowID_List = search_matchrow(newIDs, pre_proteinIDs, groupDictionary)
            afterIDs = proteinIDs
            newIDs = list(set(afterIDs)-set(pre_proteinIDs))
            for value in rowID_List:
                assignDict[groupCounter].append(value)

            pre_proteinIDs = afterIDs

            for id in rowID_List:
                if groupDictionary.get(id,'No Value Found') != 'No Value Found':
                    groupDictionary.pop(id,None)
                    kOrder.remove(id)

        groupCounter += 1

    for key in assignDict:
        assignDict[key] = list(set(assignDict[key]))

    return assignDict

# Helper function for grouping analysis
def search_matchrow(nIDs, preproteinID, groupDict):
    groupDictionary = groupDict
    rowID_List = []
    newIDsFound = []
    
    for key in groupDictionary:
        for value in nIDs:
            if value in groupDictionary[key]:
                rowID_List.append(key)
                newIDsFound = newIDsFound + list(set(groupDictionary[key]) - set(newIDsFound))


    proteinID = newIDsFound + list(set(preproteinID) - set(newIDsFound))
    
    return proteinID,rowID_List

## THIS METHOD IS UNFINISHED CODE - DO NOT USE AS IS
## Provides uniform descisions for queries in clusters
#def clusterUniformity(finalDict, helperDict, assignedGroupDict, keyOrder, inputDatabits):
#
#    lengthOfClusterDict = {}
#    topPickDict = {}
#    print("START OF DEEP COPY")
#    fDict = copy.deepcopy(finalDict)
#    hDict = copy.deepcopy(helperDict)
#    print("END OF DEEP COPY")
#    
#    try:
#        a = inputDatabits.index("subject acc.")
#    except:
#        print >> sys.stderr, ("Error: Missing Subject Accessifier (subject acc.), skipping grouping analysis")
#        return finalDict,helperDict
#
#    for groupID in assignedGroupDict:
#        lengthOfClusterDict[groupID] = 0
#        for key in keyOrder:
#            if key in assignedGroupDict[groupID]:
#                lengthOfClusterDict[groupID] += 1
#    
#    for groupID in assignedGroupDict:
#        percentageDictionary = {}
#        for spec in userPreferences['species']:
#            for key in keyOrder:
#                if key in assignedGroupDict[groupID]:
#                    for value in finalDict[key][spec]:
#
#                        if percentageDictionary.has_key(value[inputDatabits.index("subject acc.")]):
#                            percentageDictionary[value[inputDatabits.index("subject acc.")]] += value[-1]
#                        else:
#                            percentageDictionary[value[inputDatabits.index("subject acc.")]] = value[-1]
#
#        for access in percentageDictionary:
#            percentageDictionary[access] = percentageDictionary[access] / float(lengthOfClusterDict[groupID])
##        print(percentageDictionary)
##        print(fDict)
#
#        best = []
#        topNumb = 0
#        for accesser in percentageDictionary:
#            if percentageDictionary[accesser] >= topNumb:
#                topNumb = percentageDictionary[accesser]
#        for accesser in percentageDictionary:
#            if percentageDictionary[accesser] == topNumb:
#                best.append(accesser)
#
##        print(best)
#
#        for spec in userPreferences['species']:
#            for key in keyOrder:
#                if key in assignedGroupDict[groupID]:
#                    
#                    for value in fDict[key][spec]:
#                        if percentageDictionary.has_key(value[inputDatabits.index("subject acc.")]) and value[inputDatabits.index("subject acc.")] in best:
#                            fDict[key][spec] = [value]
#                            hDict[key][spec] = [value[inputDatabits.index("subject acc.")]]
#
##        print("\n\n\n")
##        print(fDict)
#        speciesPercentageDictionary = {}
#
#
#        for spec in userPreferences['species']:
#            numbMatchesPerSpecies = 0
#            for key in keyOrder:
#                if key in assignedGroupDict[groupID]:
#                    for match in fDict[key][spec]:
#                        numbMatchesPerSpecies += 1
#                        if speciesPercentageDictionary.has_key(spec):
#                            speciesPercentageDictionary[spec] += match[-1]
#                        else:
#                            speciesPercentageDictionary[spec] = match[-1]
##            speciesPercentageDictionary[spec] = speciesPercentageDictionary[spec] / float(numbMatchesPerSpecies)
#
##        print(speciesPercentageDictionary)
#            if numbMatchesPerSpecies == 0:
#                speciesPercentageDictionary[spec] = 0
#            else:
#                speciesPercentageDictionary[spec] = speciesPercentageDictionary[spec] / float(numbMatchesPerSpecies)
#        print(groupID)
#        print(speciesPercentageDictionary)
#
#        maxSpec = 0
#        maxSpecArray= []
#        for spec in userPreferences['species']:
#            if speciesPercentageDictionary[spec] >= maxSpec:
#                maxSpec = speciesPercentageDictionary[spec]
#        for spec in userPreferences['species']:
#            if speciesPercentageDictionary[spec] == maxSpec:
#                maxSpecArray.append(spec)
#        print(maxSpecArray)
#        topPickArray = []
#        for spec in userPreferences['species']:
#            for key in keyOrder:
#                if key in assignedGroupDict[groupID] and spec in maxSpecArray:
#                    topPickArray += fDict[key][spec]
#
##        topPickArray = list(set(topPickArray))
#        print(topPickArray)
#
##        print("\n\nPrinting final Dict")
##        print(finalDict)
##        print("Finished printing final dict\n\n")
#
##        print(fDict)
#
#        bestTopPickSpeciesArray = applyOtherSpeciesFilter(topPickArray, inputDatabits)
#        for key in keyOrder:
#            if key in assignedGroupDict[groupID]:
##                finalDict[key]['Top Match'] = []
##                helperDict[key]['Top Match'] = []
#
#                #print(finalDict[key]['Top Match'])
#                accumulativeArray = []
#                for spec in userPreferences['species']:
#                    accumulativeArray += finalDict[key][spec]
#
#                for element in accumulativeArray:
#                    if element in bestTopPickSpeciesArray:
#                        print("yes")
#                        
#                        finalDict[key]['Top Match'] = [element]
#                        helperDict[key]['Top Match'] = []
#                        for match in finalDict[key]['Top Match']:
#                            helperDict[key]['Top Match'].append(match[inputDatabits.index("subject acc.")])
#                    else:
#                        print("no")
##                print(finalDict[key]['Top Match'])
##        print("\nFinish groupID\n")
#
##    print(fDict)
##    finalDict = fDict
##    helperDict = hDict
##    print(lengthOfClusterDict)
##    print(percentageDictionary)
#
#
#    return finalDict,helperDict

# Write all information to .csv file
def writeToFile(finalDict,helperDict, outputFile, keyOrder, inputDatabits,assignedGroupDict):

    try:
        with open(outputFile,'w') as f:
            
            for spec in reversed(userPreferences.get('species',[]) + ['Top Pick']):
                    inputDatabits.insert(1, spec)
            inputDatabits.insert(1,'Cluster ID')
            headerDesc = ",".join(inputDatabits)

            f.write(headerDesc)
            f.write('\r')
            for k in keyOrder:
                v = finalDict[k]
                if not v:
                    f.write(k[0:k.find("~~^$%@^&@*`&<>")] + ",No Matches Found From Givin Criterion\r\r")

                else:
                    row = []
                    queryID = k[0:k.find("~~^$%@^&@*`&<>")]
                    queryID = queryID.replace(',',';')
                    row.append(queryID)
                    
                    for ky, ve in assignedGroupDict.items():
                        if k in ve:
                            groupID = ky
                        if userPreferences.get('turnOffClustering',False):
                            groupID = 'Not Requested'
                    row.append(groupID)
                    
                    for spec in userPreferences['species']+['Top Match']:
                        ID = str(helperDict.get(k,'n/a').get(spec,'n/a'))
                        ID = ID.replace(',',';')
                        ID = ID.replace("'",'').replace("[",'').replace("]",'')
                            ###
                        percentID = ''
                        for match in finalDict.get(k,'n/a').get(spec,'n/a'):
                            percentID = percentID + str(match[-1]) + "%; "
                            
                        ID = ID + " (" + percentID +")"
                        if len(helperDict.get(k,'n/a').get(spec,'n/a')) == 0:
                            ID = 'n/a'
                        ID = ID.replace("; )",")")
                        ###
                        row.append(ID)
                    row = ','.join(map(str, row))
                    f.write(row + ',')
                    matchFound = False
                    if helperDict[k]['Top Match']:
                        matchFound = True
                        for element in range(len(v['Top Match'][0][:-1])):
                            if element > 0:
                                dataPiece = str(v['Top Match'][0][element])
                                dataPiece = dataPiece.replace(',',';')
                                f.write(dataPiece+',')
                        f.write('\r')
                    if matchFound == False:
                        f.write('\r')
            f.close()
        print("\nProgram End\n")
    except:
        raise RuntimeError("Error: Write to file. Please re-check all arguments and conditions and try again.")

# Usage/help page
def usage():
    print("")
    print("BLASTCustomOutputScript.py")
    print ("-" * 15)
    print("This script selects the top matches for each entry query from a text file generated by standalone/command-line blast.")
    print("Top matches are selected based on a a variety of filters, including species, percent identity and alignment length, bitscore, evalue, reviewed/unreviewed, isoforms, and naming conventions")
    print("")
    print("Usage: ")
    print("python BLASTCustomOutputScript.py -i results.txt -s <species>,...")
    print("")
    print("Arguments are as follows:")
    print("-i     Input blast-results file generated by earlier run of standalone blast")
    print("-o     (optional)Output file name, such as results.csv. Default: input file name appended with .csv and placed under the directory of your input file")
    print("-p     (optional)Percent Identity Cutoff, as a number 0-100. To be used to discard all matches below cutoff. Default: 75")
    print("-a     (optional)Query Coverage Cutoff, as a positive number. To be used to discard all matches below cutoff. Default: 75")
    print("-m     (optional)Mismatch tolerance, as positive integer. To be used to tolerate reviewed matches that have a certain number of mismatches/missing bits in alignment. Default: 0")
    print("-e     (optional)E-Value Cutoff, as a positive number. To be used to discard all matches greater than cutoff. Default: 10")
    print("-s     (optional)Species preference, in all caps. To be used to display top matches for each species listed. Example format: HUMAN,RAT,MOUSE. Default: Display top match regardless of species")
    print("-r     (optional)Flag that specifies to skip reviewed/unreviewed filter. To be used if one does not wish to distinguish between reviewed and unreviewed matches.")
    print("-c     (optional)Flag that makes Species more important than Scoring in the final pick. To be used if one wishes to recieve matches of preferred species over those of less preferred species but that more closely resemble the sample.")
    print("-t     (optional)CLUSTERING UNIFORMITY NOT YET IMPLEMENTED, THIS ARG IS NOT CURRENTLY RELEVENT. Turn off clustering analysis. To be used if one wishes each query to be evaluated independently rather than in similarity clusters.")
    print("-d     (optional)Debug mode: turns on wrapper debugging to stdout")
    print("-h     (optional)Return to this page.")
    print("")


def exception_handler(exception_type, exception, traceback):
    print "%s: %s" % (exception_type.__name__, exception)


userPreferences = {}

# Start program here
if __name__ == "__main__":
    
    sys.excepthook = exception_handler
    inputFile = ''
    outputFile = ''

    pidentCutoff = 0
    alengthCutoff = 0
    mismatchTolerance = 0
    evalueCutoff = float("inf")
    skipReviewedFilter = False
    speciesOverScore = False
    species = []
    debugMode = False
    turnOffClustering = False
    
    
    if len(sys.argv)<=1:
        usage()
        exit(0)
    else:
        opts, _ = getopt.getopt(sys.argv[1:], "i:o:p:a:m:e:s:rctdh")

        for opt,arg in opts:
            
            if opt == "-i" :
                inputFile = arg
            elif opt == "-o" :
                outputFile = arg
            elif opt == "-p" :
                pidentCutoff = float(arg)
                userPreferences['pident'] = pidentCutoff
            elif opt == "-a" :
                alengthCutoff = float(arg)
                userPreferences['alength'] = alengthCutoff
            elif opt == "-m" :
                mismatchTolerance = int(arg)
                userPreferences['mismatch'] = mismatchTolerance
            elif opt == "-e":
                evalueCutoff = float(arg)
                userPreferences['evalue'] = evalueCutoff
            elif opt == "-s" :
                species = arg.split(",")
                userPreferences['species'] = species
            elif opt == "-r":
                skipReviewedFilter = True
                userPreferences['skipReviewed'] = skipReviewedFilter
            elif opt == "-c":
                speciesOverScore = True
                userPreferences['speciesOverScore'] = speciesOverScore
            elif opt == "-t":
                turnOffClustering = True
                userPreferences['turnOffClustering'] = turnOffClustering
            elif opt == "-d":
                debugMode = True
                userPreferences['debugMode'] = debugMode
            elif opt == "-h":
                usage()
                exit(0)
            else:
                print >> sys.stderr, ("Warning! Command-line argument: %s not recognized. Exiting...")
                exit(1)

    try:
        assert isinstance(inputFile,str) and len(inputFile) > 0
    except:
        print >> sys.stderr, "Error: No input (-i) file detected"
        exit(1)

    try:
        if pidentCutoff:
            assert isinstance(pidentCutoff,float) and pidentCutoff >= 0 and pidentCutoff <= 100
    except:
        print >> sys.stderr, "Error: Percent Identity (-p) argument does not match its parameters"
        exit(1)
        
    try:
        if alengthCutoff:
            assert isinstance(alengthCutoff,float) and alengthCutoff >= 0 
    except:
        print >> sys.stderr, "Error: Alignment Length (-a) argument does not match its parameters"
        exit(1)

    try:
        if mismatchTolerance:
            assert isinstance(mismatchTolerance,int) and mismatchTolerance >= 0
    except:
        print >> sys.stderr, "Error: Mismatch Tolerance (-m) argument does not match its parameters"
        exit(1)
    
    try:
        if evalueCutoff:
            assert isinstance(evalueCutoff,float) and evalueCutoff >= 0 
    except:
        print >> sys.stderr, "Error: E-Value (-e) argument does not match its parameters"
        exit(1)

    if len(species) == 0:
        print >> sys.stderr, ("Error: Species (-s) argument must be specified")
        exit(1)
    try:
        if species:
            assert isinstance(species,list) and len(species) >= 1 and len(species) <= 8
    except:
        print >> sys.stderr, ("Error: Species (-s) argument does not match its parameters")
        exit(1)

    if userPreferences.get('debugMode'):
        sys.stdout.write("User Arguments: ")
        print(userPreferences)

    # Begin algorithm here

    print("\nBLAST_BestHit - Program Start...Version 2.5\n")

    matchDict, keyOrder, inputDatabits = readTabular(inputFile)
    finalDict = matchDict.fromkeys(matchDict,[])
    helperDict = matchDict.fromkeys(matchDict,[])
    groupDict = matchDict.fromkeys(matchDict,[])
    for key in finalDict:
        finalDict[key] = {}
        helperDict[key] = {}
        groupDict[key] = []
        topPickArray = []

        for spec in userPreferences['species']:
            finalDict[key][spec] = bestMatchPerEntry(matchDict, key, inputDatabits, spec)
            helperDict[key][spec] = []
            for match in finalDict[key][spec]:
                helperDict[key][spec].append(match[inputDatabits.index("subject acc.")])
            topPickArray += finalDict[key][spec]
            if helperDict[key][spec]:
                groupDict[key].append(helperDict[key][spec])

        finalDict[key]['Top Match'] = finalPickOverSpecies(topPickArray, inputDatabits)
        helperDict[key]['Top Match'] = []
        for match in finalDict[key]['Top Match']:
            helperDict[key]['Top Match'].append(match[inputDatabits.index("subject acc.")])

    for key in groupDict:
        groupDict[key] = list(itertools.chain.from_iterable(groupDict[key]))

    keyOrderBackup = list(keyOrder)
    assignedGroupDict = {}  #Form: {1: [PALMD_PIG,FI_PIG],2: [TREX_HUMAN,Y_PIG,KRO_MOUSE]}
    assignedGroupDict = groupingFunction(groupDict,keyOrder)

    #CLUSTERING ANALYSIS TO BE ADDED HERE: CURRENTLY THE FOLLOWING 5 LINES SHOULD NOT IMPLEMENTED
    if not turnOffClustering and False:
        try:
            finalDict,helperDict = clusterUniformity(finalDict, helperDict, assignedGroupDict, keyOrderBackup, inputDatabits)
        except:
            print >> sys.stderr, ("Error: Handeling cluster uniformity. Skipping this step.")

    if len(outputFile) == 0:
        outputFile = inputFile[:-3] + "csv"

    writeToFile(finalDict, helperDict, outputFile, keyOrderBackup, inputDatabits,assignedGroupDict)

