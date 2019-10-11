import numpy as np
import pandas as pd
import math, sys, getopt, gzip, random
from sklearn import svm
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from imblearn.over_sampling import SMOTENC

def read(path):
    reval = []
    header = []
    variantMap = []
    if path[-3:] == ".gz":
        with gzip.open(path, mode='rt') as f:
            variantIndex = 0
            for line in f:
                if len(line.strip()) == 0:
                    continue
                if (line[:14] == "##INFO=<ID=CSQ"):
                    info = line[line.find("Format: ") + 8:-3]
                    header = info.strip().split('|')
                    continue
                elif (line[:2] == "##"):
                    continue
                elif (line[:1] == "#" ):
                    #manage header here
                    #header = line.strip().split('\t')
                    continue
                s = line.strip().split('\t')
                #Keeping track of which lines correspond to which variants in order. Some annotations have multiple entries so we dont want to mix em up.
                variantMap.append(variantIndex)
                INFOField = s[7].split(';')
                for field in INFOField:
                    if field[:4] == "CSQ=":
                        temp = field[4:].split(',')
                        for annofield in temp:
                            reval.append(np.array(annofield.split('|')))
                            #print(annofield.split('|'))
                    variantIndex += 1

    else:
        with open(path, encoding='utf8') as f:
            for line in f:
                if len(line.strip()) == 0:
                    continue
                if (line[:14] == "##INFO=<ID=CSQ"):
                    info = line[line.find("Format: ") + 8:-3]
                    header = info.strip().split('|')
                    continue
                elif (line[:2] == "##"):
                    continue
                elif (line[:1] == "#" ):
                    #manage header here
                    #header = line.strip().split('\t')
                    continue
                s = line.strip().split('\t')
                #vector = [v for v in s]
                INFOField = s[7].split(';')
                #VEPField = []
                for field in INFOField:
                    if field[:4] == "CSQ=":
                        temp = field[4:].split(',')
                        for annofield in temp:
                            reval.append(np.array(annofield.split('|')))
                            #print(annofield.split('|'))
                        
                #reval.append(VEPField)
    return np.array(header), np.array(reval), np.array(variantMap)

def readTrainingLabels(path):
    labels = []
    if path[-3:] == ".gz":
        with gzip.open(path, mode='rt') as f:
            for line in f:
                if len(line.strip()) == 0:
                    continue
                labels.append(int(line))
    else:
        with open(path, encoding='utf8') as f:
            for line in f:
                if len(line.strip()) == 0:
                    continue
                labels.append(int(line))

    return labels
                

def generateLabels(impact):
    result = []
    for label in impact:
        if (label == "HIGH"):
            result.append(1)
        if (label == "MODERATE"):
            temp = random.random()
            if (temp < .5):
                result.append(1)
            else:
                result.append(0)
        if (label == "LOW"):
            temp = random.random()
            if (temp < .1):
                result.append(1)
            else:
                result.append(0)
        if (label == "MODIFIER"):
            result.append(0)
    return result


def printWeights(weights, outFile):
    weights.tofile(outFile, sep="\n", format="%2.6f")
    #with open(outFile, mode="w+", encoding='utf8') as f:
        #f.write(weights)
        #for w in range(weights.size()):
            #f.write(w, weights[w,:])
    
def splitConsequenceColumn(data):
    parsedConsequnce = []
    names = set()
    categories = ['3_prime_UTR_variant', '5_prime_UTR_variant', 'NMD_transcript_variant',
 'coding_sequence_variant', 'downstream_gene_variant', 'feature_elongation',
 'feature_truncation', 'frameshift_variant', 'inframe_deletion',
 'inframe_insertion', 'intron_variant', 'missense_variant',
 'non_coding_transcript_exon_variant', 'non_coding_transcript_variant',
 'splice_acceptor_variant', 'splice_donor_variant', 'splice_region_variant',
 'start_lost', 'stop_gained', 'stop_lost', 'synonymous_variant',
 'upstream_gene_variant']
    for i in data:
        variantConsequence = np.zeros(len(categories), dtype=np.int8)
        elements = i.split('&')
        index = 0
        for j in elements:
            names.add(j)
            for k in categories[index:]:
                if(j == k):
                    variantConsequence[index] = 1
                index += 1
        parsedConsequnce.append(variantConsequence)
    return parsedConsequnce, names


def main(argv):
    learningMode = False
    w = None
    opts, args = getopt.getopt(argv, "hlwpb:o:d:", ['--labels', "--weight", "--output", "--data", "--back"])
    for opt, arg in opts:
        if opt == '-h':
            print("I need help too.")
            sys.exit(1)
        elif opt in ('-l', '--labels'):
            learningMode = True
            labelsPath = arg
        elif opt in ('-w', '--weight'):
            w = arg
        elif opt in ('-o', '--output'):
            outFile = arg
        elif opt in ('-d', '--data'):
            annotationFile = arg
        elif opt in ('-p', '--pheno'):
            phenoFile = arg
        elif opt in ('-b', '--back'):
            backFile = arg
            print(arg)


    header, rawData, userVariantMap = read(annotationFile)
    print(header)
    header, rawBackData, backVariantMap = read(backFile)
    if learningMode:
        Y = readTrainingLabels(labelsPath)


    print(rawData.shape)
    tempUserConsq, categories = splitConsequenceColumn(rawData[:, 1])
    userConsequenceData = pd.DataFrame(data=tempUserConsq, columns=categories)

    tempBackConsq, categories = splitConsequenceColumn(rawBackData[:, 1])
    print(categories, len(categories))
    backConsequenceData = pd.DataFrame(data=tempBackConsq, columns=categories)

    data = pd.DataFrame(data=rawData, columns=header)
    data = pd.concat([data, userConsequenceData], axis=1)

    backData = pd.DataFrame(data=rawBackData, columns=header)
    backData = pd.concat([backData, backConsequenceData], axis=1)
    print(backData["IMPACT"])
    Y = generateLabels(backData["IMPACT"])
    backData = pd.get_dummies(data=backData, columns=['IMPACT'], prefix=['IMPACT'])
    
    #Y = pd.DataFrame(data=backData.loc[:, 'IMPACT_HIGH'], columns=['IMPACT_HIGH'])
    #print("This is Y:", Y)
    #print(np.sum(Y))  
    backData = backData.drop(columns=['IMPACT_HIGH', 'Consequence'])
    backData = backData.drop(columns=['Allele', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'DISTANCE', 'STRAND', 'FLAGS', 'SYMBOL_SOURCE', 'HGNC_ID'])
    
    
    
    smote_nc = SMOTENC(categorical_features=np.arange(0,24), random_state=0)
    X_resampled, y_resampled = smote_nc.fit_resample(backData, Y)

    print("Array Shapes", np.shape(X_resampled), np.shape(y_resampled))
    
    print("X_resampled, y_resampled", X_resampled, y_resampled)
    
    #model = LogisticRegression(random_state=0, solver='liblinear', fit_intercept=True, C=10.0).fit(backConsequenceData, np.ravel(Y), sample_weight=w)
    model = LogisticRegression(random_state=0, solver='liblinear', fit_intercept=True, C=10.0).fit(X_resampled, np.ravel(y_resampled), sample_weight=w)
    #data = pd.get_dummies(data=data, columns=['IMPACT'],  prefix=['IMPACT'])
    #data = data.drop(columns=['IMPACT', 'Consequence'])
    #data = data.drop(columns=['Allele', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'DISTANCE', 'STRAND', 'FLAGS', 'SYMBOL_SOURCE', 'HGNC_ID'])
    #print(list(data.columns.values))

    

    causitiveVariants = model.predict(userConsequenceData.iloc[:,0:23])
    print(causitiveVariants, np.sum(causitiveVariants))
    scores = model.predict_proba(userConsequenceData.iloc[:,0:23])
    print(scores)
    index = 0
    nonzeroColumns = []
    for var in causitiveVariants:
        if var == 1:
            nonzeroColumns.append(index)
            print(index, scores[index,:])
        index += 1

    print(userConsequenceData.iloc[nonzeroColumns, [7,14,15,17,18,19]])

    #print(model.predict_log_proba(userConsequenceData))

    printWeights(scores[:,1], outFile)

    #model.get_params() #If i want the final weights 

    #printWeights(np.concatenate(model.intercept_, model.coef_), outFile)






if __name__ == '__main__':
    main(sys.argv[1:])  
