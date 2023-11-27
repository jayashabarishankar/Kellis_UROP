from Bio import SeqIO
import numpy as np
import argparse
from scipy.linalg import fractional_matrix_power




def prob_matrix(filename, outFile, additive_smoothing=0.1, divergence=True, rescale=None):
    sequences= [i for i in SeqIO.parse(filename, "fasta")] #list of SeqIO objects
    assert len(sequences)==2 #2 sequences will be compared
    assert len(sequences[0]) == len(sequences[1]) #sequences must be same length

    #ACGT order will be followed, sequence 1 is the one on top/ corresponding to columns
    count_matrix = np.zeros(shape=(4, 4))
    prob_matrix = np.zeros(shape=(4, 4))

    for currPos in range(0, len(sequences[0])):
        if validBase(sequences[0][currPos]) and validBase(sequences[1][currPos]):
            count_matrix[baseToInt(sequences[0][currPos])][baseToInt(sequences[1][currPos])] += 1
    for currRow in range(0, 4):
        currRowSum = 0
        for currCol in range(0, 4):
            currRowSum += count_matrix[currRow][currCol]
        for currCol in range(0, 4):
            prob_matrix[currRow][currCol] = count_matrix[currRow][currCol] / currRowSum


    count_matrix += additive_smoothing

    #Normalize rows to ensure the sum of each row is 1
    row_sums = count_matrix.sum(axis=1, keepdims=True)
    prob_matrix = count_matrix / row_sums

    # If rescale is provided, rescale the mutation spectrum
    if rescale!=None:
        rescaled_matrix = rescale_spectrum(prob_matrix, rescale)
        prob_matrix=rescaled_matrix

    print(count_matrix)
    print(prob_matrix)

    if divergence:
        divergence_value = calculate_divergence(prob_matrix)
        print("Divergence:", divergence_value)


def validBase(base):
    if base == "A":
        return True
    if base == "C":
        return True
    if base == "G":
        return True
    if base == "T":
        return True
    return False

def baseToInt(base):
        if base == "A":
              return 0
        if base == "C":
              return 1
        if base == "G":
              return 2
        if base == "T":
              return 3
        return -1


def calculate_divergence(prob_matrix):
    # Calculate divergence as the sum of off-diagonal elements divided by the sum of all elements, the number varies between 0 and 1
    #  how much two species have "diverged" or changed from each other..
    divergence_value = np.sum(prob_matrix) - np.trace(prob_matrix)
    divergence_value = divergence_value/ np.sum(prob_matrix)
    return divergence_value

def rescale_spectrum(prob_matrix, t):
    #current_divergence = calculate_divergence(prob_matrix)
    #rescale = target_divergence / current_divergence
    #rescaled_matrix = prob_matrix * rescale
    rescaled_matrix = fractional_matrix_power(prob_matrix, t) #t is time

    return rescaled_matrix


def main(inFa, outFile, additive_smoothing, divergence, rescale):


    prob_matrix(inFa, outFile)

if __name__ == "__main__":
    # Set default values
    default_inFa = "test.fa"
    default_outFile = "outFile.py"

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Make empirical mutation spectrum from multiFa.")
    parser.add_argument("inFa", nargs="?", default=default_inFa, help="A multiFa alignment file as input.")
    parser.add_argument("outFile", nargs="?", default=default_outFile, help="The output mutation spectrum.")
    parser.add_argument("--additive_smoothing", type=float, help="ensures probability is not zero")
    parser.add_argument("--divergence", action="store_true", help="Calculate empirical divergence from mutation spectrum")
    parser.add_argument("--rescale", type=float, help="Rescale the mutation spectrum based on user-specified divergence")
    args = parser.parse_args()
    # Call the main function with the provided or default values
    main(args.inFa, args.outFile, args.additive_smoothing, args.divergence, args.rescale)


"""
if __name__ == "__main__":
    inFa="test.fa"
    outFile="outFile.py"
    parser = argparse.ArgumentParser(description="Make empirical mutation spectrum from multiFa.")
    parser.add_argument("inFa", help="A multiFa alignment file as input.")
    parser.add_argument("outFile", help="The output mutation spectrum.")
    args = parser.parse_args()
    main(args.inFa, args.outFile)
"""
