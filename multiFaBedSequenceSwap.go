package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

// SubsequenceSwap swaps a series of specified regions from a bed file between two sequences in a MultiFa object
// The sequence called "background" is taken, and parts of "foreground" are swapped with the corresponding parts in background
// input of this function includes the multiFa file name, the index of the background and foreground sequences from the multiFa, and a bed object with the relevant
// start and stop indices of the swap regions
// it will return a fasta slice (which includes several sequences, akin to a multifa), with one more sequence that the sequences in the original
// multiFa; this extra sequence is the background sequence with relevant foreground spots swapped as dictated by the bed file.
func SubsequenceSwap(inFile string, background, foreground int, swapRegions []bed.Bed) []fasta.Fasta {
	var currStart, currEnd, currPos int //these variables will be used to swap regions later in the code

	// Load the original sequences from the multiFa file.
	originalSeqs := fasta.Read(inFile) //Converting a multiFa to a Fasta file

	// Ensure the specified background and foregound sequence indices are within the valid range.
	if background < 0 || background >= len(originalSeqs) || foreground < 0 || foreground >= len(originalSeqs) {
		log.Fatalf("Error: Invalid sequence index.\n")
		return []fasta.Fasta{}
	}
	answerSeq := fasta.Copy(originalSeqs[background]) //copy of the background sequence is made
	for _, swapRegion := range swapRegions {          //loop over every single swap region in the bed
		currStart = swapRegion.ChromStart
		currEnd = swapRegion.ChromEnd

		// Validate the swap region.
		if currStart < 0 || currEnd > len(answerSeq.Seq) || currStart >= currEnd {
			log.Fatalf("Error: Invalid swap region. \n")
			return []fasta.Fasta{}
		}

		// Swap the corresponding portion of answerSeq with the foreground sequence.
		for currPos = currStart; currPos < currEnd; currPos++ {
			answerSeq.Seq[currPos] = originalSeqs[foreground].Seq[currPos]
		}
	}
	var finalSeqs []fasta.Fasta
	finalSeqs = append(originalSeqs, answerSeq)
	// Create and return the new Fasta with one more sequence than the input.
	return finalSeqs
}

func main() {
	// Example usage of SubsequenceSwap function

	var inFile string = "Documents/GitHub/gonomics/cmd/multFaVisualizer/testdata/test.fa"

	// Example swap region
	swapRegion := []bed.Bed{
		{ChromStart: 1,
			ChromEnd: 3},
	}

	// Specify the background and foreground sequences here: this is an example usage
	backgroundSeqIndex := 0
	foregroundSeqIndex := 1

	result := SubsequenceSwap(inFile, backgroundSeqIndex, foregroundSeqIndex, swapRegion)

	fmt.Println(result)
}
