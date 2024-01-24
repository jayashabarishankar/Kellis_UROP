package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"strconv"
	"testing"
)

var SubsequenceSwapTestCases = []struct {
	InFile          string
	Background      int
	Foreground      int
	SwapRegions     []bed.Bed
	ExpectedNumSeqs int //expected number of sequences in the output fasta
}{
	{
		InFile:          "https://github.com/vertgenlab/gonomics/blob/main/cmd/multFaVisualizer/testdata/test.fa",
		Background:      0,
		Foreground:      1,
		SwapRegions:     []bed.Bed{{ChromStart: 1, ChromEnd: 3}, {ChromStart: 5, ChromEnd: 8}},
		ExpectedNumSeqs: 7,
	},
}

func TestSubsequenceSwap(t *testing.T) {
	for _, tc := range SubsequenceSwapTestCases {
		//loop through all test cases
		result := SubsequenceSwap(tc.InFile, tc.Background, tc.Foreground, tc.SwapRegions)

		// Check the number of sequences in the output, make sure it is one more than input
		if len(result) != tc.ExpectedNumSeqs {
			t.Errorf("Expected" + strconv.Itoa(tc.ExpectedNumSeqs) + "sequences but got" + strconv.Itoa(len(result)))
		}

		// Check the correctness of the modified sequence
		resultSeq := result[len(result)-1] //last sequence is the resultSeq

		// Verify that the specified regions are correctly swapped; we are looping through every swapped region
		//in the result and then checking if it matches the corresponding foreground part
		for i, swapRegion := range tc.SwapRegions {
			resultSlice := resultSeq.Seq[swapRegion.ChromStart:swapRegion.ChromEnd]
			expectSlice := result[tc.Foreground].Seq[swapRegion.ChromStart:swapRegion.ChromEnd]

			// compare foreground sequence and result sequence to make sure they match
			if dna.CompareSeqsIgnoreCase(expectSlice, resultSlice) != 0 {
				t.Errorf("Incorrect swap in specified region " + strconv.Itoa(i))
			}
		}
	}
}
